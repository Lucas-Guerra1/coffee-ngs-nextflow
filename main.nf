nextflow.enable.dsl=2

process FASTP {
  tag "$sample_id"
  publishDir "${params.outdir}/01_fastp", mode: 'copy'
  conda "envs/bioinfo.yml"

  input:
    tuple val(sample_id), path(r1), path(r2)

  output:
    tuple val(sample_id), path("${sample_id}_R1.trim.fq.gz"), path("${sample_id}_R2.trim.fq.gz")

  script:
  """
  fastp -i ${r1} -I ${r2} \
        -o ${sample_id}_R1.trim.fq.gz -O ${sample_id}_R2.trim.fq.gz \
        --detect_adapter_for_pe \
        --thread ${task.cpus}
  """
}

process BWA_MEM {
  tag "$sample_id"
  publishDir "${params.outdir}/02_bam", mode: 'copy'
  conda "envs/bioinfo.yml"

  input:
    tuple val(sample_id), path(r1), path(r2), path(ref_fa), path(ref_idx)

  output:
    tuple val(sample_id), path("${sample_id}.sorted.bam"), path("${sample_id}.sorted.bam.bai")

  script:
  """
  bwa mem -t ${task.cpus} \
    -R "@RG\\tID:${sample_id}\\tSM:${sample_id}\\tPL:ILLUMINA" \
    ${ref_fa} ${r1} ${r2} |
    samtools sort -@ ${task.cpus} -o ${sample_id}.sorted.bam

  samtools index ${sample_id}.sorted.bam
  """
}

process BCFTOOLS_CALL_SAMPLE {
  tag "$sample_id"
  publishDir "${params.outdir}/03_vcf", mode: 'copy'
  conda "envs/bioinfo.yml"

  input:
    tuple val(sample_id), path(bam), path(bai), path(ref_fa), path(ref_idx)

  output:
    tuple val(sample_id), path("${sample_id}.vcf.gz"), path("${sample_id}.vcf.gz.csi")

  script:
  """
  bcftools mpileup -f ${ref_fa} -Ou -s ${sample_id} ${bam} |
    bcftools call -mv -Oz -o ${sample_id}.vcf.gz

  bcftools index -c ${sample_id}.vcf.gz
  """
}

process VCF_FILTER_SNPS {
  tag "$sample_id"
  publishDir "${params.outdir}/04_snps", mode: 'copy'
  conda "envs/bioinfo.yml"

  input:
    tuple val(sample_id), path(vcf_gz), path(vcf_csi)

  output:
    tuple val(sample_id), path("${sample_id}.snps.vcf.gz"), path("${sample_id}.snps.vcf.gz.csi")

  script:
  """
  bcftools view -m2 -M2 -v snps -Oz -o ${sample_id}.snps.vcf.gz ${vcf_gz}
  bcftools index -c ${sample_id}.snps.vcf.gz
  """
}

process CONSENSUS_FASTA {
  tag "$sample_id"
  publishDir "${params.outdir}/05_fasta", mode: 'copy'
  conda "envs/bioinfo.yml"

  input:
    tuple val(sample_id), path(snp_vcf_gz), path(snp_vcf_csi), path(ref_fa)

  output:
    tuple val(sample_id), path("${sample_id}.consensus.fa")

  script:
  """
  bcftools consensus -f ${ref_fa} ${snp_vcf_gz} > ${sample_id}.consensus.fa
  """
}

process FASTTREE {
  tag "tree"
  publishDir "${params.outdir}/06_tree", mode: 'copy'
  conda "envs/bioinfo.yml"

  input:
    path fasta_files

  output:
    path("tree.nwk")

  script:
  """
  # With 1 sample the tree is trivial; with multiple samples this becomes meaningful
  FastTree -nt ${fasta_files} > tree.nwk
  """
}

workflow {
  def ref_fa  = file(params.ref_fasta)
  def ref_idx = file("${params.ref_fasta}.*")

  reads_ch = Channel
    .fromPath(params.samplesheet)
    .splitCsv(header:true)
    .map { row -> tuple(row.sample, file(row.fastq_1), file(row.fastq_2)) }

  trimmed_ch = FASTP(reads_ch)

  map_input_ch = trimmed_ch.map { sid, r1, r2 -> tuple(sid, r1, r2, ref_fa, ref_idx) }

  bam_ch = BWA_MEM(map_input_ch)

  call_in_ch = bam_ch.map { sid, bam, bai -> tuple(sid, bam, bai, ref_fa, ref_idx) }

  vcf_ch = BCFTOOLS_CALL_SAMPLE(call_in_ch)

  snps_ch = VCF_FILTER_SNPS(vcf_ch)

  consensus_ch = CONSENSUS_FASTA(snps_ch.map { sid, vcf, csi -> tuple(sid, vcf, csi, ref_fa) })

  // collect FASTAs across samples for a single tree
  fasta_list = consensus_ch.map { sid, fa -> fa }.collect()

  FASTTREE(fasta_list)
}
