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
  bwa mem -t ${task.cpus} ${ref_fa} ${r1} ${r2} |
    samtools sort -@ ${task.cpus} -o ${sample_id}.sorted.bam
  samtools index ${sample_id}.sorted.bam
  """
}

workflow {
  // stage reference + indices into each task workdir
  def ref_fa  = file(params.ref_fasta)
  def ref_idx = file("${params.ref_fasta}.*")

  reads_ch = Channel
    .fromPath(params.samplesheet)
    .splitCsv(header:true)
    .map { row -> tuple(row.sample, file(row.fastq_1), file(row.fastq_2)) }

  trimmed_ch = FASTP(reads_ch)

  map_input_ch = trimmed_ch.map { sid, r1, r2 -> tuple(sid, r1, r2, ref_fa, ref_idx) }

  BWA_MEM(map_input_ch)
}
