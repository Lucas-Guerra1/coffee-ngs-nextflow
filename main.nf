nextflow.enable.dsl=2

/*
 * PARÂMETROS UNIVERSAIS
 */
params.reads     = "data/*_{1,2}.fastq.gz"
params.ref_fasta = "data/ref.fa"
params.outdir    = "results"

// Validação de parâmetros
if (!params.reads) error "Parâmetro --reads não especificado"
if (!params.ref_fasta) error "Parâmetro --ref_fasta não especificado"

process BWA_INDEX {
    tag "$ref_fa"
    publishDir "${params.outdir}/00_ref_index", mode: 'copy'
    conda "envs/bioinfo.yml"

    input:
    path ref_fa

    output:
    tuple path(ref_fa), path("${ref_fa}*"), emit: reference_bundle

    script:
    """
    bwa index ${ref_fa}
    """
}

process FASTP {
    tag "$sample_id"
    publishDir "${params.outdir}/01_fastp", mode: 'copy'
    conda "envs/bioinfo.yml"

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_R1.trim.fq.gz"), path("${sample_id}_R2.trim.fq.gz"), emit: trimmed_reads

    script:
    """
    fastp -i ${reads[0]} -I ${reads[1]} \
          -o ${sample_id}_R1.trim.fq.gz -O ${sample_id}_R2.trim.fq.gz \
          --detect_adapter_for_pe --thread ${task.cpus}
    """
}

process BWA_MEM {
    tag "$sample_id"
    publishDir "${params.outdir}/02_bam", mode: 'copy'
    conda "envs/bioinfo.yml"

    input:
    tuple val(sample_id), path(r1), path(r2), path(ref_fa), path(indices)

    output:
    tuple val(sample_id), path("${sample_id}.sorted.bam"), path("${sample_id}.sorted.bam.bai"), emit: bam_files

    script:
    """
    bwa mem -t ${task.cpus} \
      -R "@RG\\tID:${sample_id}\\tSM:${sample_id}\\tPL:ILLUMINA" \
      ${ref_fa} ${r1} ${r2} | \
      samtools sort -@ ${task.cpus} -o ${sample_id}.sorted.bam

    samtools index ${sample_id}.sorted.bam
    """
}

process BCFTOOLS_CALL_SAMPLE {
    tag "$sample_id"
    publishDir "${params.outdir}/03_vcf", mode: 'copy'
    conda "envs/bioinfo.yml"

    input:
    tuple val(sample_id), path(bam), path(bai), path(ref_fa), path(indices)

    output:
    tuple val(sample_id), path("${sample_id}.vcf.gz"), path("${sample_id}.vcf.gz.csi"), emit: vcf_files

    script:
    """
    bcftools mpileup -f ${ref_fa} -Ou -s ${sample_id} ${bam} | \
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
    tuple val(sample_id), path("${sample_id}.snps.vcf.gz"), path("${sample_id}.snps.vcf.gz.csi"), emit: snps_files

    script:
    """
    bcftools view -m2 -M2 -v snps -Oz -o ${sample_id}.snps.vcf.gz ${vcf_gz}
    bcftools index -c ${sample_id}.snps.vcf.gz
    """
}

process MERGE_VCFS {
    tag "merge_all_samples"
    publishDir "${params.outdir}/05_merged_vcf", mode: 'copy'
    conda "envs/bioinfo.yml"

    input:
    path vcf_files
    path csi_files

    output:
    tuple path("all_samples.vcf.gz"), path("all_samples.vcf.gz.csi"), emit: merged_vcf

    script:
    """
    bcftools merge -m none -Oz -o all_samples.vcf.gz ${vcf_files}
    bcftools index -c all_samples.vcf.gz
    """
}

process VCF_TO_FASTA {
    tag "vcf_to_fasta"
    publishDir "${params.outdir}/06_alignment", mode: 'copy'
    conda "envs/bioinfo.yml"

    input:
    tuple path(merged_vcf), path(merged_csi)

    output:
    path("snp_alignment.fasta"), emit: alignment

    script:
    """
    # Converte VCF para formato FASTA contendo apenas as posições SNP
    bcftools query -f '%CHROM\\t%POS\\t%REF[\\t%TGT]\\n' ${merged_vcf} | \
    awk 'BEGIN {
        header_printed = 0
    }
    NR == 1 {
        # Primeira linha - imprime cabeçalhos das amostras
        for (i = 4; i <= NF; i++) {
            samples[i-3] = "sample_" (i-3)
            seq[i-3] = ""
        }
    }
    {
        # Para cada linha (SNP), adiciona o genótipo de cada amostra
        for (i = 4; i <= NF; i++) {
            genotype = \$i
            # Remove fases (/) ou (|) e pega apenas o primeiro alelo
            gsub(/[\\/|]/, "", genotype)
            if (length(genotype) > 0) {
                seq[i-3] = seq[i-3] substr(genotype, 1, 1)
            } else {
                seq[i-3] = seq[i-3] "N"
            }
        }
    }
    END {
        # Imprime FASTA
        for (i = 1; i <= length(seq); i++) {
            print ">" samples[i]
            print seq[i]
        }
    }' > snp_alignment.fasta
    
    # Pega os nomes reais das amostras do VCF
    bcftools query -l ${merged_vcf} > sample_names.txt
    
    # Substitui os nomes genéricos pelos nomes reais
    i=1
    while read sample; do
        sed -i "s/>sample_\$i/>\$sample/" snp_alignment.fasta
        i=\$((i+1))
    done < sample_names.txt
    """
}

process FASTTREE {
    tag "tree"
    publishDir "${params.outdir}/07_tree", mode: 'copy'
    conda "envs/bioinfo.yml"

    input:
    path alignment

    output:
    path("tree.nwk"), emit: tree

    script:
    """
    # Agora com apenas SNPs, FastTree será muito mais rápido
    FastTree -nt -gtr ${alignment} > tree.nwk
    """
}

workflow {
    // 1. Setup da Referência
    ref_ch = Channel.fromPath(params.ref_fasta)
    reference_bundle = BWA_INDEX(ref_ch)
    
    // 2. Setup das Reads
    reads_ch = Channel.fromFilePairs(params.reads, checkIfExists: true)

    // 3. Processamento
    trimmed_ch = FASTP(reads_ch)

    // 4. Alinhamento combinando reads com o bundle da referência
    bwa_input_ch = trimmed_ch.combine(reference_bundle)
    bam_ch = BWA_MEM(bwa_input_ch)

    // 5. Variant Calling
    call_input_ch = bam_ch.combine(reference_bundle)
    vcf_ch = BCFTOOLS_CALL_SAMPLE(call_input_ch)
    
    // 6. Filtro de SNPs
    snps_ch = VCF_FILTER_SNPS(vcf_ch)

    // 7. Merge de todos os VCFs
    all_vcfs = snps_ch.map{ it[1] }.collect()
    all_csis = snps_ch.map{ it[2] }.collect()
    merged_vcf = MERGE_VCFS(all_vcfs, all_csis)

    // 8. Converte VCF merged para FASTA de SNPs
    alignment = VCF_TO_FASTA(merged_vcf)

    // 9. Construção da Árvore com apenas SNPs
    FASTTREE(alignment)
}
