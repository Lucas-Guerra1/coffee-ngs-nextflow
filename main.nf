nextflow.enable.dsl=2

// Default parameters (avoid undefined warnings)
params.samplesheet = 'samplesheet.csv'
params.outdir      = 'results/00_test'

process SAMPLE_PING {
  tag "$sample_id"
  publishDir params.outdir, mode: 'copy'

  input:
    tuple val(sample_id), path(r1), path(r2)

  output:
    path("${sample_id}.txt")

  script:
  """
  echo "Sample: ${sample_id}" > ${sample_id}.txt
  echo "R1: ${r1}" >> ${sample_id}.txt
  echo "R2: ${r2}" >> ${sample_id}.txt
  """
}

workflow {
  Channel
    .fromPath(params.samplesheet)
    .splitCsv(header:true)
    .map { row -> tuple(row.sample, file(row.fastq_1), file(row.fastq_2)) }
    | SAMPLE_PING
}
