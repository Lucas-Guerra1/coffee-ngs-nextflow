# SNP Phylogenomics Pipeline (Nextflow DSL2)

An automated pipeline for **SNP detection, genome alignment construction, and phylogenetic inference** from paired-end sequencing data (Illumina), implemented using **Nextflow DSL2**.

The workflow performs:

* read quality control
* reference genome alignment
* per-sample variant calling
* SNP filtering
* multi-sample VCF merging
* FASTA alignment generation
* phylogenetic tree inference

The pipeline is fully reproducible, modular, and scalable.

---

## Workflow Overview

The pipeline executes the following steps:

```
Reads → QC → Alignment → Variant Calling → SNP Filter → Merge → SNP Alignment → Phylogeny
```

### Detailed Steps

1. **Reference indexing**

   * `bwa index`

2. **Quality control and trimming**

   * `fastp`
   * automatic adapter detection

3. **Genome alignment**

   * `bwa mem`
   * sorting and indexing with `samtools`

4. **Per-sample variant calling**

   * `bcftools mpileup`
   * `bcftools call`

5. **SNP filtering**

   * retains only biallelic variants

6. **Multi-sample VCF merging**

   * `bcftools merge`

7. **Conversion of SNPs to FASTA alignment**

   * generates concatenated SNP matrix
   * replaces generic sample names with actual sample identifiers

8. **Phylogenetic inference**

   * `FastTree`
   * GTR model

---

## Output Structure

```
results/
├── 00_ref_index/        # BWA reference indices
├── 01_fastp/            # Filtered reads
├── 02_bam/              # Aligned and indexed BAM files
├── 03_vcf/              # Per-sample VCF files
├── 04_snps/             # SNP-only VCF files
├── 05_merged_vcf/       # Multi-sample VCF
├── 06_alignment/        # SNP FASTA alignment
└── 07_tree/             # Phylogenetic tree (Newick)
```

---

## Requirements

### Software

* Nextflow ≥ 22
* Conda or Mamba

### Tools installed via Conda environment

* bwa
* fastp
* samtools
* bcftools
* FastTree

The pipeline uses:

```
envs/bioinfo.yml
```

to manage dependencies.

---

## Input Format

### Reads

The pipeline expects paired-end FASTQ files with naming pattern:

```
data/sample_1.fastq.gz
data/sample_2.fastq.gz
```

or

```
data/SAMPLEID_R1.fastq.gz
data/SAMPLEID_R2.fastq.gz
```

Configurable via:

```
--reads
```

Default:

```
data/*_{1,2}.fastq.gz
```

---

### Reference Genome

FASTA file:

```
data/ref.fa
```

Configurable via:

```
--ref_fasta
```

---

## Execution

### Basic run

```bash
nextflow run main.nf
```

---

### Run with explicit parameters

```bash
nextflow run main.nf \
  --reads "data/*_{1,2}.fastq.gz" \
  --ref_fasta data/ref.fa \
  --outdir results
```

---

### Run with multiple CPUs

```bash
nextflow run main.nf -process.cpus 8
```

---

## Parameters

| Parameter     | Description                             | Default                 |
| ------------- | --------------------------------------- | ----------------------- |
| `--reads`     | Glob pattern for paired-end FASTQ files | `data/*_{1,2}.fastq.gz` |
| `--ref_fasta` | Reference genome                        | `data/ref.fa`           |
| `--outdir`    | Output directory                        | `results`               |

The pipeline automatically validates required parameters.

---

## Generated Files

### Variant Calling Outputs

* `*.sorted.bam`
* `*.vcf.gz`
* `*.snps.vcf.gz`

### Phylogenetic Outputs

* `snp_alignment.fasta` — concatenated SNP matrix
* `tree.nwk` — phylogenetic tree

---

## Technical Description

### SNP Alignment Strategy

The pipeline:

* extracts only SNP positions
* removes phasing symbols (`/` and `|`)
* uses the first allele per position
* replaces missing data with `N`
* concatenates SNPs per sample
* generates multi-sequence FASTA alignment

This approach is suitable for:

* bacterial phylogenomics
* molecular epidemiology
* SNP-based population analysis

---

## Pipeline Architecture

Implemented in **Nextflow DSL2** with:

* independent modular processes
* typed channels
* automatic parallelization
* reproducibility via Conda environments
* compatibility with local and HPC execution

---

## Reproducibility

The pipeline ensures reproducibility through:

* environment versioning
* Conda isolation
* Nextflow execution tracking
* deterministic outputs

---

## Use Cases

* SNP-based bacterial phylogeny
* comparative genomics
* epidemiological surveillance
* population studies
* evolutionary analysis

---

## Limitations

* assumes Illumina paired-end data
* retains only biallelic SNPs
* does not perform advanced variant quality filtering
* does not mask repetitive regions
* does not include recombination filtering

---

## Suggested Future Improvements

* SNP quality filtering (DP, MQ, QUAL thresholds)
* repetitive region masking
* recombination detection
* long-read data support
* additional phylogenetic models
* automated QC reporting

---

## License

MIT License

Copyright (c) 2026
