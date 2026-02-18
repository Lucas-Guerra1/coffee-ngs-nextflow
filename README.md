# SNP Phylogenomics Pipeline (Nextflow DSL2)

Pipeline automatizado para **detecção de SNPs, construção de alinhamento genômico e inferência filogenética** a partir de dados de sequenciamento paired-end (Illumina), utilizando **Nextflow DSL2**.

O workflow realiza:

* controle de qualidade de reads
* alinhamento ao genoma de referência
* chamada de variantes por amostra
* filtragem de SNPs
* merge multi-amostras
* conversão para alinhamento FASTA
* inferência de árvore filogenética

O pipeline é totalmente reprodutível, modular e escalável.

---

## Visão Geral do Workflow

O pipeline executa as seguintes etapas:

```
Reads → QC → Alinhamento → Variant Calling → SNP Filter → Merge → SNP Alignment → Phylogeny
```

### Etapas detalhadas

1. **Indexação da referência**

   * `bwa index`

2. **Controle de qualidade e trimming**

   * `fastp`
   * detecção automática de adaptadores

3. **Alinhamento ao genoma**

   * `bwa mem`
   * ordenação e indexação com `samtools`

4. **Chamada de variantes por amostra**

   * `bcftools mpileup`
   * `bcftools call`

5. **Filtragem de SNPs**

   * apenas variantes bialélicas

6. **Merge de VCFs multi-amostras**

   * `bcftools merge`

7. **Conversão de SNPs para alinhamento FASTA**

   * gera matriz SNP concatenada
   * substitui nomes genéricos pelos nomes reais das amostras

8. **Inferência filogenética**

   * `FastTree`
   * modelo GTR

---

## Estrutura de Saída

```
results/
├── 00_ref_index/        # Índices BWA da referência
├── 01_fastp/            # Reads filtrados
├── 02_bam/              # BAMs alinhados e indexados
├── 03_vcf/              # VCFs por amostra
├── 04_snps/             # VCFs contendo apenas SNPs
├── 05_merged_vcf/       # VCF multi-amostras
├── 06_alignment/        # Alinhamento SNP em FASTA
└── 07_tree/             # Árvore filogenética (Newick)
```

---

## Requisitos

### Software

* Nextflow ≥ 22
* Conda ou Mamba

### Ferramentas instaladas via ambiente Conda

* bwa
* fastp
* samtools
* bcftools
* FastTree

O pipeline utiliza:

```
envs/bioinfo.yml
```

para gerenciar dependências.

---

## Formato de Entrada

### Reads

O pipeline espera arquivos paired-end com padrão:

```
data/sample_1.fastq.gz
data/sample_2.fastq.gz
```

ou

```
data/SAMPLEID_R1.fastq.gz
data/SAMPLEID_R2.fastq.gz
```

Configurável via:

```
--reads
```

Padrão:

```
data/*_{1,2}.fastq.gz
```

---

### Genoma de Referência

Arquivo FASTA:

```
data/ref.fa
```

Configurável via:

```
--ref_fasta
```

---

## Execução

### Execução básica

```bash
nextflow run main.nf
```

---

### Execução com parâmetros explícitos

```bash
nextflow run main.nf \
  --reads "data/*_{1,2}.fastq.gz" \
  --ref_fasta data/ref.fa \
  --outdir results
```

---

### Execução com múltiplos CPUs

```bash
nextflow run main.nf -process.cpus 8
```

---

## Parâmetros

| Parâmetro     | Descrição                                 | Padrão                  |
| ------------- | ----------------------------------------- | ----------------------- |
| `--reads`     | Padrão glob dos arquivos FASTQ paired-end | `data/*_{1,2}.fastq.gz` |
| `--ref_fasta` | Genoma de referência                      | `data/ref.fa`           |
| `--outdir`    | Diretório de saída                        | `results`               |

O pipeline valida automaticamente parâmetros obrigatórios.

---

## Arquivos Gerados

### Variant Calling

* `*.sorted.bam`
* `*.vcf.gz`
* `*.snps.vcf.gz`

### Filogenia

* `snp_alignment.fasta` — matriz SNP concatenada
* `tree.nwk` — árvore filogenética

---

## Descrição Técnica

### Estratégia de SNP Alignment

O pipeline:

* extrai apenas posições SNP
* remove fases (`/` e `|`)
* usa primeiro alelo por posição
* substitui dados ausentes por `N`
* concatena SNPs para cada amostra
* gera FASTA multi-sequência

Esse método é adequado para:

* filogenômica bacteriana
* epidemiologia molecular
* análise populacional baseada em SNP

---

## Arquitetura do Pipeline

Implementado em **Nextflow DSL2** com:

* processos modulares independentes
* canais tipados
* paralelização automática
* reprodutibilidade via Conda
* execução local ou HPC compatível

---

## Reprodutibilidade

O pipeline garante reprodutibilidade via:

* versionamento de ambiente
* isolamento Conda
* rastreamento de execução Nextflow
* outputs determinísticos

---

## Casos de Uso

* filogenia bacteriana baseada em SNP
* análise comparativa de genomas
* vigilância epidemiológica
* estudos populacionais
* análise evolutiva

---

## Limitações

* assume dados Illumina paired-end
* utiliza apenas SNPs bialélicos
* não realiza filtragem avançada de qualidade de variantes
* não executa masking de regiões repetitivas
* não inclui filtragem de recombinação

---

## Melhorias Futuras Sugeridas

* filtros de qualidade de SNP (DP, MQ, QUAL)
* masking de regiões repetitivas
* detecção de recombinação
* suporte a dados long-read
* modelos filogenéticos adicionais
* relatório de QC automatizado

---

## Licença

Definir conforme necessidade do projeto.
