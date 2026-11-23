# Tbr1 project

## Summary

This repository contains the code used to analyze the data reported in the
study **“A transcriptional switch directs olfactory subsystem identity.”** The
analyses include 10x Genomics single-cell sequencing (scRNA-seq), bulk RNA-seq, ATAC-seq, and CUT&Tag.

## Data availability

The sequencing datasets generated for this study have been deposited in the
Gene Expression Omnibus (GEO) under accession number **GSE326062**.

## Software and package requirements

No non-standard hardware is required. Memory, storage, and runtime requirements
depend on sequencing depth and the number of cells or genomic features being
analyzed.

### 10x scRNA-seq

#### Preprocessing software

Raw reads were processed with Cell Ranger using the `cellranger count`
pipeline and the Single Cell 3' v3 chemistry. Intronic reads were not included
(`include-introns = false`). Reads were aligned to a custom GRCm38/mm10-based
mouse reference named `Mouse_EGFP_Cre`, which included EGFP and Cre sequences.

| Software | Version | Purpose |
|---|---:|---|
| Cell Ranger | 7.0.1 | Read alignment, barcode and UMI processing, and generation of gene–cell count matrices |

#### R packages

| Package | Version |
|---|---:|
| PRROC | 1.4 |
| rlang | 1.1.7 |
| pROC | 1.19.0.1 |
| LiblineaR | 2.10-24 |
| glmnet | 4.1-10 |
| Matrix | 1.6-4 |
| openxlsx | 4.2.6.1 |
| ComplexUpset | 1.3.3 |
| RColorBrewer | 1.1-3 |
| ggalt | 0.6.1 |
| patchwork | 1.2.0 |
| ggplot2 | 3.5.1 |
| tibble | 3.2.1 |
| stringr | 1.5.1 |
| purrr | 1.2.2 |
| tidyr | 1.3.1 |
| dplyr | 1.1.4 |
| Seurat | 5.2.0 |
| SeuratObject | 5.0.2 |
| sp | 2.1-4 |

`SeuratObject`, `rlang`, and `sp` are supporting dependencies. The code also
uses `sctransform` through `Seurat::SCTransform()` and formatting functions
from `scales`; their version numbers remain to be added.

### Bulk RNA-seq

#### Preprocessing software

Raw paired-end reads were assessed with FastQC before and after trimming.
Trimmomatic was used with TruSeq3 paired-end adapters and the parameters
`LEADING:15`, `TRAILING:15`, `SLIDINGWINDOW:4:15`, and `MINLEN:50`. Reads were
aligned to the GRCm38/mm10 primary assembly with HISAT2 using
`--rna-strandness RF`.

Alignments to chr1–chr19, chrX, and chrY were retained. Unmapped reads and
alignments with mapping quality below 20 were removed with SAMtools. Gene-level
counts were generated with featureCounts using the GENCODE vM25 comprehensive
gene annotation. Exonic fragments were assigned to `gene_id` with paired-end,
reverse-stranded counting (`-p -s 2`), exclusion of chimeric fragments (`-C`),
and fractional assignment of reads overlapping multiple features
(`-O --fraction`).

| Software | Version | Purpose |
|---|---:|---|
| FastQC | 0.12.1 | Read-quality assessment before and after trimming |
| Trimmomatic | 0.39 | Adapter and quality trimming |
| HISAT2 | 2.2.1 | Alignment to the GRCm38/mm10 reference genome |
| SAMtools | 1.16.1 | Alignment filtering, BAM sorting, and indexing |
| Subread/featureCounts | 2.0.1 | Gene-level read counting |

#### R environment and packages

Downstream analysis was performed in **R 4.2.0**. Gene-expression counts were
normalized and tested for differential expression with DESeq2.

| Package | Version |
|---|---:|
| DESeq2 | 1.38.3 |
| dplyr | 1.1.4 |
| tidyr | 1.3.0 |
| stringr | 1.5.1 |
| tibble | 3.2.1 |
| ggplot2 | 3.4.4 |
| ggrepel | 0.9.5 |
| ggpubr | 0.6.0 |
| RColorBrewer | 1.1.3 |
| scales | 1.3.0 |

### ATAC-seq

#### Preprocessing, peak calling, and visualization software

Raw paired-end reads were assessed with FastQC. fastp was used to remove
adapters, PolyG sequences, and low-quality reads using a base-quality threshold
of 20, a maximum of 40% low-quality bases per read, and a minimum read length
of 30 bp. Reads were aligned to the GRCm38/mm10 primary assembly with Bowtie2
in end-to-end, very-sensitive mode, with a maximum fragment length of 2,000 bp.

Alignments to chr1–chr19, chrX, and chrY were retained. Unmapped reads and
alignments with mapping quality below 20 were removed, and PCR duplicates were
removed with Picard. Peaks were called with Genrich in ATAC mode (`-j`), while
retaining unpaired alignments (`-y`), removing PCR duplicates (`-r`), excluding
chrM, and applying P-value and Q-value thresholds of 0.05.

deepTools was used to generate BPM-normalized bigWig files and signal
heatmaps. `bamCoverage` used an effective genome size of 2,652,783,500, a
10-bp bin size, and a 30-bp smoothing length. Signal matrices were calculated
around region centers using 2-kb upstream and downstream windows and 10-bp
bins.

| Software | Version | Purpose |
|---|---:|---|
| fastp | 0.23.2 | Adapter, PolyG, and quality trimming |
| FastQC | 0.12.1 | Read-quality assessment |
| Bowtie2 | 2.5.1 | Alignment to the GRCm38/mm10 reference genome |
| SAMtools | 1.16.1 | Alignment filtering, BAM sorting, and indexing |
| Picard MarkDuplicates | 2.18.29-SNAPSHOT | Removal of PCR duplicates |
| Genrich | 0.6.1 | ATAC-seq peak calling |
| deepTools | 3.5.2 | bigWig generation, signal matrices, and heatmaps |

### CUT&Tag

#### Preprocessing and visualization software

Raw paired-end reads were assessed with FastQC and trimmed with Trimmomatic
using Nextera paired-end adapters and the parameters `LEADING:15`,
`TRAILING:15`, `SLIDINGWINDOW:4:15`, and `MINLEN:50`. Reads were aligned to the
GRCm38/mm10 primary assembly with Bowtie2 in end-to-end, very-sensitive mode.
Mixed and discordant alignments were suppressed, and the permitted fragment
length was 10–700 bp.

Unmapped reads and alignments with mapping quality below 20 were removed, and
PCR duplicates were removed with Picard. deepTools was used to generate
BPM-normalized bigWig files and heatmaps using an effective genome size of
2,652,783,500, a 10-bp bin size, and a 30-bp smoothing length.

| Software | Version | Purpose |
|---|---:|---|
| FastQC | 0.11.7 | Raw-read quality assessment |
| Trimmomatic | 0.39 | Adapter and quality trimming |
| Java/JDK | 12.0.2_10 | Runtime for Trimmomatic |
| Bowtie2 | 2.5.1 | Alignment to the GRCm38/mm10 reference genome |
| SAMtools | 1.16.1 | Alignment filtering, BAM sorting, indexing, and regional depth calculation |
| Picard MarkDuplicates | 2.18.29-SNAPSHOT | Removal of PCR duplicates |
| deepTools | 3.5.2 | bigWig generation, signal matrices, and heatmaps |

#### Peak calling

CUT&Tag peaks were called with MACS2 in paired-end BAM mode (`-f BAMPE`) using
the corresponding IgG sample as the control. Peaks were identified with a
P-value threshold of 1 × 10⁻³ (`-p 1e-3`), and summit positions were reported
with `--call-summits`. Because PCR duplicates had already been removed, all
remaining fragments were retained with `--keep-dup all`.

| Software | Version | Purpose |
|---|---:|---|
| MACS2 | 2.2.6 | CUT&Tag peak calling |

#### Peak annotation: R environment and packages

Peak annotation was performed in **R 4.4.1** with ChIPseeker. A transcript
database was built from the *Mus musculus* Ensembl release 102 GTF annotation
(`Mus_musculus.GRCm38.102.gtf.gz`) for GRCm38. BED coordinates were converted
from 0-based to 1-based GRanges coordinates, and chromosome naming styles were
harmonized. Peaks were annotated relative to transcription start sites using a
window of −3 kb to +3 kb, with `org.Mm.eg.db` used for mouse gene annotation.

| Package | Version |
|---|---:|
| ChIPseeker | 1.42.0 |
| txdbmaker | 1.2.0 |
| AnnotationDbi | 1.68.0 |
| GenomicRanges | 1.58.0 |
| IRanges | 2.40.0 |
| GenomeInfoDb | 1.42.0 |
| org.Mm.eg.db | 3.20.0 |
