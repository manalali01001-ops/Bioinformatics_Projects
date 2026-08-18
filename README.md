# Bioinformatics & Computational Biology Projects

This repository contains two data analysis pipelines: a **Single-Cell RNA-Seq Quality Control & Preprocessing** pipeline and a **Genomic Sequence Pattern Mining & Analysis** toolkit.

---

## 🧬 Project 1: Single-Cell RNA-Seq Quality Control & Preprocessing

### Overview
An end-to-end single-cell RNA sequencing (scRNA-seq) workflow analyzing PBMC/T-cell datasets using **Scanpy** and **AnnData**. The pipeline filters low-quality cells, isolates Naive T cells, performs mitochondrial QC, and executes total-count normalization.

### 🛠️ Tech Stack
* **Language:** Python 3.x
* **Core Libraries:** `scanpy`, `anndata`, `pandas`, `numpy`, `matplotlib`
* **Data Format:** `.h5ad`

### Key Workflow Steps
1. **Data Ingestion:** Loaded raw count matrices (5,016 cells × 20,953 genes).
2. **Subpopulation Selection:** Filtered for Naive T-cell subsets (`cell.type == 'Naive'`).
3. **Quality Control:** Filtered doublets (`n_genes_by_counts < 2000`) and mitochondrial contamination (`pct_counts_mt < 5%`).
4. **Normalization:** Applied CP10K total-count normalization across $2,068$ passing cells.

---

## 🔬 Project 2: Genomic Sequence Pattern Mining Engine

### Overview
A computational biology toolkit built in Python to perform Open Reading Frame (ORF) identification, $k$-mer frequency profiling, sequence alignment matching, and FASTQ quality score evaluation across viral and multi-FASTA genomic datasets.

### 🛠️ Tech Stack
* **Language:** Python 3.x
* **Core Libraries:** `Biopython` (`Bio.SeqIO`), `numpy`, `collections`
* **Formats:** `.fasta`, `.fa`, `.fastq`

### Core Analyses
1. **ORF Finder:** Scans forward reading frames (Frames 1–3) for start (`ATG`) and stop codons.
2. **$k$-mer Frequency:** Sliding-window extraction ($k = 5, 6, 7, 12$) to quantify sequence repeats.
3. **Pattern Matching:** Mismatch-tolerant (up to 2 mm) pattern finder on the Bacteriophage Lambda genome.
4. **FASTQ Quality Control:** Converts Phred+33 ASCII scores to per-cycle mean quality profiles.
