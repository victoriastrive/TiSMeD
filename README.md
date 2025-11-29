# TiSMeD: Tissue-Specific Methylation and Expression Analysis Pipeline

This repository contains the source code for the **TiSMeD** database analysis pipeline. It provides a robust framework for identifying tissue-specific markers and housekeeping patterns from large-scale DNA methylation, transcriptomics, and proteomics datasets.

The pipeline utilizes a Specificity Measure (SPM_adjust) combined with a sigmoid-weighted scoring system (T-Score) to assess the confidence of tissue specificity, ensuring high-confidence biomarker discovery.

---

## Table of Contents

* [Key Features](#key-features)
* [Requirements](#requirements)
* [Directory Structure](#directory-structure)
* [Input Data Format](#input-data-format)

  * [Expression / Methylation Matrix (`.csv`)](#1-expression--methylation-matrix-csv)
  * [Metadata File (`count.txt`)](#2-metadata-file-counttxt)
* [Usage](#usage)

  * [1. Clone the Repository](#1-clone-the-repository)
  * [2. Prepare Data](#2-prepare-data)
  * [3. Run the Pipeline](#3-run-the-pipeline)
* [License](#license)

---

## Key Features

The pipeline is divided into two main analytical streams.

### 1. Methylation Analysis (ratio-based, 0–1)

* **TSHM (Tissue-Specific hypermethylation sites)**
  Detects genomic regions that are hypermethylated in specific tissues.

* **TSUM (Tissue-Specific hypomethylation sites)**
  Detects genomic regions that are hypomethylated in specific tissues.

* **Housekeeping Methylation Patterns**
  Identifies constitutively hypermethylated or low hypomethylated sites across all tissues.

### 2. Expression Analysis (abundance-based)

* **TSG/TSP (Tissue-Specific Genes/Proteins)**
  Detects tissue-specific expression patterns in transcriptomic or proteomic data.


---

## Requirements

* **R** (≥ 4.0)

* **Java Runtime Environment (JRE)**
  Required to run `SPM.jar`.

* **R packages**:

  ```r
  install.packages(c("data.table", "tidyverse", "readxl", "tidyr"))
  ```

---

## Directory Structure

```text
├── functions.R       # Core algorithm definitions (TSHM, TSLM, TSG/TSP, etc.)
├── main.R            # Main execution script (Run this file)
├── SPM.jar           # Java executable for SPM calculation (Required)
├── example
│   └── data/             # Input directory (See Data Format below)
│      └── methylome/
│      └── transcriptome/
└   └── results/          # Output directory (Auto-generated)
```

You may adapt the directory layout as needed, but the default scripts assume the structure above.

---

## Input Data Format

The pipeline expects:

* One or more `.csv` matrices (methylation, transcriptome, proteome)
* A metadata text file (`count.txt` )

### 1. Expression / Methylation Matrix (`.csv`)

* **Rows:** Probes, genes, or proteins

* **Columns:**

  * Column 1: Feature ID
    Examples: `probeID`, `gene_symbol`, etc.
  * Columns 2...N: Tissue samples (numeric values)

### 2. Metadata File (`count.txt`)

A tab-separated or space-separated file containing sample and tissue counts per dataset. It is used for T-Score weighting.

Example:

| dataset        | sample | tissue |
| -------------- | ------ | ------ |
| dataset_name_A | 100    | 15     |
| dataset_name_B | 50     | 8      |

---

## Usage

### 1. Clone the Repository

```bash
git clone https://github.com/victoriastrive/TiSMeD.git
cd TiSMeD
````

### 2. Prepare Data

Place your input data into the corresponding folders under the `example/` directory (or your own custom path):

* Methylation data:

  ```text
  example/data/methylome/
  ```

* Transcriptome data:

  ```text
  example/data/transcriptome/
  ```

* Proteome data (optional, same format as transcriptome):

  ```text
  example/data/proteome/
  ```

Ensure the metadata file (`count.txt`) is present in the corresponding data directory and matches the datasets used in the matrices.

### 3. Run the Pipeline

The pipeline is designed to be modular. Instead of running everything at once, you can choose which part to run (methylation only, expression only, or selected modules).

#### 3.1 Load Functions in R

In an R or RStudio session:

```r
setwd("path/to/TiSMeD")

source("functions.R")
source("main.R")
```

#### 3.2 Run Methylation Analysis (TSHM / TSUM / Housekeeping Methylation)

Run the full methylation pipeline (TSHM, TSUM, HKHM, HKUM):

```r
run_methylation_pipeline(
  base_dir      = getwd(),
  spm_cut       = 0.8,
  meth_data_dir = "example/data/methylome",
  spm_jar       = file.path(getwd(), "SPM.jar"),
  run_tshm      = TRUE,   # Tissue-specific hypermethylation sites
  run_tsum      = TRUE,   # Tissue-specific hypomethylation sites
  run_hkhm      = TRUE,   # Housekeeping hypermethylated sites
  run_hkum      = TRUE    # Housekeeping low hypomethylated sites
)
```

If you only want part of the methylation analysis, you can disable modules by setting the corresponding flags to `FALSE`.
For example, only TSHM and TSUM:

```r
run_methylation_pipeline(
  base_dir      = getwd(),
  spm_cut       = 0.8,
  meth_data_dir = "example/data/methylome",
  spm_jar       = file.path(getwd(), "SPM.jar"),
  run_tshm      = TRUE,
  run_tsum      = TRUE,
  run_hkhm      = FALSE,
  run_hkum      = FALSE
)
```

#### 3.3 Run Expression Analysis (TSG / TSP)

To run only the transcriptome/proteome analysis (TSG/TSP):

```r
run_expression_pipeline(
  base_dir       = getwd(),
  spm_cut        = 0.8,
  trans_data_dir = "example/data/transcriptome",
  spm_jar        = file.path(getwd(), "SPM.jar")
)
```

## License

This project is licensed under the **MIT License**.

See the [`LICENSE`](LICENSE) file for full license terms.

