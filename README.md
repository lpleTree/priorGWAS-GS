# priorGWAS-GS: Multi-Model GWAS and Genomic Selection Integration for Poplar

## Overview

This repository contains the source code and input data for a multi-model Genome-Wide Association Study (GWAS) and Genomic Selection (GS) integration framework applied to 849 poplar (*Populus*) samples. The pipeline evaluates four GWAS models and three GS prediction models under a repeated cross-validation framework, incorporating GWAS-identified significant SNPs into GS predictions through multiple strategies to improve prediction accuracy.

## Background

Genomic selection has become a powerful tool in plant breeding, but its accuracy depends heavily on the genetic architecture of target traits. Integrating GWAS results—specifically, trait-associated SNPs—into GS models can potentially improve prediction performance. This study systematically evaluates:

1. Which GWAS model best captures the genetic architecture of poplar traits
2. How GWAS-derived SNP subsets at varying significance thresholds affect GS accuracy
3. Whether incorporating GWAS-significant SNPs as fixed effects improves prediction over standard GS models
4. The interaction between different GS models (sBLUP, gBLUP, cBLUP) and GWAS-assisted strategies

## Data Description

### Genotype Data
- **File**: `sample849_012tab.txt` (~90 MB, Git LFS tracked)
- **Format**: 012 genotype coding (0 = homozygous reference, 1 = heterozygous, 2 = homozygous alternative)
- **Samples**: 849 poplar individuals (rows)
- **Markers**: SNP markers (columns), with sample IDs in the first column (`Taxa`)

### SNP Map
- **File**: `GM.txt`
- **Format**: Tab-separated, three columns
  - `SNP`: SNP marker name (e.g., `chr01a_42236`)
  - `Chromosome`: Chromosome number
  - `Position`: Physical position (bp)

### Phenotype Data
- **File**: `杨树849个样本表型数据.txt` (Phenotype data for 849 poplar samples)
- **Format**: Tab-separated
  - `Taxa`: Sample ID
  - `trait1`–`trait6`: Six phenotypic traits (the analysis processes one trait at a time)

## Methods Workflow

### Phase 1: Multi-Model GWAS

Four GWAS models are run simultaneously on the full dataset:

| Model | Type | Description |
|-------|------|-------------|
| **FarmCPU** | Mixed Linear Model | Fixed and random model Circulating Probability Unification — controls false positives and negatives |
| **BLINK** | Bayesian | Bayesian-information and Linkage-disequilibrium Iteratively Nested Keyway |
| **MLM** | Mixed Linear Model | Standard mixed linear model with PCA + kinship |
| **GLM** | General Linear Model | General linear model with PCA only |

The model detecting the most significant associations is automatically selected as the **best model** for downstream GS analysis.

### Phase 2: GWAS-Assisted Genomic Selection

The core function `GWAS_assisted_GS()` implements a cross-validation framework with three GS models:

| Model | Full Name | Description |
|-------|-----------|-------------|
| **sBLUP** | Subset BLUP | BLUP using a subset of markers |
| **gBLUP** | Genomic BLUP | Standard genomic BLUP with kinship matrix |
| **cBLUP** | Compressed BLUP | BLUP with compressed relationship matrix |

#### Cross-Validation Design
- **K-fold**: 5-fold or 10-fold cross-validation
- **Repeats**: Configurable (e.g., 1–5 repeats)
- **Evaluation metric**: Pearson correlation coefficient (*r*) between predicted and observed phenotypes

#### GWAS Integration Strategies

Within each cross-validation fold, GWAS is re-run on the training set, and the following strategies are applied:

| Strategy | Description |
|----------|-------------|
| **Strategy 2** | Subset SNPs at different P-value thresholds (0.01, 0.03, 0.05, 0.07, 0.1, 0.2, 0.3, 0.4, 0.5, 1.0) from the best GWAS model; use selected SNPs for GS prediction |
| **Strategy 4** | The most significant SNP from the best GWAS model is used as a fixed-effect covariate; simultaneously, SNP subsets are built at different P-value thresholds for GS prediction (MAS + GS integration) |
| **Strategy 6** | Significant SNPs detected by **all four GWAS models** (union set) are used as fixed-effect covariates; SNP subsets from the best model at different P-value thresholds are used for GS prediction |

## Repository Structure

```
priorGWAS-GS/
├── multi_model_GWAS_GS_integration_cv.R  # Main analysis pipeline (187 lines)
├── gapit_functionsw.txt                   # Modified GAPIT R package functions
├── GM.txt                                 # SNP map (SNP, Chromosome, Position)
├── sample849_012tab.txt                   # Genotype matrix (012 coding) — Git LFS
├── 杨树849个样本表型数据.txt               # Phenotype data (849 samples × 6 traits)
└── .gitattributes                         # Git LFS tracking configuration
```

## Dependencies

### R Packages

```r
library(openxlsx)      # Excel file I/O
library(tidyverse)     # Data manipulation and visualization
library(data.table)    # Fast data reading
library(ggplot2)       # Plotting
library(caret)         # Cross-validation fold generation
library(dplyr)         # Data manipulation
library(Matrix)        # Sparse matrix operations
```

### GAPIT (Genome Association and Prediction Integrated Tool)

The pipeline uses a **modified version** of GAPIT functions loaded from `gapit_functionsw.txt`. Key modifications include:
- Custom output format handling
- Adjusted model parameter defaults
- Modified GWAS result file naming and parsing

> **Note**: The standard GAPIT package is available at [zzlab.net/GAPIT](http://zzlab.net/GAPIT/) or via the [GAPIT GitHub repository](https://github.com/jiabowang/GAPIT).

## Usage

### 1. Clone the Repository

```bash
git clone https://github.com/lpleTree/priorGWAS-GS.git
cd priorGWAS-GS
```

### 2. Ensure Git LFS is Installed

The genotype file `sample849_012tab.txt` is tracked with Git LFS. Install Git LFS before cloning or pull LFS files after cloning:

```bash
git lfs install
git lfs pull
```

### 3. Run the Analysis

Open `multi_model_GWAS_GS_integration_cv.R` in R or RStudio and modify the trait selection:

```r
myY = read.table(file = "杨树849个样本表型数据.txt", header = TRUE)[ ,c(1,2)]  # columns 1-2 for trait1
# For trait2, change to: [ ,c(1,3)]
# For trait3, change to: [ ,c(1,4)]
# ...
```

Then run the script. The pipeline will:
1. Execute multi-model GWAS (FarmCPU, BLINK, MLM, GLM)
2. Identify the best-performing GWAS model
3. Run GWAS-assisted GS with 5-fold / 10-fold cross-validation
4. Output prediction accuracy results as CSV files

### 4. Key Parameters

Within the function call at the bottom of the script:

```r
result <- GWAS_assisted_GS(GD, GM, myY,
    models = c("sBLUP", "gBLUP", "cBLUP"),
    p.levels = c(0.01, 0.03, 0.05, 0.07, 0.1, 0.2, 0.3, 0.4, 0.5, 1.0),
    n_repeats = 1,
    n_folds = 10)
```

- `models`: GS models to evaluate
- `p.levels`: P-value thresholds for SNP subsetting
- `n_repeats`: Number of cross-validation repeats
- `n_folds`: Number of cross-validation folds

### 5. Output

- **Phase 1 output**: `GAPIT.Association.Filter_GWAS_results.csv` — Significant SNPs from all four GWAS models
- **Phase 2 output**: `differ-stra/trait*.csv` — Prediction accuracy (*r*) for each strategy × model × P-value combination

## Notes

- The analysis processes **one trait at a time**; modify column indices in `myY` to switch traits
- The genotype file is ~90 MB and uses Git LFS; standard `git clone` will download a pointer file unless `git lfs pull` is run
- The pipeline automatically creates a `differ-stra/` subdirectory and changes working directory to it during GS analysis
- Random seed is set (`set.seed(123)`) for reproducibility

## Citation

If you use this code or data in your research, please cite the corresponding publication (to be updated upon acceptance) and the GAPIT software:

> Wang, J., & Zhang, Z. (2021). GAPIT Version 3: Boosting Power and Accuracy for Genomic Association and Prediction. *Genomics, Proteomics & Bioinformatics*, 19(4), 629–640.

## License

This project is shared for academic and research purposes. Please contact the authors for permissions regarding commercial use.
