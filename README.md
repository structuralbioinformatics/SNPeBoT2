# SNPeBoT2

SNPeBoT2 is a downloadable package for predicting the effect of SNPs on transcription factor binding and disease association.

---

## Table of Contents

- [Installation](#installation)
- [Usage](#usage)
  - [Predict SNP Effect](#predict_snp_effectsh)
  - [Run Association Analysis](#runassocsh)

---

## Installation

1. **Download this repository**

2. **Install ModCRElib** inside the SNPeBoT2 folder by following the instructions at:
   [https://github.com/structuralbioinformatics/ModCRElib](https://github.com/structuralbioinformatics/ModCRElib)
   - After installation, a `ModCRElib/` folder containing all necessary files should be present in the current directory
   - Ensure all ModCRElib dependencies are met

3. **Download required data files** into the `AssociationFilter/` folder and decompress them:
   - [TFLink_Homo_sapiens_interactions_All_simpleFormat_v1.0.tsv.gz](https://cdn.netbiol.org/tflink/download_files/TFLink_Homo_sapiens_interactions_All_simpleFormat_v1.0.tsv.gz)
   - [AllPredictions.AvgHiC.ABC0.015.minus150.ForABCPaperV3.txt.gz](https://mitra.stanford.edu/engreitz/oak/public/Nasser2021/AllPredictions.AvgHiC.ABC0.015.minus150.ForABCPaperV3.txt.gz)
   **Store the following in AssociationFilter/disgenet/**
   - [disgenet/all_gene_disease_pmid_associations.tsv.gz](https://www.disgenet.org) **An account is required to access disgenet (Free for academic purposes)**

---

## Usage

### `Predict_SNP_effect.sh`

Generates allele-specific binding (ASB) predictions for a list of reference and alternate sequences for a given transcription factor.

**Arguments:**

| # | Argument | Description |
|---|----------|-------------|
| 1 | `TF name` | Name or ID of the transcription factor (should match the ID used in arguments 3 and 4) |
| 2 | `job_id` | A unique job ID used to name output folders |
| 3 | `ASB file path` | Path to the file containing the SNPs to be tested (see format below) |
| 4 | `TF model` | Path to the model of the TF bound to DNA (can be generated with ModCRElib or another tool) |
| 5 | `template ID` | ID of the template used in the model. If ModCRElib was used, this can be found in the model filename. Otherwise, use `general` |

**ASB input file format** (4 columns per SNP, tab-separated):

| Column | Description |
|--------|-------------|
| 1 | TF name |
| 2 | Reference sequence (51 nt, mutation at position 26) |
| 3 | Alternate sequence (51 nt, mutation at position 26) |
| 4 | SNP ID (recommended format: `{chr}-{position}-{ref allele}-{alt allele}-{TF name}`) |

Example row:
```
FIGLA   GGATTTATCCATGTTTTTGCATGTAGAGATGGCTTGAAAAACAAACTATAT   GGATTTATCCATGTTTTTGCATGTACAGATGGCTTGAAAAACAAACTATAT   chr7-133924858-G-C-FIGLA
```

**Output:**
- Models submitted by the user are stored in `Output_{job_id}/models_{job_id}/`
- SNP binding effect predictions are stored in `Output_{job_id}/Predictions_{job_id}.tsv`

---

### `RunAssoc.sh`

Predicts whether a SNP is associated with a disease, based on the output of `Predict_SNP_effect.sh`.

**Arguments:**

| # | Argument | Description |
|---|----------|-------------|
| 1 | `TF name` | UniProt ID of the transcription factor |
| 2 | `Chrom` | Chromosome on which the SNP is located |
| 3 | `Prediction` | Predicted binding effect from `Predict_SNP_effect.sh` (`gain`, `loss`, or `no-change`) |
| 4 | `Position` | Genomic position of the SNP |
| 5 | `job_id` | The job ID used in `Predict_SNP_effect.sh` |
| 6 | `Assembly` | Human genome assembly used for the SNP location (e.g. `hg38`) |
| 7 | `SNP ID` | The SNP ID as used in `Predict_SNP_effect.sh` |

**Output:**
- Results are stored in `Output_{job_id}/associations/assoc_{SNP_ID}.csv`
