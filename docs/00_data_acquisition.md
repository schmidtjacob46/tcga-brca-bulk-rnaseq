# 00 - Data Acquisition

## Overview

This project uses publicly available bulk RNA-sequencing data from **The Cancer Genome Atlas Breast Invasive Carcinoma (TCGA-BRCA)** cohort.

Raw gene expression count files were obtained from the **National Cancer Institute (NCI) Genomic Data Commons (GDC)** using the GDC Data Portal and the GDC Data Transfer Tool.

This document describes how the input data were obtained before beginning the computational analysis.

---

# Data Source

Database:
- National Cancer Institute (NCI) Genomic Data Commons (GDC)

Project:
- TCGA-BRCA (Breast Invasive Carcinoma)

Portal:
https://portal.gdc.cancer.gov/

---

# Dataset Selection

The following filters were applied within the GDC Data Portal.

| Category | Selection |
|-----------|-----------|
| Project | TCGA-BRCA |
| Data Category | Transcriptome Profiling |
| Data Type | Gene Expression Quantification |
| Workflow Type | STAR - Counts |

This returns harmonized RNA-seq gene count files generated using the GDC STAR Counts workflow.

---

# Files Downloaded

The following files were downloaded from the GDC Portal:

- GDC Manifest
- Repository metadata (JSON)
- STAR Counts gene expression files

The manifest contains the UUIDs used by the GDC Data Transfer Tool to download the expression files.

The repository metadata contains sample-level annotations including:

- TCGA sample barcode
- Sample type
- File UUID
- File name
- Case information

These metadata are later parsed to identify matched tumor-normal samples.

---

# Download Procedure

After selecting the files within the GDC Portal:

1. Add files to the cart.
2. Download the manifest file.
3. Download the accompanying repository metadata (JSON).
4. Install the GDC Data Transfer Tool.
5. Download all STAR count files using:

```bash
gdc-client download -m gdc_manifest.txt
```

The downloaded directory contains one folder per UUID, each containing a STAR Counts TSV file.

---

# Why STAR Counts?

The GDC STAR Counts workflow provides harmonized gene-level read counts generated using a standardized alignment and quantification pipeline.

Raw counts are appropriate for differential expression analysis using DESeq2 because DESeq2 models count data directly and performs its own normalization internally.

---

# Project Workflow

The remainder of this repository follows the workflow below:

```
GDC Portal
      │
      ▼
Download STAR Count Files
      │
      ▼
Prepare Metadata
      │
      ▼
Quality Control
      │
      ▼
Construct Count Matrix
      │
      ▼
Paired DESeq2 Analysis
      │
      ▼
Tumor-only VST Transformation
      │
      ▼
CAF / Cytotoxic Signature Analysis
      │
      ▼
PAM50 Subtype Analysis
```

---

# Repository Organization

The subsequent notebooks correspond to the following analysis steps:

| Notebook | Description |
|-----------|-------------|
| 01 | Prepare GDC STAR Count files |
| 02 | Prepare TCGA metadata |
| 03 | Sample quality control |
| 04 | Construct gene count matrix |
| 05 | Paired DESeq2 differential expression |
| 06 | Generate tumor-only VST expression matrix |
| 07 | Integrate PAM50 subtype annotations |
| 08 | Tumor signature analysis |

---

# Reproducibility Notes

Large raw sequencing files and intermediate expression matrices are **not included** in this GitHub repository because of their size.

This repository instead provides:

- analysis notebooks
- analysis scripts
- processed summary tables
- publication-quality figures
- documentation describing each analysis step

allowing the complete workflow to be reproduced from the publicly available TCGA data.