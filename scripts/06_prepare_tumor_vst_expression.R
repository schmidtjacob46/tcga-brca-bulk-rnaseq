```r
# TCGA-BRCA tumor-only variance-stabilized expression preparation
#
# Purpose:
#   Prepare a tumor-only variance-stabilized expression matrix for downstream
#   exploratory analyses, including PCA, heatmaps, correlations, and gene
#   signature scoring.
#
# Inputs:
#   1. Raw unstranded STAR gene-count matrix
#   2. Metadata containing one primary tumor sample per patient
#
# Outputs:
#   1. Tumor-only VST expression matrix as CSV
#   2. Tumor-only DESeqDataSet object as RDS
#   3. Tumor-only varianceStabilizingTransformation object as RDS
#
# Notes:
#   - Raw counts are used as input to DESeq2.
#   - One tumor per patient is retained to avoid pseudoreplication.
#   - A design of ~ 1 is used because no differential-expression comparison
#     is performed in this script.
#   - The VST matrix is intended for exploratory and signature-based analyses,
#     not for differential-expression testing.

library(DESeq2)

# 1. Read inputs ------------------------------------------------------------

counts <- read.csv(
  "data/processed/tcga_brca_raw_unstranded_counts_sample_level.csv",
  check.names = FALSE
)

metadata_tumor <- read.csv(
  "data/metadata/tcga_brca_tumor_one_per_patient.csv",
  check.names = FALSE
)

# 2. Validate cohort and sample identifiers ---------------------------------

cat("Tumor cohort summary:\n")
print(table(metadata_tumor$sample_group))

cat(
  "Unique patients:",
  length(unique(metadata_tumor$patient_barcode)),
  "\n"
)

cat(
  "Unique tumor samples:",
  length(unique(metadata_tumor$sample_barcode)),
  "\n"
)

# Confirm that each metadata sample is represented in the count matrix.
stopifnot(
  all(metadata_tumor$sample_barcode %in% colnames(counts))
)

# Confirm that the metadata contains one tumor sample per patient.
stopifnot(
  !anyDuplicated(metadata_tumor$patient_barcode),
  !anyDuplicated(metadata_tumor$sample_barcode)
)

# 3. Prepare the numeric count matrix ---------------------------------------

# DESeq2 expects genes as rows, samples as columns, and numeric counts as values.
# Move gene identifiers from the first column into the row names.
stopifnot("gene_id" %in% colnames(counts))

rownames(counts) <- counts$gene_id
counts$gene_id <- NULL

# Retain and order count-matrix columns according to the tumor metadata.
counts_tumor <- counts[, metadata_tumor$sample_barcode, drop = FALSE]

# Confirm exact agreement between the count matrix and metadata sample order.
stopifnot(
  identical(colnames(counts_tumor), metadata_tumor$sample_barcode)
)

cat(
  "Unfiltered tumor count matrix:",
  nrow(counts_tumor), "genes x",
  ncol(counts_tumor), "samples\n"
)

# 4. Construct the tumor-only DESeq2 object ---------------------------------

# The design contains only an intercept because this object is being used for
# normalization and transformation rather than differential-expression testing.
dds_tumor <- DESeqDataSetFromMatrix(
  countData = round(as.matrix(counts_tumor)),
  colData = metadata_tumor,
  design = ~ 1
)

# 5. Filter low-expression genes --------------------------------------------

# Keep genes with at least 10 counts in at least 10 tumor samples.
# This removes genes dominated by sparse count noise while retaining genes
# expressed in biologically relevant tumor subsets.
keep <- rowSums(counts(dds_tumor) >= 10) >= 10

cat(
  "Genes retained after filtering:",
  sum(keep), "of", length(keep), "\n"
)

dds_tumor <- dds_tumor[keep, ]

cat(
  "Filtered DESeq2 object:",
  nrow(dds_tumor), "genes x",
  ncol(dds_tumor), "samples\n"
)

# 6. Normalize counts and apply VST -----------------------------------------

# Estimate sample-specific size factors to account for differences in library
# depth and RNA composition between tumors.
dds_tumor <- estimateSizeFactors(dds_tumor)

# Apply the variance-stabilizing transformation.
# blind = TRUE is appropriate because there is no experimental design variable
# whose biological differences need to be preserved during dispersion fitting.
vsd_tumor <- vst(
  dds_tumor,
  blind = TRUE
)

vst_matrix <- assay(vsd_tumor)

stopifnot(
  identical(dim(vst_matrix), dim(dds_tumor)),
  identical(colnames(vst_matrix), metadata_tumor$sample_barcode)
)

cat(
  "Final VST matrix:",
  nrow(vst_matrix), "genes x",
  ncol(vst_matrix), "samples\n"
)

# 7. Save outputs ------------------------------------------------------------

write.csv(
  vst_matrix,
  "data/processed/tcga_brca_tumor_vst_expression.csv",
  quote = FALSE
)

saveRDS(
  dds_tumor,
  "data/processed/tcga_brca_tumor_only_dds.rds"
)

saveRDS(
  vsd_tumor,
  "data/processed/tcga_brca_tumor_only_vsd.rds"
)

cat("Tumor-only VST preparation completed successfully.\n")
```
