# TCGA-BRCA paired tumor-normal DESeq2 analysis
# Input: raw GDC STAR Counts matrix + matched tumor-normal metadata
# Output: paired tumor-vs-normal differential expression results
#
# This analysis uses only patients with both primary tumor and solid normal samples.
# The paired design controls for patient-specific baseline expression differences.
# Positive log2FoldChange values indicate higher expression in tumor relative to normal.
# Negative log2FoldChange values indicate higher expression in normal relative to tumor.

library(DESeq2)
library(tidyverse)
library(ggrepel)

# 1. Read inputs -----------------------------------------------------------

counts <- read.csv(
  "data/processed/tcga_brca_raw_unstranded_counts_sample_level.csv",
  row.names = 1,
  check.names = FALSE
)

metadata <- read.csv(
  "data/metadata/tcga_brca_matched_tumor_normal_one_per_group.csv",
  check.names = FALSE
)

gene_annot <- read.csv(
  "data/processed/tcga_brca_gene_annotation_from_star_counts.csv",
  check.names = FALSE
)

# 2. Subset to matched samples --------------------------------------------

# DESeq2 requires that count matrix columns and metadata rows are in the exact same sample order.
matched_samples <- metadata$sample_barcode
counts_matched <- counts[, matched_samples]

stopifnot(all(colnames(counts_matched) == metadata$sample_barcode))

# 3. Set factor levels -----------------------------------------------------

# Set normal as the reference level so the model coefficient is tumor vs normal.
metadata$sample_group <- factor(metadata$sample_group, levels = c("normal", "tumor"))
metadata$patient_barcode <- factor(metadata$patient_barcode)

print(table(metadata$sample_group))
print(length(unique(metadata$patient_barcode)))

# 4. Create DESeq2 object --------------------------------------------------

# Paired design:
# patient_barcode accounts for matched patient identity.
# sample_group estimates the tumor-vs-normal expression difference.
dds <- DESeqDataSetFromMatrix(
  countData = round(as.matrix(counts_matched)),
  colData = metadata,
  design = ~ patient_barcode + sample_group
)

# 5. Filter low-count genes ------------------------------------------------

# Remove genes with very low counts across samples to reduce noise and improve model stability.
keep <- rowSums(counts(dds) >= 10) >= 10
dds <- dds[keep, ]

print(dds)

# 6. Run DESeq2 ------------------------------------------------------------

dds <- DESeq(dds)

# Save the fitted DESeq2 object so the model does not need to be rerun later.
dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)

saveRDS(
  dds,
  "data/processed/tcga_brca_paired_tumor_vs_normal_deseq2_dds.rds"
)

# 7. Extract tumor-vs-normal results --------------------------------------

# Extract the tumor-vs-normal coefficient from the fitted model.
# This coefficient estimates tumor expression relative to normal after accounting for patient identity.
res <- results(dds, name = "sample_group_tumor_vs_normal")

summary(res)

# betaConv records whether the model coefficient converged for each gene.
# Non-converged genes are retained in the output but excluded from summary counts and interpretation.
res_df <- as.data.frame(res) %>%
  rownames_to_column("gene_id") %>%
  left_join(
    data.frame(
      gene_id = rownames(dds),
      betaConv = mcols(dds)$betaConv
    ),
    by = "gene_id"
  ) %>%
  left_join(gene_annot, by = "gene_id") %>%
  arrange(padj)

# 8. Summarize differential expression ------------------------------------

sig_summary <- tibble(
  threshold = c(
    "padj < 0.1",
    "padj < 0.05",
    "padj < 0.05 and abs(log2FC) > 1",
    "padj < 0.05 and abs(log2FC) > 2"
  ),
  n_genes = c(
    sum(res_df$padj < 0.1 & res_df$betaConv == TRUE, na.rm = TRUE),
    sum(res_df$padj < 0.05 & res_df$betaConv == TRUE, na.rm = TRUE),
    sum(res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1 & res_df$betaConv == TRUE, na.rm = TRUE),
    sum(res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 2 & res_df$betaConv == TRUE, na.rm = TRUE)
  )
)

sig_summary

# Save result tables.
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)

write.csv(
  res_df,
  "results/tables/tcga_brca_paired_tumor_vs_normal_deseq2_results.csv",
  row.names = FALSE
)

write.csv(
  sig_summary,
  "results/tables/tcga_brca_paired_tumor_vs_normal_deseq2_sig_summary.csv",
  row.names = FALSE
)

# Save top genes as a quick biological sanity-check table.
top20_genes <- res_df %>%
  filter(betaConv == TRUE, padj < 0.05) %>%
  arrange(padj) %>%
  select(gene_id, gene_name, gene_type, baseMean, log2FoldChange, padj) %>%
  head(20)

top20_genes

write.csv(
  top20_genes,
  "results/tables/tcga_brca_paired_tumor_vs_normal_top20_genes.csv",
  row.names = FALSE
)

# 9. Volcano plot ----------------------------------------------------------

# Volcano plot:
# x-axis shows tumor-vs-normal effect size.
# y-axis shows statistical significance.
# padj is floored at 1e-300 only for plotting to avoid infinite -log10 values.
volcano_df <- res_df %>%
  filter(betaConv == TRUE, !is.na(padj)) %>%
  mutate(
    padj_plot = pmax(padj, 1e-300),
    neglog10_padj = -log10(padj_plot),
    sig_category = case_when(
      padj < 0.05 & log2FoldChange > 1  ~ "Higher in tumor",
      padj < 0.05 & log2FoldChange < -1 ~ "Higher in normal",
      TRUE ~ "Other"
    )
  )

top_labels <- volcano_df %>%
  filter(padj < 0.05, abs(log2FoldChange) > 2) %>%
  arrange(padj) %>%
  slice_head(n = 15)

volcano_plot <- ggplot(volcano_df, aes(x = log2FoldChange, y = neglog10_padj, color = sig_category)) +
  geom_point(alpha = 0.5, size = 1) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_text_repel(
    data = top_labels,
    aes(label = gene_name),
    size = 3,
    max.overlaps = 20,
    show.legend = FALSE
  ) +
  labs(
    title = "TCGA-BRCA paired tumor vs normal",
    x = "Log2 fold change (tumor vs normal)",
    y = "-Log10 adjusted p-value",
    color = "Category"
  ) +
  theme_minimal()

volcano_plot

dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)

ggsave(
  "results/figures/tcga_brca_paired_tumor_vs_normal_volcano.png",
  plot = volcano_plot,
  width = 8,
  height = 6,
  dpi = 300
)
