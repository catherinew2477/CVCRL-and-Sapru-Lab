# download limma package

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("limma")

library(limma)

# prepare design matrix
sampleinfo <- read.table("siRNA_sampleinfo copy.txt", head = TRUE, row.names = 1)
group <- factor(sampleinfo$siRNA, levels = c("untreated", "treated"))
design <- model.matrix(~ 0 + group)
colnames(design) <- c("untreated", "treated")

# read in normalized counts, I think it isn't log2-transformed
norm_counts <- read.table("siRNA_norm_counts.txt", head = TRUE, sep = '\t', row.names = 1)
log2_counts <- log2(norm_counts+1)
write.csv(log2_counts,"log2transformed_knockdown.csv")

# fit the linear model
fit <- lmFit(log2_counts, design)

# specify contrasts --> treated vs untreated
contrast.matrix <- makeContrasts(
  Treatment_vs_Control = treated - untreated,
  levels = design
)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# extract DEGs
deg_results <- topTable(fit2, adjust="fdr", number=Inf)
significant_genes <- deg_results[deg_results$adj.P.Val < 0.05, ]
write.csv(significant_genes, "sig_genes_knockdown.csv")

# make volcano plot
library(org.Hs.eg.db)
library(tidyverse)
library(ggplot2)
library(AnnotationDbi)
library(ggrepel)


deg_results_colorcode <- deg_results %>%
  mutate(Change = case_when(
    adj.P.Val < 0.05 & logFC > 1 ~ "Increased",
    adj.P.Val < 0.05 & logFC <= -1 ~ "Decreased",
    TRUE ~ "No change"
  ))

genes_to_label <- c("COL8A1", "CCL2")

ggplot(res, aes(x = logFC, y = -log10(adj.P.Val))) + geom_point(aes(color = Change), alpha = 0.5) +
  scale_color_manual(values = c("Increased" = "blue", "Decreased" = "red", "No change" = "darkgrey")) +
  theme_minimal() + labs(x ="Log2 Fold Change", y = "Significance (-Log10)") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed")  +
  geom_text_repel(data = subset(res, X %in% genes_to_label), aes(label = X), size = 5, box.padding = 1, point.padding = 1)
  

# label genes, convert from Ensembl ID
gene_symbols <- mapIds(org.Hs.eg.db,
                       keys = rownames(deg_results_colorcode),
                       column = "SYMBOL",   # The column to get (e.g., gene symbol)
                       keytype = "ENSEMBL", # Type of the input key (e.g., Ensembl ID)
                       multiVals = "first")

gene_symbols_unique <- gene_symbols[!duplicated(gene_symbols)]
write.csv(deg_results_colorcode, "deg_results.csv")
res <- read.table("deg_results.csv", head = TRUE, sep = ',')
res$X <- gene_symbols_unique[match(res$X, names(gene_symbols_unique))]


write.csv(res, "labeled_res.csv")

# making heatmap
library(ComplexHeatmap)
significant_genes <- read.table("sig_genes_knockdown.csv", header = TRUE, sep = ',')
heatmap_genes <- significant_genes[(abs(significant_genes$logFC))>1.6,]
heatmap_genes_counts <- log2_counts[heatmap_genes$X,]
heatmap_genes_counts <- as.matrix(heatmap_genes_counts)
colnames(heatmap_genes_counts) <- c("siCTRL1", "siCTRL2", "siCTRL3", "siCOL8A1_1", "siCOL8A1_2", "siCOL8A1_3")
rownames(heatmap_genes_counts) <- gene_symbols_unique[match(rownames(heatmap_genes_counts), names(gene_symbols_unique))]
pheatmap(heatmap_genes_counts,
         scale = "row",                  # Scale data by row (z-score)
         cluster_rows = FALSE,           # Disable row clustering if needed
         cluster_cols = FALSE,           # Disable column clustering if needed
         color = colorRampPalette(c("red", "black", "green"))(50),
         show_rownames = TRUE,           # Show row names
         show_colnames = TRUE,
         name = "Z-score")

