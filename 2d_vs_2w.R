setwd("/Users/catherinewang/Desktop/CVCRL/Emory mouse studies/")

#### LIBRARIES ####
library(Seurat)
library(dplyr)
library(ggplot2)
library(ggpubr)


# For GO analysis
library(tibble)
library(dplyr)
library(clusterProfiler)
library(org.Mm.eg.db)
library(ComplexHeatmap)
library(circlize)
library(forcats)
library(patchwork)


#### CREATE SEURAT OBJECTS FOR EACH SAMPLE ####

# load data
lca_2d_data <- Read10X("CellReports_scRNAseq/lca_2d_pcl_filtered_feature_bc_matrix/")
lca_2wk_data <- Read10X("CellReports_scRNAseq/lca_2wk_pcl_filtered_feature_bc_matrix/")
rca_2d_data <- Read10X("CellReports_scRNAseq/rca_2d_pcl_filtered_feature_bc_matrix/")
rca_2wk_data <- Read10X("CellReports_scRNAseq/rca_2wk_pcl_filtered_feature_bc_matrix/")

# create
seurat_lca_2d <- CreateSeuratObject(counts = lca_2d_data, project = "lca_2d")
seurat_lca_2wk <- CreateSeuratObject(counts = lca_2wk_data, project = "lca_2wk")
seurat_rca_2d <- CreateSeuratObject(counts = rca_2d_data, project = "rca_2d")
seurat_rca_2wk <- CreateSeuratObject(counts = rca_2wk_data, project = "rca_2wk")

#### ADD METADATA ####

# type of flow
seurat_lca_2d$condition <- "disturbed"
seurat_lca_2wk$condition <- "disturbed"
seurat_rca_2d$condition <- "normal"
seurat_rca_2wk$condition <- "normal"

# 2 days or 2 weeks?
seurat_lca_2d$time <- "days"
seurat_lca_2wk$time <- "weeks"
seurat_rca_2d$time <- "days"
seurat_rca_2wk$time <- "weeks"

#### COMBINE DATASETS ####
combined_seurat <- merge(seurat_lca_2d, y = list(seurat_lca_2wk, seurat_rca_2d, seurat_rca_2wk))

#### QC AND UMAP ####

combined_seurat[["percent.mt"]] <- PercentageFeatureSet(combined_seurat, pattern = "^MT-")
combined_seurat <- subset(combined_seurat, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 10)

# Normalize and find variable features
combined_seurat <- NormalizeData(combined_seurat)
combined_seurat <- FindVariableFeatures(combined_seurat)

# Scaling 
combined_seurat <- ScaleData(combined_seurat)

# Run PCA (dimensionality reduction)
combined_seurat <- RunPCA(combined_seurat, npcs = 100)

# UMAP and clustering
set.seed(42)
combined_seurat <- RunUMAP(combined_seurat, dims = 1:100)
combined_seurat <- FindNeighbors(combined_seurat, dims = 1:100)
combined_seurat <- FindClusters(combined_seurat, resolution = 0.5)

# Plot UMAP (6375 cells)
umap <- DimPlot(combined_seurat, reduction = "umap", label = TRUE) 
umap

# Number of cells in each cluster
table(Idents(combined_seurat))


# Plot UMAP distinguishing between diseased and healthy cells
combined_seurat$sample <- paste(combined_seurat$condition, combined_seurat$time, sep = "_")

UMAP_samples <- DimPlot(combined_seurat, reduction = "umap", group.by = "sample") +
        ggtitle("UMAP of cell clusters from each sample")
UMAP_samples

UMAP_condition <- DimPlot(combined_seurat, reduction = "umap", group.by = "condition") +
  ggtitle("UMAP of disturbed vs. normal conditions")
UMAP_condition

#### CLUSTER 2 SUBSET ANALYSIS (725 cells) ####
cluster2_seurat_subset <- subset(combined_seurat, idents = 2)

cluster2_counts_lca_2d <- GetAssayData(cluster2_seurat_subset, layer = "counts.lca_2d")
cluster2_counts_lca_2wk <- GetAssayData(cluster2_seurat_subset, layer = "counts.lca_2wk")
cluster2_counts_rca_2d <- GetAssayData(cluster2_seurat_subset, layer = "counts.rca_2d")
cluster2_counts_rca_2wk <- GetAssayData(cluster2_seurat_subset, layer = "counts.rca_2wk")

cluster2_subset_counts <- cbind(cluster2_counts_lca_2d, cluster2_counts_lca_2wk,
                                cluster2_counts_rca_2d, cluster2_counts_rca_2wk)

cluster2_subset_metadata <- cluster2_seurat_subset@meta.data

cluster2_seurat_subset <- CreateSeuratObject(
  counts = cluster2_subset_counts,
  meta.data = cluster2_subset_metadata,
  project = "SMC_Subset_UMAP"
)

# Normalize, find variable features, scale, run PCA
cluster2_seurat_subset <- NormalizeData(cluster2_seurat_subset)
cluster2_seurat_subset <- FindVariableFeatures(cluster2_seurat_subset, selection.method = "vst", nfeatures = 2000)
cluster2_seurat_subset <- ScaleData(cluster2_seurat_subset, features = rownames(cluster2_seurat_subset))
cluster2_seurat_subset <- RunPCA(cluster2_seurat_subset, features = VariableFeatures(object = cluster2_seurat_subset))


# Find clusters
ElbowPlot(combined_seurat)
cluster2_seurat_subset <- FindNeighbors(cluster2_seurat_subset, dims = 1:15)
cluster2_seurat_subset <- FindClusters(cluster2_seurat_subset, resolution = 1.5)

# Run UMAP
cluster2_seurat_subset <- RunUMAP(cluster2_seurat_subset, dims = 1:10)

cluster2_seurat_subset_umap <- DimPlot(cluster2_seurat_subset, reduction = "umap", label = TRUE) + ggtitle("UMAP of Cluster 2 Subset")
cluster2_seurat_subset_umap

#### GO ANALYSIS - Cluster 2 vs ECs (0,1,8,5,11) ####
cluster2_ECs_subset <- subset(combined_seurat, idents = (c(0,1,8,5,11,2)))

cluster2_ECs_counts_lca_2d <- GetAssayData(cluster2_ECs_subset, layer = "counts.lca_2d")
cluster2_ECs_counts_lca_2wk <- GetAssayData(cluster2_ECs_subset, layer = "counts.lca_2wk")
cluster2_ECs_counts_rca_2d <- GetAssayData(cluster2_ECs_subset, layer = "counts.rca_2d")
cluster2_ECs_counts_rca_2wk <- GetAssayData(cluster2_ECs_subset, layer = "counts.rca_2wk")

cluster2_ECs_subset_counts <- cbind(cluster2_ECs_counts_lca_2d, cluster2_ECs_counts_lca_2wk, 
                                    cluster2_ECs_counts_rca_2d, cluster2_ECs_counts_rca_2wk)

# Extract the metadata for the subset cells
cluster2_ECs_subset_metadata <- combined_seurat@meta.data

# Create a new Seurat object
cluster2_ECs_seurat_subset <- CreateSeuratObject(
  counts = cluster2_ECs_subset_counts,
  meta.data = cluster2_ECs_subset_metadata,
  project = "Subset_UMAP"
)

# Normalize, find variable features, scale, run PCA
cluster2_ECs_seurat_subset <- NormalizeData(cluster2_ECs_seurat_subset)
cluster2_ECs_seurat_subset <- FindVariableFeatures(cluster2_ECs_seurat_subset, selection.method = "vst", nfeatures = 2000)
cluster2_ECs_seurat_subset <- ScaleData(cluster2_ECs_seurat_subset, features = rownames(cluster2_ECs_seurat_subset))
cluster2_ECs_seurat_subset <- RunPCA(cluster2_ECs_seurat_subset, features = VariableFeatures(object = cluster2_ECs_seurat_subset))

# Find clusters
cluster2_ECs_seurat_subset <- FindNeighbors(cluster2_ECs_seurat_subset, dims = 1:10)
cluster2_ECs_seurat_subset <- FindClusters(cluster2_ECs_seurat_subset, resolution = 0.5)

# Run UMAP
set.seed(42)
cluster2_ECs_seurat_subset <- RunUMAP(cluster2_ECs_seurat_subset, dims = 1:10)

# Plot UMAP
cluster2_ECs_seurat_subset_umap <- DimPlot(cluster2_ECs_seurat_subset, reduction = "umap", label = TRUE) + ggtitle("UMAP of Cluster 2 + ECs Subset")
cluster2_ECs_seurat_subset_umap

# Show which cells are PECAM+CDH5+
threshold <- 0
pecam1_positive_cells <- WhichCells(cluster2_ECs_seurat_subset, expression = Pecam1 > threshold)
cdh5_positive_cells <- WhichCells(cluster2_ECs_seurat_subset, expression = Cdh5 > threshold)

pecam1_cdh5_positive_cells <- intersect(pecam1_positive_cells, cdh5_positive_cells)

# Recreate a metadata column to label PECAM1/CDH5 positive cells
cluster2_ECs_seurat_subset$PECAM1_CDH5_Positive <- ifelse(Cells(cluster2_ECs_seurat_subset) %in% pecam1_cdh5_positive_cells, "PECAM1+CDH5+", "Other")
Idents(cluster2_ECs_seurat_subset) <- "PECAM1_CDH5_Positive"
cluster2_ECs_positive_umap <- DimPlot(cluster2_ECs_seurat_subset, group.by = "PECAM1_CDH5_Positive", label = TRUE) +
  ggtitle("UMAP of PECAM1 and CDH5 Positive Cells")
cluster2_ECs_positive_umap 

# Find DEGS - 2 and 5 vs everything else in subset

# Merge 2 and 5 into one cluster to represent the non-ECs in this subset
cluster2_ECs_seurat_subset$merged_clusters <- cluster2_ECs_seurat_subset$seurat_clusters
cluster2_ECs_seurat_subset$merged_clusters[cluster2_ECs_seurat_subset$seurat_clusters %in% c("2","5")] <- "2"

DimPlot(cluster2_ECs_seurat_subset, group.by = "merged_clusters", label = TRUE) + ggtitle("EC and non-EC Clusters UMAP")


Idents(cluster2_ECs_seurat_subset) <- "seurat_clusters"
cluster2_ECs_degs <- FindMarkers(object = cluster2_ECs_seurat_subset, ident.1 = "2",
                          min.pct = 0.25, logfc.threshold = 0.25)
filtered_cluster2_ECs_degs <- cluster2_ECs_degs %>% filter(p_val_adj < 0.05 & avg_log2FC > 0.25)
cluster2_ECs_gene_list <- filtered_cluster2_ECs_degs %>%
  rownames_to_column("gene") %>%
  dplyr::select(gene, avg_log2FC) %>%
  arrange(desc(avg_log2FC))

# Save DEGs
write.csv(filtered_cluster2_ECs_degs, "cluster2_vs_ECs_degs.csv")

# Convert gene symbols to Entrez IDs
cluster2_ECs_gene_list$entrez <- mapIds(org.Mm.eg.db, keys = cluster2_ECs_gene_list$gene, 
                                 column = "ENTREZID", keytype = "SYMBOL", 
                                 multiVals = "first")

# Run enrichment analysis (Gene Ontology - Biological Processes)
cluster2_ECs_enrichment_res <- enrichGO(
  gene          = na.omit(cluster2_ECs_gene_list$entrez),
  OrgDb         = org.Mm.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  pAdjustMethod = "BH",
  qvalueCutoff  = 0.05
)

calculate_fold_enrichment <- function(GeneRatio, BgRatio) {
  gene_numerator <- as.numeric(sub("/.*", "", GeneRatio))
  gene_denominator <- as.numeric(sub(".*/", "", GeneRatio))
  
  bg_numerator <- as.numeric(sub("/.*", "", BgRatio))
  bg_denominator <- as.numeric(sub(".*/", "", BgRatio))
  
  # Compute (GeneRatio) / (BgRatio)
  result <- (gene_numerator / gene_denominator) / (bg_numerator / bg_denominator)
  return (result)
}

cluster2_ECs_enrichment_data <- as.data.frame(cluster2_ECs_enrichment_res) %>%
  mutate(
    FoldEnrichment = calculate_fold_enrichment(GeneRatio, BgRatio),
    logFDR = -log10(p.adjust),
    GeneCount = as.numeric(sub("/.*", "", GeneRatio))
  ) %>%
  arrange(desc(logFDR)) %>%
  head(20)

ggplot(cluster2_ECs_enrichment_data, aes(x = FoldEnrichment, y = fct_reorder(Description, FoldEnrichment))) +
  geom_segment(aes(xend = 0, yend = Description), color = "gray") +
  geom_point(aes(size = GeneCount, color = logFDR)) +
  scale_color_gradient(low = "blue", high = "red") +
  scale_size(range = c(3, 8)) +
  labs(
    x = "Fold Enrichment",
    y = "",
    color = "-log10(FDR)",
    size = "N. of Genes"
  ) +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 12))

# Run GO-MF analysis
cluster2_ECs_MF_res <- enrichGO(
  gene          = na.omit(cluster2_ECs_gene_list$entrez),
  OrgDb         = org.Mm.eg.db,
  keyType       = "ENTREZID",
  ont           = "MF",
  pAdjustMethod = "BH",
  qvalueCutoff  = 0.05
)
# number of genes that mapped to term 
length(na.omit(cluster2_ECs_gene_list$entrez))

View(as.data.frame(cluster2_ECs_MF_res@result))

cluster2_ECs_MF_data <- as.data.frame(cluster2_ECs_MF_res) %>%
  mutate(
    FoldEnrichment = calculate_fold_enrichment(GeneRatio, BgRatio),
    logFDR = -log10(p.adjust),
    GeneCount = as.numeric(sub("/.*", "", GeneRatio))
  ) %>%
  arrange(desc(logFDR)) %>%
  head(20)



ggplot(cluster2_ECs_MF_data, aes(x = FoldEnrichment, y = fct_reorder(Description, FoldEnrichment))) +
  geom_segment(aes(xend = 0, yend = Description), color = "gray") +
  geom_point(aes(size = GeneCount, color = logFDR)) +
  scale_color_gradient(low = "blue", high = "red") +
  scale_size(range = c(3, 8)) +
  labs(
    x = "Fold Enrichment",
    y = "",
    color = "-log10(FDR)",
    size = "N. of Genes"
  ) +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 12))

# Run GO-CC analysis
cluster2_ECs_CC_res <- enrichGO(
  gene          = na.omit(cluster2_ECs_gene_list$entrez),
  OrgDb         = org.Mm.eg.db,
  keyType       = "ENTREZID",
  ont           = "CC",
  pAdjustMethod = "BH",
  qvalueCutoff  = 0.05
)

cluster2_ECs_CC_data <- as.data.frame(cluster2_ECs_CC_res) %>%
  mutate(
    FoldEnrichment = calculate_fold_enrichment(GeneRatio, BgRatio),
    logFDR = -log10(p.adjust),
    GeneCount = as.numeric(sub("/.*", "", GeneRatio))
  ) %>%
  arrange(desc(logFDR)) %>%
  head(20)

ggplot(cluster2_ECs_CC_data, aes(x = FoldEnrichment, y = fct_reorder(Description, FoldEnrichment))) +
  geom_segment(aes(xend = 0, yend = Description), color = "gray") +
  geom_point(aes(size = GeneCount, color = logFDR)) +
  scale_color_gradient(low = "blue", high = "red") +
  scale_size(range = c(3, 8)) +
  labs(
    x = "Fold Enrichment",
    y = "",
    color = "-log10(FDR)",
    size = "N. of Genes"
  ) +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 12))


#### COL8A1, COL1A1, ACTA2, TAGLN, BMP2, BMP4, ELN expression across all cells ####

# Show COL8A1 expression by cluster
col8a1_positive <- WhichCells(combined_seurat, expression = Col8a1 > 0)
length(col8a1_positive) # 1939 cells that expressed COL8A1

VlnPlot(combined_seurat, features = "Col8a1", pt.size = 0) + ggtitle("COL8A1 Expression Across All Clusters")
FeaturePlot(combined_seurat, features = "Col8a1", pt.size = 0) + ggtitle("COL8A1 Expression")

# ACTA2 - smooth muscle cell marker
acta2_positive <- WhichCells(combined_seurat, expression = Acta2 > 0)
length(acta2_positive) # 3091 cells that expressed ACTA2

VlnPlot(combined_seurat, features = "Acta2", pt.size = 0) + ggtitle("ACTA2 Expression Across All Clusters")
FeaturePlot(combined_seurat, features = "Acta2", pt.size = 0) + ggtitle("ACTA2 Expression")

# TAGLN - SMC marker
tagln_positive <- WhichCells(combined_seurat, expression = Tagln > 0)
length(tagln_positive) # 2626 cells that expressed ACTA2

VlnPlot(combined_seurat, features = "Tagln", pt.size = 0) + ggtitle("TAGLN Expression Across All Clusters")
FeaturePlot(combined_seurat, features = "Tagln", pt.size = 0) + ggtitle("TAGLN Expression")

# BMP2 and BMP4 - calcification and bone formation
VlnPlot(combined_seurat, features = "Bmp2", pt.size = 0) + ggtitle("BMP2 Expression Across All Clusters")
FeaturePlot(combined_seurat, features = "Bmp2", pt.size = 0) + ggtitle("BMP2 Expression")

VlnPlot(combined_seurat, features = "Bmp4", pt.size = 0) + ggtitle("BMP4 Expression Across All Clusters")
FeaturePlot(combined_seurat, features = "Bmp4", pt.size = 0) + ggtitle("BMP4 Expression")

# ELN - ECM protein for vascular integrity and remodeling
VlnPlot(combined_seurat, features = "Eln", pt.size = 0) + ggtitle("ELN Expression Across All Clusters")
FeaturePlot(combined_seurat, features = "Eln", pt.size = 0) + ggtitle("ELN Expression")

# COL1A1 - expressed in calcified tissues
VlnPlot(combined_seurat, features = "Col1a1", pt.size = 0) + ggtitle("COL1A1 Expression Across All Clusters")
FeaturePlot(combined_seurat, features = "Col1a1", pt.size = 0) + ggtitle("COL1A1 Expression")


#### IDENTIFYING SMCs ####
acta2_tagln_positive_cells <- intersect(acta2_positive, tagln_positive)
length(acta2_tagln_positive_cells) # 2054 cells

combined_seurat$ACTA2_TAGLN_Positive <- ifelse(Cells(combined_seurat) %in% acta2_tagln_positive_cells, "ACTA2+TAGLN+", "Other")
smc_umap <- DimPlot(combined_seurat, group.by = "ACTA2_TAGLN_Positive", label = FALSE, cols = c("ACTA2+TAGLN+" = "#03bfc4", "Other" = "#f7766d")) +
  ggtitle("UMAP of ACTA2 and TAGLN Positive Cells")
smc_umap
table(Idents(combined_seurat), combined_seurat$ACTA2_TAGLN_Positive)


#### SMC SUBCLUSTER ANALYSIS ####
SMC_subset <- subset(combined_seurat, idents = c(1,4,14))

smc_counts_lca_2d <- GetAssayData(SMC_subset, layer = "counts.lca_2d")
smc_counts_lca_2wk <- GetAssayData(SMC_subset, layer = "counts.lca_2wk")
smc_counts_rca_2d <- GetAssayData(SMC_subset, layer = "counts.rca_2d")
smc_counts_rca_2wk <- GetAssayData(SMC_subset, layer = "counts.rca_2wk")

SMC_subset_counts <- cbind(smc_counts_lca_2d, smc_counts_lca_2wk, smc_counts_rca_2d, smc_counts_rca_2wk)

# Extract the metadata for the subset cells
SMC_subset_metadata <- SMC_subset@meta.data

# Create a new Seurat object
SMC_seurat_subset <- CreateSeuratObject(
  counts = SMC_subset_counts,
  meta.data = SMC_subset_metadata,
  project = "SMC_Subset_UMAP"
)

# Normalize, find variable features, scale, run PCA
SMC_seurat_subset <- NormalizeData(SMC_seurat_subset)
SMC_seurat_subset <- FindVariableFeatures(SMC_seurat_subset, selection.method = "vst", nfeatures = 2000)
SMC_seurat_subset <- ScaleData(SMC_seurat_subset, features = rownames(SMC_seurat_subset))
SMC_seurat_subset <- RunPCA(SMC_seurat_subset, features = VariableFeatures(object = SMC_seurat_subset))

# Find clusters
SMC_seurat_subset <- FindNeighbors(SMC_seurat_subset, dims = 1:10)
SMC_seurat_subset <- FindClusters(SMC_seurat_subset, resolution = 0.5)

# Run UMAP
SMC_seurat_subset <- RunUMAP(SMC_seurat_subset, dims = 1:10)

SMC_seurat_subset_umap <- DimPlot(SMC_seurat_subset, reduction = "umap", label = TRUE) + ggtitle("UMAP of SMCs Subset")
SMC_seurat_subset_umap

SMC_subset_condition_UMAP <- DimPlot(SMC_seurat_subset, reduction = "umap", group.by = "condition", label = TRUE) + ggtitle("UMAP of SMCs by condition")
SMC_subset_condition_UMAP

# Merging clusters?

# Look at COL8A1 expression in SMCs
FeaturePlot(SMC_seurat_subset, features = "Col8a1", pt.size = 0) + ggtitle("COL8A1 Expression in SMC Clusters")

# Violin plots for COL8A1 expression across all "SMC" cell clusters
VlnPlot(SMC_seurat_subset, features = "Col8a1", pt.size = 0) + ggtitle("COL8A1 Expression Across TAGLN/ACTA2-Positive Clusters")


#### IDENTIFYING ENDOTHELIAL CELLS ####
threshold <- 0
features <- rownames(combined_seurat)
"Pecam1" %in% features
"Cdh5" %in% features

# Identify cells expressing PECAM1
pecam1_positive_cells <- WhichCells(combined_seurat, expression = Pecam1 > threshold)

# Identify cells expressing CDH5
cdh5_positive_cells <- WhichCells(combined_seurat, expression = Cdh5 > threshold)

# Find the intersection of PECAM1 and CDH5 positive cells to find which cells 
# express both
pecam1_cdh5_positive_cells <- intersect(pecam1_positive_cells, cdh5_positive_cells)
length(pecam1_cdh5_positive_cells)
ncol(combined_seurat)
# 2935 assumed endothelial cells out of 6375 cells --> about 46%

# Create a new metadata column to label PECAM1/CDH5 positive cells
combined_seurat$PECAM1_CDH5_Positive <- ifelse(Cells(combined_seurat) %in% pecam1_cdh5_positive_cells, "PECAM1+CDH5+", "Other")

# Plot UMAP and color cells based on PECAM1/CDH5 positivity
pc_umap <- DimPlot(combined_seurat, group.by = "PECAM1_CDH5_Positive", label = TRUE) +
  ggtitle("UMAP of PECAM1 and CDH5 Positive Cells")
pc_umap
table(Idents(combined_seurat), combined_seurat$PECAM1_CDH5_Positive)
# clusters with PECAM1+CDH5 are in clusters 0, 1, 5, 8, 11 (a little less), 14

#### ENDOTHELIAL SUBCLUSTER ANALYSIS ####
seurat_subset <- subset(combined_seurat, idents = c(0, 1, 5, 8, 11, 14))

# Extract counts for each dataset-specific layer
counts_lca_2d <- GetAssayData(seurat_subset, layer = "counts.lca_2d")
counts_lca_2wk <- GetAssayData(seurat_subset, layer = "counts.lca_2wk")
counts_rca_2d <- GetAssayData(seurat_subset, layer = "counts.rca_2d")
counts_rca_2wk <- GetAssayData(seurat_subset, layer = "counts.rca_2wk")

subset_counts <- cbind(counts_lca_2d, counts_lca_2wk, counts_rca_2d, counts_rca_2wk)

# Extract the metadata for the subset cells
subset_metadata <- seurat_subset@meta.data

# Create a new Seurat object
ECs_seurat_subset <- CreateSeuratObject(
  counts = subset_counts,
  meta.data = subset_metadata,
  project = "Subset_UMAP"
)

# Normalize, find variable features, scale, run PCA
ECs_seurat_subset <- NormalizeData(ECs_seurat_subset)
ECs_seurat_subset <- FindVariableFeatures(ECs_seurat_subset, selection.method = "vst", nfeatures = 2000)
ECs_seurat_subset <- ScaleData(ECs_seurat_subset, features = rownames(ECs_seurat_subset))
ECs_seurat_subset <- RunPCA(ECs_seurat_subset, features = VariableFeatures(object = ECs_seurat_subset))

# Find clusters
ECs_seurat_subset <- FindNeighbors(ECs_seurat_subset, dims = 1:10)
ECs_seurat_subset <- FindClusters(ECs_seurat_subset, resolution = 0.5)

# Run UMAP
set.seed(42)
ECs_seurat_subset <- RunUMAP(ECs_seurat_subset, dims = 1:10)

# Plot UMAP
ECs_seurat_subset_umap <- DimPlot(ECs_seurat_subset, reduction = "umap", label = TRUE) + 
                                  ggtitle("UMAP of ECs Subset")
ECs_seurat_subset_umap

# Feature plot for COL8A1
FeaturePlot(ECs_seurat_subset, features = "Col8a1", pt.size = 0) + ggtitle("COL8A1 Expression")

# Violin plots for COL8A1 expression across all "endothelial" cell clusters
VlnPlot(ECs_seurat_subset, features = "Col8a1", pt.size = 0) + ggtitle("COL8A1 Expression Across PECAM1/CDH5-Positive Clusters")

# Merging clusters
ECs_seurat_subset$merged_clusters <- ECs_seurat_subset$seurat_clusters
ECs_seurat_subset$merged_clusters[ECs_seurat_subset$seurat_clusters %in% c("2","4")] <- "2"
ECs_seurat_subset$merged_clusters[ECs_seurat_subset$seurat_clusters %in% c("1","6","7")] <- "1"
ECs_seurat_subset$merged_clusters[ECs_seurat_subset$seurat_clusters %in% c("5")] <- "4"

# Merged EC clusters UMAP

# Make the cluster labels more legible

# Fetch coordinates of the cluster labels
coords <- FetchData(ECs_seurat_subset, vars = c("umap_1", "umap_2", "merged_clusters")) %>%
  group_by(merged_clusters) %>%
  summarize(umap_1 = mean(umap_1), umap_2 = mean(umap_2))

DimPlot(ECs_seurat_subset, group.by = "merged_clusters", label = FALSE) + 
        ggtitle(NULL) + 
        geom_text(
          data = coords,
          aes(x = umap_1, y = umap_2, label = merged_clusters),
          size = 8,
          fontface = "bold",
          color = "black",
          stroke = 0.5 # works better with ggrepel sometimes
        )

# UMAP of EC clusters color coded by condition
ECs_subset_condition_UMAP <- DimPlot(ECs_seurat_subset, group.by = "condition", label = FALSE) +
                            ggtitle("EC Clusters By Condition UMAP")
ECs_subset_condition_UMAP

# UMAP of EC clusters color coded by sample
ECs_seurat_subset$sample <- paste(ECs_seurat_subset$condition, ECs_seurat_subset$time, sep = "_")

ECs_subset_samples_UMAP <- DimPlot(ECs_seurat_subset, group.by = "sample", label = FALSE) +
  ggtitle("EC Clusters By Sample UMAP")
ECs_subset_samples_UMAP

# COL8A1 expression across merged clusters
FeaturePlot(ECs_seurat_subset, features = "Col8a1", pt.size = 0) + ggtitle("COL8A1 Expression - Merged Clusters")
Idents(ECs_seurat_subset) <- "merged_clusters"
VlnPlot(ECs_seurat_subset, features = "Col8a1", pt.size = 0) +
        ggtitle(NULL) +
        theme(plot.title = element_text(size = 12))

custom_colors <- c("disturbed" = "#E63940", "normal" = "gray")
VlnPlot(ECs_seurat_subset, features = "Col8a1", group.by = "merged_clusters", split.by = "condition") +
  ggtitle("COL8A1 Expression Among Merged EC Clusters by Condition") +
  theme_minimal() +  scale_fill_manual(values = custom_colors)


#### FINDING MARKERS FOR EC SUBCLUSTERS ####
Idents(ECs_seurat_subset) <- "merged_clusters"
table(Idents(ECs_seurat_subset))

# Find all markers for each cluster
ec_markers <- FindAllMarkers(
  ECs_seurat_subset, 
  only.pos = TRUE,           # only keep positive markers
  min.pct = 0.25,            # gene must be expressed in at least 25% of cells
  logfc.threshold = 0.25     # default log fold-change threshold
)

# Top 10 markers for each cluster
top10 <- ec_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
write.csv(top10, "top10_EC_subcluster_markers.csv")

DoHeatmap(ECs_seurat_subset, features = top10$gene) + NoLegend()


# Violin plot for COL8A1 expression - disturbed/normal flow 2 days vs 2 weeks

# Set the order of the violin plots (normal-days, disturbed days, normal weeks, disturbed weeks)
ECs_seurat_subset$sample <- factor(
  ECs_seurat_subset$sample,
  levels = c("normal_days", "disturbed_days", "normal_weeks", "disturbed_weeks")
)

# Flip legend order by changing factor levels
ECs_seurat_subset$condition <- factor(
  ECs_seurat_subset$condition,
  levels = c("disturbed", "normal")
)

VlnPlot(ECs_seurat_subset, features = "Col8a1", group.by = "sample", split.by = "condition") + 
  theme_minimal() + 
  ggtitle(NULL) +
  guides(fill = guide_legend(reverse = TRUE))

# Violin plot with statistical bars
data <- FetchData(ECs_seurat_subset, vars = c("Col8a1", "sample", "condition"))
data$sample <- factor(data$sample, levels = c("normal_days", "disturbed_days", 
                                                    "normal_weeks", "disturbed_weeks"))

# Convert previous violin plot to ggplot to add stat comparison
p <- ggplot(data, aes(x = sample, y = Col8a1, fill = condition)) +
  geom_violin(position = position_dodge(0.9), trim = FALSE) +
  stat_compare_means(
    aes(group = condition),
    comparisons = comparisons,
    method = "wilcox.test",
    label = "p.signif",
    tip.length = 0.01,
    position = position_dodge(0.9)
  ) +
  coord_cartesian(ylim = c(0.25, NA)) +  # Hide expression values below 0
  theme_minimal() +
  guides(fill = guide_legend(reverse = TRUE))
p

# Expression level difference is statistically significant between conditions?
ECs_seurat_2days <- subset(ECs_seurat_subset, subset = time == "days")

# Perform differential expression analysis comparing disturbed vs normal conditions for 2 days
col8a1_diff_expression_2days <- FindMarkers(ECs_seurat_2days, features = "Col8a1", 
                                            group.by = "condition", ident.1 = "disturbed", ident.2 = "normal")

ECs_seurat_2weeks <- subset(ECs_seurat_subset, subset = time == "weeks")
col8a1_diff_expression_2weeks <- FindMarkers(ECs_seurat_2weeks, features = "Col8a1", 
                                            group.by = "condition", ident.1 = "disturbed", ident.2 = "normal")

custom_colors <- c("disturbed" = "#E63940", "normal" = "gray")
VlnPlot(ECs_seurat_subset, features = "Col8a1", group.by = "seurat_clusters", split.by = "condition") +
  ggtitle("COL8A1 Expression Among Clusters by Condition") +
  theme_minimal() +  scale_fill_manual(values = custom_colors)

#### GETTING PERCENTAGES ####

# Count number of endothelial cells in total
cell_counts <- ECs_seurat_subset@meta.data %>%
  group_by(condition) %>%
  summarise(cell_count = n())

# Count number of ECs that express COL8A1 above threshold of 1
COL8A1_expression_by_cell <- FetchData(ECs_seurat_subset, vars = "Col8a1")
COL8A1_expression_by_cell

COL8A1_expression_by_cell$condition <- ECs_seurat_subset@meta.data$condition
COL8A1_expression_by_cell$time <- ECs_seurat_subset@meta.data$time

COL8A1_positive_by_condition <- COL8A1_expression_by_cell %>%
  group_by(condition) %>%
  summarise(COL8A1_positive = sum(Col8a1 > 1))

COL8A1_positive_by_time <- COL8A1_expression_by_cell %>%
  group_by(time) %>%
  summarise(disturbed_COL8A1_positive = sum(Col8a1 > 0 & condition == "disturbed"), 
            normal_COL8A1_positive = sum(Col8a1 > 0 & condition == "normal"))

percentages <- cbind(COL8A1_positive_by_condition, cell_counts)[, c(1,2,4)]
percentages$percent <- round(percentages$COL8A1_positive/percentages$cell_count, 3)


#### GO ANALYSIS - Disturbed vs. normal EC clusters ####


# Get DEGs and filter them by significance
# Group clusters 0 and 3 as "normal" and clusters 1, 2, 4 as "disturbed"
ECs_seurat_subset$flow_type <- ifelse(ECs_seurat_subset$seurat_clusters %in% c(0, 3), "normal", "disturbed")
Idents(ECs_seurat_subset) <- "flow_type"
ECs_subset_degs <- FindMarkers(object = ECs_seurat_subset, ident.1 = "disturbed",
                              ident.2 = "normal", min.pct = 0.25, logfc.threshold = 0.25)
filtered_ECs_degs <- ECs_subset_degs %>% filter(p_val_adj < 0.05 & avg_log2FC > 0.25)
ECs_gene_list <- filtered_ECs_degs %>%
  rownames_to_column("gene") %>%
  dplyr::select(gene, avg_log2FC) %>%
  arrange(desc(avg_log2FC))

# Save DEGs to csv
write.csv(filtered_ECs_degs, "ECs_disturbed_vs_normal_degs.csv")

# Convert gene symbols to Entrez IDs
ECs_gene_list$entrez <- mapIds(org.Mm.eg.db, keys = ECs_gene_list$gene, 
                               column = "ENTREZID", keytype = "SYMBOL", 
                               multiVals = "first")

# Remove genes that converted to NA
ECs_gene_list <- ECs_gene_list[!is.na(ECs_gene_list$entrez), ]

# Run enrichment analysis (Gene Ontology - Biological Processes)
ECs_enrichment_results <- enrichGO(
  gene          = na.omit(ECs_gene_list$entrez),
  OrgDb         = org.Mm.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  pAdjustMethod = "BH",
  qvalueCutoff  = 0.05
)

length(na.omit(ECs_gene_list$entrez))

# Calculate fold enrichment, generate a bar plot
ECs_enrichment_data <- as.data.frame(ECs_enrichment_results) %>%
  mutate(
    FoldEnrichment = calculate_fold_enrichment(GeneRatio, BgRatio),
    logFDR = -log10(p.adjust),
    GeneCount = as.numeric(sub("/.*", "", GeneRatio))
  ) %>%
  arrange(desc(logFDR)) %>%
  head(20)

ggplot(ECs_enrichment_data, aes(x = FoldEnrichment, y = fct_reorder(Description, FoldEnrichment))) +
  geom_segment(aes(xend = 0, yend = Description), color = "gray") +
  geom_point(aes(size = GeneCount, color = logFDR)) +
  scale_color_gradient(low = "blue", high = "red") +
  scale_size(range = c(3, 8)) +
  labs(
    x = "Fold Enrichment",
    y = "",
    color = "-log10(FDR)",
    size = "N. of Genes"
  ) +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 12))

# Run GO-MF analysis
ECs_MF_results <- enrichGO(
  gene          = na.omit(ECs_gene_list$entrez),
  OrgDb         = org.Mm.eg.db,
  keyType       = "ENTREZID",
  ont           = "MF",
  pAdjustMethod = "BH",
  qvalueCutoff  = 0.05
)

# Calculate fold enrichment, generate a bar plot
ECs_MF_data <- as.data.frame(ECs_MF_results) %>%
  mutate(
    FoldEnrichment = calculate_fold_enrichment(GeneRatio, BgRatio),
    logFDR = -log10(p.adjust),
    GeneCount = as.numeric(sub("/.*", "", GeneRatio))
  ) %>%
  arrange(desc(logFDR)) %>%
  head(20)

ggplot(ECs_MF_data, aes(x = FoldEnrichment, y = fct_reorder(Description, FoldEnrichment))) +
  geom_segment(aes(xend = 0, yend = Description), color = "gray") +
  geom_point(aes(size = GeneCount, color = logFDR)) +
  scale_color_gradient(low = "blue", high = "red") +
  scale_size(range = c(3, 8)) +
  labs(
    x = "Fold Enrichment",
    y = "",
    color = "-log10(FDR)",
    size = "N. of Genes"
  ) +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 12))

# Run GO-CC analysis
ECs_CC_results <- enrichGO(
  gene          = na.omit(ECs_gene_list$entrez),
  OrgDb         = org.Mm.eg.db,
  keyType       = "ENTREZID",
  ont           = "CC",
  pAdjustMethod = "BH",
  qvalueCutoff  = 0.05
)

# Calculate fold enrichment, generate a bar plot
ECs_CC_data <- as.data.frame(ECs_CC_results) %>%
  mutate(
    FoldEnrichment = calculate_fold_enrichment(GeneRatio, BgRatio),
    logFDR = -log10(p.adjust),
    GeneCount = as.numeric(sub("/.*", "", GeneRatio))
  ) %>%
  arrange(desc(logFDR)) %>%
  head(20)

ggplot(ECs_CC_data, aes(x = FoldEnrichment, y = fct_reorder(Description, FoldEnrichment))) +
  geom_segment(aes(xend = 0, yend = Description), color = "gray") +
  geom_point(aes(size = GeneCount, color = logFDR)) +
  scale_color_gradient(low = "blue", high = "red") +
  scale_size(range = c(3, 8)) +
  labs(
    x = "Fold Enrichment",
    y = "",
    color = "-log10(FDR)",
    size = "N. of Genes"
  ) +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 12))



#### GO ANALYSIS - Long-term disturbed (cluster 1) vs. rest of EC subset ####

Idents(ECs_seurat_subset) <- "seurat_clusters"
ECs_degs_1 <- FindMarkers(object = ECs_seurat_subset, ident.1 = "1",
                          min.pct = 0.25, logfc.threshold = 0.25)
filtered_ECs_degs_1 <- ECs_degs_1 %>% filter(p_val_adj < 0.05 & avg_log2FC > 0.25)
ECs_gene_list_1 <- filtered_ECs_degs_1 %>%
  rownames_to_column("gene") %>%
  dplyr::select(gene, avg_log2FC) %>%
  arrange(desc(avg_log2FC))

# save DEGs to csv
write.csv(filtered_ECs_degs_1, "cluster1_vs_ECs_DEGs.csv")

# Convert gene symbols to Entrez IDs
ECs_gene_list_1$entrez <- mapIds(org.Mm.eg.db, keys = ECs_gene_list_1$gene, 
                               column = "ENTREZID", keytype = "SYMBOL", 
                               multiVals = "first")

# Remove genes that converted to NA
ECs_gene_list_1 <- ECs_gene_list_1[!is.na(ECs_gene_list_1$entrez), ]

# Run enrichment analysis (Gene Ontology - Biological Processes)
ECs_enrichment_res_1 <- enrichGO(
  gene          = na.omit(ECs_gene_list_1$entrez),
  OrgDb         = org.Mm.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  pAdjustMethod = "BH",
  qvalueCutoff  = 0.05
)

# Generate a bar plot 
ECs_enrichment_data_1 <- as.data.frame(ECs_enrichment_res_1) %>%
  mutate(
    FoldEnrichment = calculate_fold_enrichment(GeneRatio, BgRatio),
    logFDR = -log10(p.adjust),
    GeneCount = as.numeric(sub("/.*", "", GeneRatio))
  ) %>%
  arrange(desc(logFDR)) %>%
  head(20)

ggplot(ECs_enrichment_data_1, aes(x = FoldEnrichment, y = fct_reorder(Description, FoldEnrichment))) +
  geom_segment(aes(xend = 0, yend = Description), color = "gray") +
  geom_point(aes(size = GeneCount, color = logFDR)) +
  scale_color_gradient(low = "blue", high = "red") +
  scale_size(range = c(3, 8)) +
  labs(
    x = "Fold Enrichment",
    y = "",
    color = "-log10(FDR)",
    size = "N. of Genes"
  ) +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 12))

# Run GO-MF analysis
ECs_1_MF_res <- enrichGO(
  gene          = na.omit(ECs_gene_list_1$entrez),
  OrgDb         = org.Mm.eg.db,
  keyType       = "ENTREZID",
  ont           = "MF",
  pAdjustMethod = "BH",
  qvalueCutoff  = 0.05
)

# Generate a bar plot 
ECs_1_MF_data <- as.data.frame(ECs_1_MF_res) %>%
  mutate(
    FoldEnrichment = calculate_fold_enrichment(GeneRatio, BgRatio),
    logFDR = -log10(p.adjust),
    GeneCount = as.numeric(sub("/.*", "", GeneRatio))
  ) %>%
  arrange(desc(logFDR)) %>%
  head(20)

ggplot(ECs_1_MF_data, aes(x = FoldEnrichment, y = fct_reorder(Description, FoldEnrichment))) +
  geom_segment(aes(xend = 0, yend = Description), color = "gray") +
  geom_point(aes(size = GeneCount, color = logFDR)) +
  scale_color_gradient(low = "blue", high = "red") +
  scale_size(range = c(3, 8)) +
  labs(
    x = "Fold Enrichment",
    y = "",
    color = "-log10(FDR)",
    size = "N. of Genes"
  ) +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 12))

# Run GO-CC analysis
ECs_1_CC_res <- enrichGO(
  gene          = na.omit(ECs_gene_list_1$entrez),
  OrgDb         = org.Mm.eg.db,
  keyType       = "ENTREZID",
  ont           = "CC",
  pAdjustMethod = "BH",
  qvalueCutoff  = 0.05
)

# Generate a bar plot 
ECs_1_CC_data <- as.data.frame(ECs_1_CC_res) %>%
  mutate(
    FoldEnrichment = calculate_fold_enrichment(GeneRatio, BgRatio),
    logFDR = -log10(p.adjust),
    GeneCount = as.numeric(sub("/.*", "", GeneRatio))
  ) %>%
  arrange(desc(logFDR)) %>%
  head(20)

ggplot(ECs_1_CC_data, aes(x = FoldEnrichment, y = fct_reorder(Description, FoldEnrichment))) +
  geom_segment(aes(xend = 0, yend = Description), color = "gray") +
  geom_point(aes(size = GeneCount, color = logFDR)) +
  scale_color_gradient(low = "blue", high = "red") +
  scale_size(range = c(3, 8)) +
  labs(
    x = "Fold Enrichment",
    y = "",
    color = "-log10(FDR)",
    size = "N. of Genes"
  ) +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 12))

#### GO ANALYSIS - Week vs. days clusters for normal flow cells ####

normal_ECs_seurat_subset <- subset(ECs_seurat_subset, subset = flow_type == "normal")
table(normal_ECs_seurat_subset$flow_type) # 1246 endothelial cells exposed to normal flow

# Find DEGs between weeks and days groups
Idents(normal_ECs_seurat_subset) <- normal_ECs_seurat_subset$time
norm_ECs_subset_degs <- FindMarkers(object = normal_ECs_seurat_subset, ident.1 = "weeks",
                               ident.2 = "days", min.pct = 0.25, logfc.threshold = 0.25)
filtered_norm_ECs_degs <- norm_ECs_subset_degs %>% filter(p_val_adj < 0.05 & avg_log2FC > 0.25)
norm_ECs_gene_list <- filtered_norm_ECs_degs %>%
  rownames_to_column("gene") %>%
  dplyr::select(gene, avg_log2FC) %>%
  arrange(desc(avg_log2FC))

# save DEGs to csv
write.csv(filtered_norm_ECs_degs, "normalflow_ECs_2week_vs_2day.csv")

# Convert gene symbols to Entrez IDs
norm_ECs_gene_list$entrez <- mapIds(org.Mm.eg.db, keys = norm_ECs_gene_list$gene, 
                               column = "ENTREZID", keytype = "SYMBOL", 
                               multiVals = "first")

# Remove genes that converted to NA
norm_ECs_gene_list <- norm_ECs_gene_list[!is.na(norm_ECs_gene_list$entrez), ]

# Run enrichment analysis (Gene Ontology - Biological Processes)
norm_ECs_enrichment_results <- enrichGO(
  gene          = na.omit(norm_ECs_gene_list$entrez),
  OrgDb         = org.Mm.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  pAdjustMethod = "BH",
  qvalueCutoff  = 0.05
)

norm_ECs_enrichment_data <- as.data.frame(norm_ECs_enrichment_results) %>%
  mutate(
    FoldEnrichment = calculate_fold_enrichment(GeneRatio, BgRatio),
    logFDR = -log10(p.adjust),
    GeneCount = as.numeric(sub("/.*", "", GeneRatio))
  ) %>%
  arrange(desc(logFDR)) %>%
  head(20)

ggplot(norm_ECs_enrichment_data, aes(x = FoldEnrichment, y = fct_reorder(Description, FoldEnrichment))) +
  geom_segment(aes(xend = 0, yend = Description), color = "gray") +
  geom_point(aes(size = GeneCount, color = logFDR)) +
  scale_color_gradient(low = "blue", high = "red") +
  scale_size(range = c(3, 8)) +
  labs(
    x = "Fold Enrichment",
    y = "",
    color = "-log10(FDR)",
    size = "N. of Genes"
  ) +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 12))


# Run GO-MF analysis
norm_ECs_MF_results <- enrichGO(
  gene          = na.omit(norm_ECs_gene_list$entrez),
  OrgDb         = org.Mm.eg.db,
  keyType       = "ENTREZID",
  ont           = "MF",
  pAdjustMethod = "BH",
  qvalueCutoff  = 0.05
)

norm_ECs_MF_data <- as.data.frame(norm_ECs_MF_results) %>%
  mutate(
    FoldEnrichment = calculate_fold_enrichment(GeneRatio, BgRatio),
    logFDR = -log10(p.adjust),
    GeneCount = as.numeric(sub("/.*", "", GeneRatio))
  ) %>%
  arrange(desc(logFDR)) %>%
  head(20)

ggplot(norm_ECs_MF_data, aes(x = FoldEnrichment, y = fct_reorder(Description, FoldEnrichment))) +
  geom_segment(aes(xend = 0, yend = Description), color = "gray") +
  geom_point(aes(size = GeneCount, color = logFDR)) +
  scale_color_gradient(low = "blue", high = "red") +
  scale_size(range = c(3, 8)) +
  labs(
    x = "Fold Enrichment",
    y = "",
    color = "-log10(FDR)",
    size = "N. of Genes"
  ) +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 12))

# Run GO-CC analysis
norm_ECs_CC_results <- enrichGO(
  gene          = na.omit(norm_ECs_gene_list$entrez),
  OrgDb         = org.Mm.eg.db,
  keyType       = "ENTREZID",
  ont           = "CC",
  pAdjustMethod = "BH",
  qvalueCutoff  = 0.05
)

norm_ECs_CC_data <- as.data.frame(norm_ECs_CC_results) %>%
  mutate(
    FoldEnrichment = calculate_fold_enrichment(GeneRatio, BgRatio),
    logFDR = -log10(p.adjust),
    GeneCount = as.numeric(sub("/.*", "", GeneRatio))
  ) %>%
  arrange(desc(logFDR)) %>%
  head(20)

ggplot(norm_ECs_CC_data, aes(x = FoldEnrichment, y = fct_reorder(Description, FoldEnrichment))) +
  geom_segment(aes(xend = 0, yend = Description), color = "gray") +
  geom_point(aes(size = GeneCount, color = logFDR)) +
  scale_color_gradient(low = "blue", high = "red") +
  scale_size(range = c(3, 8)) +
  labs(
    x = "Fold Enrichment",
    y = "",
    color = "-log10(FDR)",
    size = "N. of Genes"
  ) +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 12))

#### GO ANALYSIS - Week vs. days clusters for disturbed flow cells ####
disturbed_ECs_seurat_subset <- subset(ECs_seurat_subset, subset = flow_type == "disturbed")
table(disturbed_ECs_seurat_subset$flow_type) # 1597 endothelial cells exposed to disturbed flow

# Find DEGs between weeks and days groups
Idents(disturbed_ECs_seurat_subset) <- disturbed_ECs_seurat_subset$time
dist_ECs_subset_degs <- FindMarkers(object = disturbed_ECs_seurat_subset, ident.1 = "weeks",
                                    ident.2 = "days", min.pct = 0.25, logfc.threshold = 0.25)
filtered_dist_ECs_degs <- dist_ECs_subset_degs %>% filter(p_val_adj < 0.05 & avg_log2FC > 0.25)
dist_ECs_gene_list <- filtered_dist_ECs_degs %>%
  rownames_to_column("gene") %>%
  dplyr::select(gene, avg_log2FC) %>%
  arrange(desc(avg_log2FC))

# save DEGs to csv
write.csv(filtered_dist_ECs_degs, "disturbedECs_2week_vs_2day_DEGs.csv")

# Convert gene symbols to Entrez IDs
dist_ECs_gene_list$entrez <- mapIds(org.Mm.eg.db, keys = dist_ECs_gene_list$gene, 
                                    column = "ENTREZID", keytype = "SYMBOL", 
                                    multiVals = "first")

# Run enrichment analysis (Gene Ontology - Biological Processes)
dist_ECs_enrichment_results <- enrichGO(
  gene          = na.omit(dist_ECs_gene_list$entrez),
  OrgDb         = org.Mm.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  pAdjustMethod = "BH",
  qvalueCutoff  = 0.05
)

dist_ECs_enrichment_data <- as.data.frame(dist_ECs_enrichment_results) %>%
  mutate(
    FoldEnrichment = calculate_fold_enrichment(GeneRatio, BgRatio),
    logFDR = -log10(p.adjust),
    GeneCount = as.numeric(sub("/.*", "", GeneRatio))
  ) %>%
  arrange(desc(logFDR)) %>%
  head(20)

ggplot(dist_ECs_enrichment_data, aes(x = FoldEnrichment, y = fct_reorder(Description, FoldEnrichment))) +
  geom_segment(aes(xend = 0, yend = Description), color = "gray") +
  geom_point(aes(size = GeneCount, color = logFDR)) +
  scale_color_gradient(low = "blue", high = "red") +
  scale_size(range = c(3, 8)) +
  labs(
    x = "Fold Enrichment",
    y = "",
    color = "-log10(FDR)",
    size = "N. of Genes"
  ) +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 12))

# Run GO-MF analysis
dist_ECs_MF_results <- enrichGO(
  gene          = na.omit(dist_ECs_gene_list$entrez),
  OrgDb         = org.Mm.eg.db,
  keyType       = "ENTREZID",
  ont           = "MF",
  pAdjustMethod = "BH",
  qvalueCutoff  = 0.05
)

dist_ECs_MF_data <- as.data.frame(dist_ECs_MF_results) %>%
  mutate(
    FoldEnrichment = calculate_fold_enrichment(GeneRatio, BgRatio),
    logFDR = -log10(p.adjust),
    GeneCount = as.numeric(sub("/.*", "", GeneRatio))
  ) %>%
  arrange(desc(logFDR)) %>%
  head(15)

dist_ECs_MF_data <- dist_ECs_MF_data %>%
  mutate(Category = case_when(
    grepl("MHC|CD8|antigen", Description, ignore.case = TRUE) ~ "MHC/antigen binding",
    grepl("RNA polymerase|transcription|DNA-binding|SMAD|receptor", Description, ignore.case = TRUE) ~ "Transcription regulation",
    grepl("GTP|G protein|GTPase|guanyl", Description, ignore.case = TRUE) ~ "GTP/G-protein activity",
    grepl("ubiquitin", Description, ignore.case = TRUE) ~ "Ubiquitin ligase activity",
    grepl("collagen|extracellular matrix", Description, ignore.case = TRUE) ~ "Extracellular matrix",
    grepl("peptide|amide|heterodimer", Description, ignore.case = TRUE) ~ "Peptide/rotein binding",
    grepl("growth factor", Description, ignore.case = TRUE) ~ "Growth factor binding"
  ))

collapsed_data_MF <- dist_ECs_MF_data %>%
  group_by(Category) %>%
  summarise(
    FoldEnrichment = mean(FoldEnrichment),
    logFDR = mean(logFDR),
    GeneCount = sum(GeneCount)
  ) %>%
  ungroup() %>%
  arrange(desc(logFDR))


ggplot(collapsed_data_MF, aes(x = FoldEnrichment, y = fct_reorder(Category, FoldEnrichment))) +
  geom_segment(aes(xend = 0, yend = Category), color = "gray") +
  geom_point(aes(size = GeneCount, color = logFDR)) +
  scale_color_gradient(low = "blue", high = "red") +
  scale_size(range = c(3, 8)) +
  labs(
    x = "Fold Enrichment",
    y = "",
    color = "-log10(FDR)",
    size = "N. of Genes"
  ) +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 12))


# Run GO-CC analysis
dist_ECs_CC_results <- enrichGO(
  gene          = na.omit(dist_ECs_gene_list$entrez),
  OrgDb         = org.Mm.eg.db,
  keyType       = "ENTREZID",
  ont           = "CC",
  pAdjustMethod = "BH",
  qvalueCutoff  = 0.05
)

dist_ECs_CC_data <- as.data.frame(dist_ECs_CC_results) %>%
  mutate(
    FoldEnrichment = calculate_fold_enrichment(GeneRatio, BgRatio),
    logFDR = -log10(p.adjust),
    GeneCount = as.numeric(sub("/.*", "", GeneRatio))
  ) %>%
  arrange(desc(logFDR)) %>%
  head(15)

# remove actomyosin
dist_ECs_CC_data <- dist_ECs_CC_data %>% filter(Description != "actomyosin")

# consolidate GO terms
dist_ECs_CC_data <- dist_ECs_CC_data %>%
  mutate(Category = case_when(
    grepl("MHC", Description, ignore.case = TRUE) ~ "MHC complex",
    grepl("actin", Description, ignore.case = TRUE) ~ "Actin filament",
    grepl("stress fiber", Description, ignore.case = TRUE) ~ "Actin filament",
    grepl("Golgi", Description, ignore.case = TRUE) ~ "Golgi membrane",
    grepl("endoplasmic reticulum", Description, ignore.case = TRUE) ~ "Endoplasmic reticulum",
    grepl("extracellular matrix", Description, ignore.case = TRUE) ~ "Extracellular matrix",
    grepl("basement", Description, ignore.case = TRUE) ~ "Extracellular matrix",
    grepl("vacuo", Description, ignore.case = TRUE) ~ "Vacuolar membrane",
    TRUE ~  "RNA polymerase II transcription regulator complex"
  ))

collapsed_data_CC <- dist_ECs_CC_data %>%
  group_by(Category) %>%
  summarise(
    FoldEnrichment = mean(FoldEnrichment),
    logFDR = mean(logFDR),
    GeneCount = sum(GeneCount)
  ) %>%
  ungroup() %>%
  arrange(desc(logFDR))

ggplot(collapsed_data_CC, aes(x = FoldEnrichment, y = fct_reorder(Category, FoldEnrichment))) +
  geom_segment(aes(xend = 0, yend = Category), color = "gray") +
  geom_point(aes(size = GeneCount, color = logFDR)) +
  scale_color_gradient(low = "blue", high = "red") +
  scale_size(range = c(3, 8)) +
  labs(
    x = "Fold Enrichment",
    y = "",
    color = "-log10(FDR)",
    size = "N. of Genes"
  ) +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 12))


