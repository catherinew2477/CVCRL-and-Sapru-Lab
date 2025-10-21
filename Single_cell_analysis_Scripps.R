# CODE FROM KYLE

setwd("/Users/catherinewang/Desktop/CVCRL/Scripps single-cell")

# BiocManager::install("monocle3")
# install.packages("circlize")

# Install Seurat v4.3.0.1 (last v4 release before v5)
# remotes::install_version("Seurat", version = "4.3.0.1")

library(Seurat)
library(dplyr)
library(ggplot2)
#library(monocle3)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ComplexHeatmap)
library(tibble)
library(circlize)
library(tidyr)

# CREATE SEURAT OBJECTS --------------
# Load the datasets
p1_ac_data <- Read10X("GSE159677_RAW/p1_ac/filtered_feature_bc_matrix")
p1_pa_data <- Read10X("GSE159677_RAW/p1_pa/filtered_feature_bc_matrix")
p2_ac_data <- Read10X("GSE159677_RAW/p2_ac/filtered_feature_bc_matrix")
p2_pa_data <- Read10X("GSE159677_RAW/p2_pa/filtered_feature_bc_matrix")
p3_ac_data <- Read10X("GSE159677_RAW/p3_ac/filtered_feature_bc_matrix")
p3_pa_data <- Read10X("GSE159677_RAW/p3_pa/filtered_feature_bc_matrix")

# Create Seurat objects
seurat_p1_ac <- CreateSeuratObject(counts = p1_ac_data, project = "p1_ac")
seurat_p1_pa <- CreateSeuratObject(counts = p1_pa_data, project = "p1_pa")
seurat_p2_ac <- CreateSeuratObject(counts = p2_ac_data, project = "p2_ac")
seurat_p2_pa <- CreateSeuratObject(counts = p2_pa_data, project = "p2_pa")
seurat_p3_ac <- CreateSeuratObject(counts = p3_ac_data, project = "p3_ac")
seurat_p3_pa <- CreateSeuratObject(counts = p3_pa_data, project = "p3_pa")

## METADATA-----------------------
# Add metadata for disease status
seurat_p1_ac$condition <- "diseased"
seurat_p1_pa$condition <- "healthy"
seurat_p2_ac$condition <- "diseased"
seurat_p2_pa$condition <- "healthy"
seurat_p3_ac$condition <- "diseased"
seurat_p3_pa$condition <- "healthy"

# Add metadata for patient ID
seurat_p1_ac$patient <- "p1"
seurat_p1_pa$patient <- "p1"
seurat_p2_ac$patient <- "p2"
seurat_p2_pa$patient <- "p2"
seurat_p3_ac$patient <- "p3"
seurat_p3_pa$patient <- "p3"

## COMBINE THE DATASETS ---------
# Combine datasets
combined_seurat <- merge(seurat_p1_ac, y = list(seurat_p1_pa, seurat_p2_ac, seurat_p2_pa, seurat_p3_ac, seurat_p3_pa))

# Perform QC
combined_seurat[["percent.mt"]] <- PercentageFeatureSet(combined_seurat, pattern = "^MT-")
combined_seurat <- subset(combined_seurat, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 10)

# Normalize and find variable features
combined_seurat <- NormalizeData(combined_seurat)
combined_seurat <- FindVariableFeatures(combined_seurat)

# SCALING AND DIMENSIONALITY REDUCTION --------------
combined_seurat <- ScaleData(combined_seurat)
combined_seurat <- RunPCA(combined_seurat, npcs = 100)

# UMAP and clustering
combined_seurat <- RunUMAP(combined_seurat, dims = 1:100)
combined_seurat <- FindNeighbors(combined_seurat, dims = 1:100)
combined_seurat <- FindClusters(combined_seurat, resolution = 0.5)

# Plot UMAP
set.seed(123)
initial_umap <- DimPlot(combined_seurat, reduction = "umap", label = TRUE) 
initial_umap

## ENDOTHELIAL CELLS ---------
# Define a threshold for identifying positive cells (you can adjust this based on your dataset)
threshold <- 0  # This can be set depending on expression levels in your data

# Identify cells expressing PECAM1
pecam1_positive_cells <- WhichCells(combined_seurat, expression = PECAM1 > threshold)

# Identify cells expressing CDH5
cdh5_positive_cells <- WhichCells(combined_seurat, expression = CDH5 > threshold)

# Find the intersection of PECAM1 and CDH5 positive cells
pecam1_cdh5_positive_cells <- intersect(pecam1_positive_cells, cdh5_positive_cells)

# Get the number of PECAM1/CDH5 positive cells
length(pecam1_cdh5_positive_cells)

# Create a new metadata column to label PECAM1/CDH5 positive cells
combined_seurat$PECAM1_CDH5_Positive <- ifelse(Cells(combined_seurat) %in% pecam1_cdh5_positive_cells, "PECAM1+CDH5+", "Other")

# Plot UMAP and color cells based on PECAM1/CDH5 positivity
ecs_umap <- DimPlot(combined_seurat, group.by = "PECAM1_CDH5_Positive", label = TRUE) +
  ggtitle("UMAP of PECAM1 and CDH5 Positive Cells")
ecs_umap

## SUBCLUSTER ANALYSIS ------
# Subset clusters 
seurat_subset <- subset(combined_seurat, idents = c(5, 8, 12, 24))

# Normalize the data
seurat_subset <- NormalizeData(seurat_subset)

# Identify the most variable features (genes)
seurat_subset <- FindVariableFeatures(seurat_subset, selection.method = "vst", nfeatures = 2000)

# Scale the data
seurat_subset <- ScaleData(seurat_subset, features = rownames(seurat_subset))

# Run PCA
seurat_subset <- RunPCA(seurat_subset, features = VariableFeatures(object = seurat_subset))

# Find clusters
seurat_subset <- FindNeighbors(seurat_subset, dims = 1:10)
seurat_subset <- FindClusters(seurat_subset, resolution = 0.5)

# Run UMAP
set.seed(123)
seurat_subset <- RunUMAP(seurat_subset, dims = 1:10)

# Plot UMAP
subset_umap <- DimPlot(seurat_subset, reduction = "umap", label = TRUE) + ggtitle("UMAP of Endothelial Cell Subset")
print(subset_umap)


# Identify cells expressing PECAM1
pecam1_positive_cells <- WhichCells(seurat_subset, expression = PECAM1 > threshold)

# Identify cells expressing CDH5
cdh5_positive_cells <- WhichCells(seurat_subset, expression = CDH5 > threshold)

# Find the intersection of PECAM1 and CDH5 positive cells
pecam1_cdh5_positive_cells <- intersect(pecam1_positive_cells, cdh5_positive_cells)

# Get the number of PECAM1/CDH5 positive cells
length(pecam1_cdh5_positive_cells)

# Create a new metadata column to label PECAM1/CDH5 positive cells
seurat_subset$PECAM1_CDH5_Positive <- ifelse(Cells(seurat_subset) %in% pecam1_cdh5_positive_cells, "PECAM1+CDH5+", "Other")

# Plot UMAP and color cells based on PECAM1/CDH5 positivity
ecs_subset_umap <- DimPlot(seurat_subset, group.by = "PECAM1_CDH5_Positive", label = TRUE) +
  ggtitle("UMAP of PECAM1 and CDH5 Positive Cells")
ecs_subset_umap

VlnPlot(seurat_subset, features = c("BMP4"), group.by = "seurat_clusters", pt.size = 0)
FeaturePlot(seurat_subset, features = c("BMP4"))

# Calculating proportions of positive cells in each cluster in subset
meta <- seurat_subset@meta.data %>%
  select(seurat_clusters, PECAM1_CDH5_Positive)

ec_markers_per_cluster <- meta %>%
  group_by(seurat_clusters) %>%
  summarise(
    total_cells = n(),
    PECAM1_CDH5_Positive = sum(PECAM1_CDH5_Positive == "PECAM1+CDH5+"),
    proportion_positive = PECAM1_CDH5_Positive / total_cells
  )

ec_markers_per_cluster$proportion_positive <- round(ec_markers_per_cluster$proportion_positive, digits = 2)


# Calculating proportions of diseased/healthy cells that express COL8A1

# Access the raw counts
col8a1_counts <- GetAssayData(seurat_subset, slot = "counts")["COL8A1", ]

# Add a metadata column: COL8A1+
seurat_subset$COL8A1_Positive <- ifelse(col8a1_counts > 0, "COL8A1+", "COL8A1-")

col8a1_props_in_ECs <- seurat_subset@meta.data %>%
  group_by(condition) %>%
  summarise(
    total_cells = n(),
    COL8A1_Positive = sum(COL8A1_Positive == "COL8A1+"),
    proportion_positive = COL8A1_Positive / total_cells
  )

col8a1_props_in_ECs$proportion_positive <- round(col8a1_props_in_ECs$proportion_positive, 2)



# visualize COL8A1 expression across EC clusters 

FeaturePlot(seurat_subset, features = "COL8A1")
VlnPlot(seurat_subset, features = "COL8A1", group.by = "seurat_clusters")

# color cells based on belonging to AC or PA sample
DimPlot(seurat_subset, group.by = "condition", label = TRUE) +
  ggtitle("UMAP Colored by Condition (AC vs PA)")


#### GO ANALYSIS BETWEEN HIGH AND LOW COL8A1 CLUSTERS ####


# Defining high and low COL8A1 clusters
high_clusters <- c("2", "9", "10", "11") # has COl8A1 expression
low_clusters  <- setdiff(unique(Idents(seurat_subset)), high_clusters)

seurat_subset <- JoinLayers(seurat_subset, assay = "RNA")

deg_col8a1 <- FindMarkers(
  seurat_subset,
  ident.1 = high_clusters,
  ident.2 = low_clusters,
  logfc.threshold = 0.25,    # genes upregulated by at least 1.2x
  min.pct = 0.1              # expressed in at least 10% of cells
)

top10_markers <- deg_col8a1 %>%
  arrange(p_val_adj) %>%  # sort genes by adjusted p-value
  head(10) %>%            # take the first 10 rows
  rownames()  

# heatmap of top 10 markers
top10_markers_heatmap <- DoHeatmap(
  seurat_subset,
  features = top10_markers,
  group.by = "seurat_clusters"  
) +
  theme(
    axis.text.y = element_text(size = 20),  # increase gene name font size
  )
top10_markers_heatmap

deg_col8a1$gene <- rownames(deg_col8a1)
write.csv(deg_col8a1, "deg_col8a1.csv")

# Install BiocManager
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Install clusterProfiler and org.Hs.eg.db from Bioconductor
BiocManager::install("clusterProfiler")
BiocManager::install("org.Hs.eg.db")


library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(ggplot2)
library(forcats)
library(AnnotationDbi)

# Prepare gene list from FindMarkers results
# Keep only upregulated genes in high COL8A1 clusters
deg_genes <- deg_col8a1 %>%
  filter(p_val_adj < 0.05) %>%     
  rownames() 

# Map gene symbols to Entrez IDs
deg_genes_entrez <- mapIds(
  org.Hs.eg.db, 
  keys = deg_genes,
  column = "ENTREZID",
  keytype = "SYMBOL",
  multiVals = "first"
)

# Remove NAs
deg_genes_entrez <- na.omit(deg_genes_entrez)

# customize background genes list so it only includes genes expressed in EC subset

# Extract gene symbols for DEGs (upregulated in high COL8A1 clusters)
deg_genes_symbols <- rownames(deg_col8a1)

all_genes_symbols <- rownames(seurat_subset)

# Convert both DEG list and background to Entrez IDs
deg_genes_entrez <- bitr(
  deg_genes_symbols,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)$ENTREZID

all_genes_entrez <- bitr(
  all_genes_symbols,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)$ENTREZID

# Run GO analysis with customized background gene list
col8a1_GO <- enrichGO(
  gene          = deg_genes_entrez,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "CC",                 # Biological Process
  universe      = all_genes_entrez,     # custom background
  pAdjustMethod = "BH",
  qvalueCutoff  = 0.05
)

# Run GO enrichment (Biological Process) with whole genome as background oops
col8a1_GO <- enrichGO(
  gene          = deg_genes_entrez,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "CC",
  pAdjustMethod = "BH",
  qvalueCutoff  = 0.05
)

# Prepare data for plotting
calculate_fold_enrichment <- function(GeneRatio, BgRatio) {
  gene_numerator <- as.numeric(sub("/.*", "", GeneRatio))
  gene_denominator <- as.numeric(sub(".*/", "", GeneRatio))
  
  bg_numerator <- as.numeric(sub("/.*", "", BgRatio))
  bg_denominator <- as.numeric(sub(".*/", "", BgRatio))
  
  (gene_numerator / gene_denominator) / (bg_numerator / bg_denominator)
}

col8a1_GO_data <- as.data.frame(col8a1_GO) %>%
  mutate(
    FoldEnrichment = calculate_fold_enrichment(GeneRatio, BgRatio),
    logFDR = -log10(p.adjust),
    GeneCount = as.numeric(sub("/.*", "", GeneRatio))
  ) %>%
  arrange(desc(logFDR)) %>%
  head(20)

# Plot top enriched GO terms
ggplot(col8a1_GO_data, aes(x = FoldEnrichment, y = fct_reorder(Description, FoldEnrichment))) +
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
