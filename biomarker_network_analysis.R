setwd("/Users/catherinewang/Desktop/Sapru/TRAPMODS/")

#### LIBRARIES ####
library(ComplexHeatmap)
library(rstatix)
library(dplyr)
library(tibble)
library(tidyr)
library(circlize) # make sure 0 in heatmap is white
library(uwot)
library(ggplot2)
library(tidyverse)

#### GENERATING CORRELATION MATRIX/HEATMAP ####

# load CSV file of NULISA markers (log2 fold changes)
network_markers <- read.table("network_markers_cw.csv", sep = ",", header = TRUE, row.names = 1)

# Generate correlation matrix
# remove patients with NA values
filtered_markers <- network_markers[complete.cases(network_markers), ]
filtered_markers[] <- lapply(filtered_markers, as.numeric)
corr_matrix <- cor(filtered_markers, method = "spearman", use = "pairwise.complete.obs")

col_fun <- colorRamp2(
  c(-0.7, 0, 0.7), 
  c("blue", "white", "red")
)

# Generate heatmap of correlation matrix (make it prettier later)
Heatmap(corr_matrix, 
        col = col_fun,
        show_row_names = FALSE,
        show_column_names = FALSE,
        show_column_dend = TRUE,
        name = "Rho",
        column_title = "Initial correlation matrix of fold change response
        to transfusion (n = 249)"
)


# Get p-values
p_vals <- network_markers %>% cor_test(method = "spearman", )
p_vals$p_adjusted <- p.adjust(p_vals$p, method = "fdr")

# Convert cor_test() results to matrix
adj_p_vals_matrix <- p_vals %>% 
  select(var1, var2, p_adjusted) %>%
  pivot_wider(names_from = var2, values_from = p_adjusted) %>%
  column_to_rownames("var1")

# Export CSV of p-value matrix
write.csv(adj_p_vals_matrix, "adj_p_vals_matrix.csv")

# -log10 for visualization
p_vals_log_matrix <- -log10(as.matrix(adj_p_vals_matrix))

# Replace Inf (adjusted p-value = 0) with NA for visualization
p_vals_log_matrix[p_vals_log_matrix == Inf] <- NA

Heatmap(p_vals_log_matrix, 
        name = "-log10(p_adjusted)", 
        na_col = "red",
        show_row_names = FALSE,
        show_column_names = FALSE,
        show_column_dend = TRUE,
        column_title = "Log-transformed adjusted p-value matrix"
        )

#### CORRELATION MATRIX HEATMAP ####
rows <- hclust(dist(corr_matrix))
cols <- hclust(dist(corr_matrix))

# Cut dendrograms into clusters

# determine k with elbow plot
library(factoextra)
fviz_nbclust(as.matrix(corr_matrix), FUN = hcut, method = "wss")

k <- 4
row_clusters <- cutree(rows, k = k)
col_clusters <- cutree(cols, k = k)

# Create a data frame mapping each gene to its cluster
marker_names <- rownames(corr_matrix)
row_cluster_df <- data.frame(Gene = marker_names, Cluster = row_clusters)

# Check cluster sizes if needed
table(row_cluster_df$Cluster)

# --- DEFINE COLORS ---
cluster_colors <- c(
  "1" = "#E69F00", # orange
  "2" = "#56B4E9", # light blue
  "3" = "#009E73", # yellow
  "4" = "#F0E442", # bluish green
  "5" = "#0072B2", # blue
  "6" = "#D55E00", # vermillion
  "7" = "#CC79A7"  # reddish purple
)

# --- CUSTOM CLUSTER ORDER ---
custom_order <- c(1, 2, 3, 4, 5, 6, 7)  # <--- change this to whatever order you want

# Reorder rows and columns based on cluster membership
row_order <- order(factor(row_clusters, levels = custom_order))
col_order <- order(factor(col_clusters, levels = custom_order))

# Apply new order to correlation matrix
corr_matrix_reordered <- corr_matrix[row_order, col_order]

# Reorder cluster assignments accordingly
row_clusters_reordered <- row_clusters[row_order]
col_clusters_reordered <- col_clusters[col_order]

# --- ANNOTATIONS ---
library(ComplexHeatmap)
library(circlize)

row_ann <- rowAnnotation(
  Cluster = row_clusters_reordered,
  col = list(Cluster = cluster_colors),
  show_annotation_name = FALSE,
  width = unit(0.5, "cm")
)

col_ann <- HeatmapAnnotation(
  Cluster = col_clusters_reordered,
  col = list(Cluster = cluster_colors),
  show_annotation_name = FALSE,
  height = unit(0.5, "cm"),
  show_legend = FALSE
)

# --- COLOR SCALE ---
col_fun <- colorRamp2(
  c(-0.7, 0, 0.7),
  c("blue", "white", "red")
)

# --- PLOT HEATMAP ---
Heatmap(
  corr_matrix_reordered,
  col = col_fun,
  name = "Rho",  # Spearman correlation coefficient
  cluster_rows = FALSE,   # disable clustering since we ordered manually
  cluster_columns = FALSE,
  show_row_names = FALSE,
  show_column_names = FALSE,
  left_annotation = row_ann,
  top_annotation = col_ann,
  column_title = "Correlation matrix of log2FC response\n to transfusion (n = 249) with 146 non-bypass patients",
  column_title_gp = gpar(fontsize = 13)
)

# cluster 3 is most prominent - let's look into the 28 markers that are included!
heatmap_cluster_3 <- row_cluster_df$Gene[row_cluster_df$Cluster == 3]

#### UMAP CLUSTER ####
set.seed(123)
umap_input_mat <- t(as.matrix(filtered_markers))
LFC_umap_res <- umap(umap_input_mat, n_neighbors = 15, min_dist = 0.1, metric = "correlation")

# UMAP coordinates
LFC_umap_df <- as.data.frame(LFC_umap_res) 
LFC_umap_df$Marker <- rownames(umap_input_mat)

# UMAP without colored clusters
ggplot(LFC_umap_df, aes(x = V1, y = V2, label = Marker)) +
  geom_point() +
  geom_text(size = 2, vjust = -0.5)

# use elbow method to determine k
wss <- vector()
kmax <- 10  # test up to 10 clusters
for (k in 1:kmax) {
  km <- kmeans(umap_input_mat, centers = k, nstart = 25)
  wss[k] <- km$tot.withinss
}

# plot elbow curve - no prominent elbow
plot(1:kmax, wss, type = "b", pch = 19, frame = FALSE,
     xlab = "Number of clusters (k)",
     ylab = "Total within-cluster sum of squares")


# Plot UMAP using k-means clustering
k <- 6
set.seed(123)
clusters <- kmeans(LFC_umap_df[, 1:2], centers = k)$cluster
LFC_umap_df$Cluster <- as.factor(clusters)

ggplot(LFC_umap_df, aes(x = V1, y = V2, color = Cluster)) +
  geom_point(size = 2) +
  scale_color_manual(values = cluster_colors, drop = FALSE) +  # <- keeps unused colors stable
  theme_minimal() + ggtitle("UMAP of LFC of 249 markers across 146 non-bypass patients (k=6)")

# looking markers in cluster 3 (green)
markers_per_umap_cluster <- split(LFC_umap_df$Marker, LFC_umap_df$Cluster)
umap_cluster_3 <- markers_per_umap_cluster[[3]]
length(markers_per_umap_cluster[[3]]) #27 markers

#### COMPARING BETWEEN GEPHI, UMAP, AND HEATMAP ####

# combine all markers into one unique list
all_cluster_markers <- unique(c(umap_cluster_3, heatmap_cluster_3, gephi_cluster0_genes))

marker_presence_df <- data.frame(
  Marker = all_cluster_markers,
  in_umap_cluster_3 = all_cluster_markers %in% umap_cluster_3,
  in_heatmap_cluster_3 = all_cluster_markers %in% heatmap_cluster_3,
  in_gephi_cluster0 = all_cluster_markers %in% gephi_cluster0_genes
)

# find genes present in all three lists
common_cluster_markers <- marker_presence_df$Marker[
  marker_presence_df$in_umap_cluster_3 &
    marker_presence_df$in_heatmap_cluster_3 &
    marker_presence_df$in_gephi_cluster0
]
# 25 NULISA markers that show up in all three clustering methods
common_cluster_markers

# other markers that are missing in at least one
missing_markers <- marker_presence_df[rowSums(marker_presence_df[, -which(names(marker_presence_df) == "Marker")]) < ncol(marker_presence_df) - 1, ]

#### STRINGdb prep for 25 markers #####
uniprot_ids <- read.csv("/Users/catherinewang/Desktop/Sapru/TRAPMODS/STRINGdb\ stuff/TRAPMODS\ UNIPROT\ -\ Sheet1.csv", header = TRUE)

# get UNIPROT IDs for the 25 markers
common_cluster_markers_df <- data.frame(id = common_cluster_markers)
matched_25_markers <- left_join(common_cluster_markers_df, uniprot_ids, by = "id")
paste(matched_25_markers$uniprot, collapse = ",")

# get UNIPROT IDs for Gephi markers
gephi_cluster_markers_df <- data.frame(id = gephi_cluster0_genes)
matched_gephi_markers <- left_join(gephi_cluster_markers_df, uniprot_ids, by = "id")
paste(matched_gephi_markers$uniprot, collapse = ",")

# get UNIPROT IDs for heatmap markers
heatmap_cluster_markers_df <- data.frame(id = heatmap_cluster_3)
matched_heatmap_markers <- left_join(heatmap_cluster_markers_df, uniprot_ids, by = "id")
paste(matched_heatmap_markers$uniprot, collapse = ",")

# get UNIPROT IDs for UMAP markers
umap_cluster_markers_df <- data.frame(id = umap_cluster_3)
matched_umap_markers <- left_join(umap_cluster_markers_df, uniprot_ids, by = "id")
paste(matched_umap_markers$uniprot, collapse = ",")


#### MAPPING TO PELOD CLASSES ####
# Class 1 means patient outcomes are improving and PELOD score decreases
# Class 2 means patient outcomes are worsening and PELOD score increases

cluster_markers_expr <- filtered_markers[, colnames(filtered_markers) %in% gephi_cluster0_genes]
pelod_classes <- read.csv("clinical data/final/pelod_gmm_class.csv")
class1_patients <- pelod_classes$ParticipantID[pelod_classes$class == 1] #214/271 patients

# add PELOD classification to markers dataset
cluster_markers_expr$PELOD_class <- ifelse(
  rownames(cluster_markers_expr) %in% class1_patients, 1, 2
)
table(cluster_markers_expr$PELOD_class)
# 101 patients are Class 1
# 45 patients are Class 2

#### HEATMAP OF HOW CLUSTER EXPRESSION DIFFERS BETWEEN PELOD CLASSES ####

range(cluster_markers_expr, na.rm = TRUE)
col_fun_pelod <- colorRamp2(
  c(-10, 0, 10),
  c("blue", "white", "red")
)

pelod_colors <- c("1" = "#66C2A5", "2" = "#FC8D62")

ht_original <- Heatmap(
  as.matrix(cluster_markers_expr[, -ncol(cluster_markers_expr)]),
  row_split = cluster_markers_expr$PELOD_class,
  cluster_rows = TRUE
)

# Extract the clustered row order
original_order <- row_order(ht_original)

# Step 2: Apply the same row order but without row_split, and add color bar
row_ha <- rowAnnotation(
  PELOD_class = cluster_markers_expr$PELOD_class,
  col = list(PELOD_class = pelod_colors),
  show_annotation_name = FALSE,
  annotation_name_side = "top"
)

Heatmap(
  as.matrix(cluster_markers_expr[, -ncol(cluster_markers_expr)]),
  name = "Log2 Fold Change",
  col = col_fun_pelod,
  cluster_rows = FALSE,       
  cluster_columns = TRUE,     
  show_row_names = FALSE,
  show_column_names = TRUE,
  column_names_rot = 45,
  left_annotation = row_ha,
  row_order = unlist(original_order),  # preserve your old order
  column_title = "Cluster marker expression (n = 33) split by PELOD class across 146 patients"
)

#### MAPPING TO LCA ####
df_lca <- read.csv("/Users/catherinewang/Desktop/Sapru/TRAPMODS/clinical\ data/trapmods_lca.csv")
cluster_markers_expr$lca <- df_lca$lca[match(rownames(cluster_markers_expr), df_lca$sample)]

#### HEATMAP OF HOW CLUSTER EXPRESSION DIFFERS BETWEEN PELOD+LATENT CLASSES ####
lca_colors <- c("1" = "#E78AC3", "2" = "skyblue")

# Step 1: Get the row order from the original clustering
ht_original <- Heatmap(
  as.matrix(cluster_markers_expr[, -ncol(cluster_markers_expr)]),
  row_split = cluster_markers_expr$PELOD_class,
  cluster_rows = TRUE
)
original_order <- row_order(ht_original)

# Step 2: Add multiple row annotations (PELOD + LCA)
row_ha <- rowAnnotation(
  PELOD_class = cluster_markers_expr$PELOD_class,
  Latent_class = cluster_markers_expr$lca,
  col = list(
    PELOD_class = pelod_colors,
    Latent_class = lca_colors
  ),
  show_annotation_name = FALSE   # hides labels on the side, keeps color bars
)

# Step 3: Plot the main heatmap
Heatmap(
  as.matrix(cluster_markers_expr[, -( (ncol(cluster_markers_expr)-1) : ncol(cluster_markers_expr) )]),
  name = "Log2 Fold Change",
  col = col_fun_pelod,
  cluster_rows = FALSE,       
  cluster_columns = TRUE,     
  show_row_names = FALSE,
  show_column_names = TRUE,
  column_names_rot = 45,
  left_annotation = row_ha,
  row_order = unlist(original_order),
  column_title = "Cluster marker expression (n = 33) split by PELOD and LCA classes across 146 patients"
)

#### BAR CHART ####
cluster_markers_expr$mean_expr <- rowMeans(
  cluster_markers_expr[, !(colnames(cluster_markers_expr) %in% c("PELOD_class", "lca"))],
  na.rm = TRUE
)

cluster_markers_expr$cluster_expr_level <- ifelse(
  cluster_markers_expr$mean_expr > median(cluster_markers_expr$mean_expr, na.rm = TRUE),
  "High", "Low"
)

pelod_colors <- c("1" = "#66C2A5", "2" = "#FC8D62")
cluster_markers_expr$PELOD_class <- factor(cluster_markers_expr$PELOD_class, levels = c(1, 2))
# Plot
ggplot(cluster_markers_expr, aes(x = cluster_expr_level, fill = PELOD_class)) +
  geom_bar(position = "stack", color = "white") +
  scale_fill_manual(values = pelod_colors, name = "PELOD Class") +
  labs(
    x = "Gephi Cluster Expression Level",
    y = "Number of Patients",
    title = "Distribution of High/Low Cluster Expression by PELOD Class"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    axis.title = element_text(face = "bold"),
    legend.position = "right"
  )
