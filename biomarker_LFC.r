#### LIBRARIES ####
library(rstatix)
library(tidyr)
library(ComplexHeatmap)
library(circlize)
library(pheatmap)
library(umap)
library(uwot) # newer version of umap library
library(dplyr)
library(tibble)
library(ggplot2)


#### LOAD AND CLEAN HEME + INFLAMMATORY MARKER DATA ####

# heme + 249 inflammatory markers
heme_inf_markers <- read.csv("heme_inflamm.csv")
heme_inf_markers <- heme_inf_markers[heme_inf_markers$scoreday %in% c(0, 1), ]


# remove CTSS because of low expression
heme_inf_markers <-heme_inf_markers %>% select(-CTSS)

# log transform heme only
heme_inf_markers$heme <- log1p(heme_inf_markers$heme)

day0_heme_inf <- heme_inf_markers %>% filter(scoreday == 0)
day1_heme_inf <- heme_inf_markers %>% filter(scoreday == 1)

# join day0 and day 1 by ParticipantID
merged_heme_inf <- day1_heme_inf %>%
  select(-scoreday) %>%
  inner_join(day0_heme_inf %>% select(-scoreday), by = "ParticipantID", suffix = c("_day1", "_day0"))

rownames(merged_heme_inf) <- merged_heme_inf$ParticipantID
merged_heme_inf <- merged_heme_inf %>% select(-ParticipantID)

# get the marker names (without _day1/_day0)
heme_inf_marker_names <- unique(sub("_day[01]$", "", colnames(merged_heme_inf)))

# calculate log2 fold changes
heme_inf_logfc <- sapply(heme_inf_marker_names, function(m) {
  log2((merged_heme_inf[[paste0(m, "_day1")]]) /   
                (merged_heme_inf[[paste0(m, "_day0")]]))
})

heme_inf_logfc <- as.data.frame(heme_inf_logfc)


# convert to data frame and keep rownames as PatientID
heme_inf_logfc <- as.data.frame(heme_inf_logfc)
rownames(heme_inf_logfc) <- rownames(merged_heme_inf)

#### LOAD AND CLEAN COAGULATION MARKER DATA ####

coag_markers <- read.csv("coagulation.csv")
coag_markers <- coag_markers[coag_markers$scoreday %in% c(0, 1), ]

# log transform everything
coag_markers[, 3:11] <- log1p(coag_markers[, 3:11])

day0_coag <- coag_markers %>% filter(scoreday == 0)
day1_coag <- coag_markers %>% filter(scoreday == 1)

# join day0 and day 1 by ParticipantID
merged_coag <- day1_coag %>%
  select(-scoreday) %>%
  inner_join(day0_coag %>% select(-scoreday), by = "ParticipantID", suffix = c("_day1", "_day0"))

rownames(merged_coag) <- merged_coag$ParticipantID
merged_coag <- merged_coag %>% select(-ParticipantID)

# get the marker names (without _day1/_day0)
coag_marker_names <- unique(sub("_day[01]$", "", colnames(merged_coag)))

# calculate log2 fold changes
coag_logfc <- sapply(coag_marker_names, function(m) {
  log2((merged_coag[[paste0(m, "_day1")]]) / 
         (merged_coag[[paste0(m, "_day0")]]))
})


# convert to data frame and keep rownames as PatientID
coag_logfc <- as.data.frame(coag_logfc)
rownames(coag_logfc) <- rownames(merged_coag)


heme_inf_logfc <- heme_inf_logfc %>% rownames_to_column(var = "ParticipantID")
coag_logfc <- coag_logfc %>% rownames_to_column(var = "ParticipantID")

# FINAL: 226 patients with 249 inflammatory markers, heme, and 9 coagulation markers
merged_markers <- heme_inf_logfc %>% left_join(coag_logfc, by = "ParticipantID")
rownames(merged_markers) <- merged_markers$ParticipantID
merged_markers$ParticipantID <- NULL

# final csv of LFC
write.csv(merged_markers, "all_TRAPMODS_markers_LFC.csv")

# Make a note of how many patients have all NULISA profiles, all coagulation profiles, and all heme
sum(complete.cases(merged_markers[, 2:250])) # NULISA - 203
sum(complete.cases(merged_markers[, 1])) # heme - 219
sum(complete.cases(merged_markers[, 251:259])) #- 222

#### CORRELATION MATRIX HEATMAP ####


#### SUBSETTING PATIENTS BASED ON MARKER PROFILES AND BYPASS GROUP ####
# Take the subset of patients that have inflammatory, thrombosis, and heme (complete cases) markers

bypass_groups <- read.csv("bypass_groups.csv")
bypass_groups$ICUAdmissionType <- ifelse(bypass_groups$ICUAdmissionType == "Surgical, cardiac", 1, 0)
sum(bypass_groups$ICUAdmissionType) #107 bypass patients 
bypass_patients <- bypass_groups$ParticipantID[bypass_groups$ICUAdmissionType == 1]
non_bypass_patients <- bypass_groups$ParticipantID[bypass_groups$ICUAdmissionType == 0]

# add bypass classification to markers dataset with complete profiles
merged_markers <- merged_markers %>%
  mutate(bypass = ifelse(rownames(merged_markers) %in% bypass_patients, 1, 0))

non_bypass_patients_NULISA <- rownames(
  merged_markers[merged_markers$bypass == 0 & complete.cases(merged_markers[, 2:250]), ]
)
# 133 non-bypass patients with complete NULISA profiles
bypass_patients_NULISA <- rownames(
  merged_markers[merged_markers$bypass == 1 & complete.cases(merged_markers[, 2:250]), ]
)
# 70 bypass patients with complete NULISA profiles
# adds up to 203 patients with complete NULISA profiles yayyy

complete_markers_subset <- merged_markers[complete.cases(merged_markers), ] #195 patients
sum(complete_markers_subset$bypass == 0) # 125


# Generate correlation matrix
corr_matrix <- cor(complete_markers_subset, method = "spearman", use = "pairwise.complete.obs")

# get p-values
p_vals <- merged_markers %>% cor_test(method = "spearman", )
p_vals$p_adjusted <- p.adjust(p_vals$p, method = "fdr")

# Convert cor_test() results to matrix
adj_p_vals_matrix <- p_vals %>% 
  select(var1, var2, p_adjusted) %>%
  pivot_wider(names_from = var2, values_from = p_adjusted) %>%
  column_to_rownames("var1")

#### GENERATE COLOR-CODED HEATMAP WITH CLUSTERS ####

# --- CLUSTERING ---
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
  column_title = "Correlation matrix of fold change response to transfusion (n = 259)"
)

# seeing what the markers are per cluster
marker_names <- rownames(corr_matrix)

# Create a data frame mapping each gene to its cluster
row_cluster_df <- data.frame(
  Gene = marker_names,
  Cluster = row_clusters
)

row_cluster_list <- split(row_cluster_df$Gene, row_cluster_df$Cluster)
row_cluster_list[["1"]]
row_cluster_list[["2"]]
row_cluster_list[["3"]]
row_cluster_list[["4"]]
table(row_cluster_df$Cluster)


#### UMAP OF CLUSTERED MARKERS ####
#### UMAP OF MARKERS ####
# convert to numeric matrix
set.seed(123)

lfc_umap_res <- umap(t(filtered_markers), n_neighbors = 15, min_dist = 0.1, metric = "correlation")

lfc_umap_df <- as.data.frame(lfc_umap_res) 
lfc_umap_df$Marker <- colnames(filtered_markers)

# UMAP without colored clusters
ggplot(lfc_umap_df, aes(x = V1, y = V2, label = Marker)) +
  geom_point() +
  geom_text(size = 2, vjust = -0.5)


# elbow plot to find k
wss <- vector()
kmax <- 10  # test up to 10 clusters
for (k in 1:kmax) {
  km <- kmeans(t(filtered_markers), centers = k, nstart = 25)
  wss[k] <- km$tot.withinss
}

# plot elbow curve
plot(1:kmax, wss, type = "b", pch = 19, frame = FALSE,
     xlab = "Number of clusters (k)",
     ylab = "Total within-cluster sum of squares")

k <- 4
set.seed(123)
clusters <- kmeans(lfc_umap_df[, 1:2], centers = k)$cluster
lfc_umap_df$Cluster <- as.factor(clusters)
#### Fixed color scheme ####
fixed_colors <- c(
  "1" = "blue",
  "2" = "#ff785f",
  "3" = "#009E73",
  "4" = "purple",
  "5" = "orange",
  "6" = "cyan"
  # extend if needed
)

ggplot(lfc_umap_df, aes(x = V1, y = V2, color = Cluster)) +
  geom_point(size = 2) +
  scale_color_manual(values = fixed_colors, drop = FALSE) +  # <- keeps unused colors stable
  theme_minimal() + ggtitle("UMAP of raw LFC with k-means clustering (k=4)")

