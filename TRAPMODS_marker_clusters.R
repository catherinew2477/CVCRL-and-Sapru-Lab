setwd("/Users/catherinewang/Desktop/Sapru/TRAPMODS/")

#### LIBRARIES ####
library(ComplexHeatmap)
library(pheatmap)
library(umap)
library(uwot) # newer version of umap library
library(dplyr)
library(ggplot2)


#### LOADING AND CLEANING ####

# markers <- read.table("network_markers_cw.csv", sep = ",", header = TRUE, row.names = 1)
# 
# # Noting patients with NAs and removing patients as needed
# rownames(markers)[apply(is.na(markers), 1, any)]
# # "26-015" "26-017" "39-002" "39-003" "44-001" "48-001" "59-007"
# NA_patients <- c("26-015","26-017","39-002","39-003","44-001","48-001","59-007")
# for (x in NA_patients) {
#   print(length(is.na(markers[x, ])))
# }
# # for patients with NAs, it's for all 249 markers --> remove these patients
# filtered_markers <- markers %>% filter(if_all(everything(), ~ !is.na(.))) # 146 patients

#### GENERATING HEATMAP OF LOG FOLD CHANGES ####

# generate heatmap of raw data after replacing NAs with 0s
markers_changedto_0 <- markers
markers_changedto_0[is.na(markers_changedto_0)] <- 0
Heatmap(
  as.matrix(markers_changedto_0),
  name = "logFC",
  cluster_columns = TRUE,
  cluster_rows = TRUE,
  row_names_gp = gpar(fontsize = 4),
  column_names_gp = gpar(fontsize = 3),
  column_names_rot = 45  # rotate column labels diagonally
)

# heatmap of raw data after filtering out NA patients
Heatmap(
  as.matrix(filtered_markers),
  name = "logFC",
  cluster_columns = TRUE,
  cluster_rows = TRUE,
  row_names_gp = gpar(fontsize = 4),
  column_names_gp = gpar(fontsize = 3),
  column_names_rot = 45  # rotate column labels diagonally
)


# way too many markers in this heatmap --> determine number of clusters to divide up
hc <- hclust(dist(t(filtered_markers))) 
# dist(t(markers)) calculates distance between markers based on their expression across patients
plot(hc)

k_num <- 5
clusters <- cutree(hc, k = k_num)  

# compute average expression of markers across patients in each cluster 
cluster_avg <- sapply(1:k_num, function(k) {
  rowMeans(filtered_markers[, clusters == k, drop = FALSE])
})

# generate heatmap with clusters
pheatmap(cluster_avg, fontsize = 5)
# first cluster has most prominent contrast...
# Each column is the average expression of markers from cluster ___ 
# Each row is a patient

# seeing what markers are in each cluster
heatmap_marker_clusters <- split(names(clusters), clusters)
heatmap_marker_clusters[[1]]
heatmap_marker_clusters[[2]]
heatmap_marker_clusters[[3]]
heatmap_marker_clusters[[4]]
heatmap_marker_clusters[[5]]

# markers in each cluster
table(clusters)
markers_per_cluster <- split(names(clusters), clusters)
markers_per_cluster[[1]]

#### UMAP OF MARKERS ####
# convert to numeric matrix
marker_mat <- t(as.matrix(filtered_markers)) # transpose to make the rows the markers
storage.mode(marker_mat) <- "numeric"

set.seed(123)
# euclidean or correlation metric?
raw_umap_result <- umap(marker_mat, n_neighbors = 15, min_dist = 0.1, metric = "correlation")

# UMAP coordinates
raw_umap_df <- as.data.frame(raw_umap_result) # add $layout if using umap package instead of uwot
raw_umap_df$Marker <- rownames(marker_mat)

# UMAP without colored clusters
ggplot(raw_umap_df, aes(x = V1, y = V2, label = Marker)) +
  geom_point() +
  geom_text(size = 2, vjust = -0.5)

# use elbow method to determine k?

mat <- t(filtered_markers)
wss <- vector()
kmax <- 10  # test up to 10 clusters
for (k in 1:kmax) {
  km <- kmeans(mat, centers = k, nstart = 25)
  wss[k] <- km$tot.withinss
}

# plot elbow curve
plot(1:kmax, wss, type = "b", pch = 19, frame = FALSE,
     xlab = "Number of clusters (k)",
     ylab = "Total within-cluster sum of squares")

# 4 or 5 clusters looks optimal?

# Example with kmeans (replace k with whatever you want)
k <- 6
set.seed(123)
clusters <- kmeans(raw_umap_df[, 1:2], centers = k)$cluster
raw_umap_df$Cluster <- as.factor(clusters)
#### Fixed color scheme ####
fixed_colors <- c(
  "1" = "blue",
  "2" = "#ff785f",
  "3" = "green",
  "4" = "purple",
  "5" = "orange",
  "6" = "cyan"
  # extend if needed
)

ggplot(raw_umap_df, aes(x = V1, y = V2, color = Cluster)) +
  geom_point(size = 2) +
  scale_color_manual(values = fixed_colors, drop = FALSE) +  # <- keeps unused colors stable
  theme_minimal() + ggtitle("UMAP of raw LFC with k-means clustering (k=6)")


# from correlation heatmap
markers_in_row_cluster3 <- rownames(corr_matrix)[row_clusters == 3]
length(markers_in_row_cluster3)

markers_per_cluster <- split(raw_umap_df$Marker, raw_umap_df$Cluster)
markers_per_cluster[[1]]
length(markers_per_cluster[[1]])
markers_per_cluster[[2]] 
length(markers_per_cluster[[2]]) #27 markers in the red UMAP cluster
markers_per_cluster[[3]]
length(markers_per_cluster[[3]])
markers_per_cluster[[4]]
length(markers_per_cluster[[4]])
markers_per_cluster[[5]]
length(markers_per_cluster[[5]])
markers_per_cluster[[6]]
length(markers_per_cluster[[6]])

#### Comparing UMAP markers clusters to red Gephi cluster ####

sum(markers_per_cluster[[3]] %in% nodes_class0$Id) #27/28
setdiff(markers_per_cluster[[2]], nodes_class0$Id)
# THPO is present in UMAP cluster and not present in Gephi cluster
setdiff(nodes_class0$Id, markers_per_cluster[[2]])
# "TNFSF13" "CCL28"   "HGF"     "TNFSF12" "IL1B"    "FLT1"
# markers that are present in Gephi but not in the UMAP cluster

#### Comparing UMAP cluster 2 to cluster 3 of the heatmap ####
sum(markers_per_cluster[[2]] %in% rownames(corr_matrix)[row_clusters == 3]) #25
setdiff(markers_per_cluster[[2]], rownames(corr_matrix)[row_clusters == 3]) 
# FTH1 and NAMPT and THPO
setdiff(rownames(corr_matrix)[row_clusters == 3], markers_per_cluster[[2]]) 
# "CCL8"    "IL1B"    "TNFSF13"

# list of markers from heatmap cluster 3
corr_heatmap_markers <- rownames(corr_matrix)[row_clusters == 3]
cat(corr_heatmap_markers, sep = ",  ")

# list of markers from red Gephi cluster)
cat(gephi_cluster0_genes, sep = ", ")
length(gephi_cluster0_genes)

# define sets
umap2 <- markers_per_cluster[[2]]
gephi_red <- nodes_class0$Id
heatmap3 <- rownames(corr_matrix)[row_clusters == 3]

# combine all unique markers from the three clusters
all_markers <- unique(c(umap2, gephi_red, heatmap3)) 
length(all_markers) # 35 markers

# make table
summary_tbl <- data.frame(
  Marker = all_markers,
  In_UMAP2 = all_markers %in% umap2,
  In_Gephi = all_markers %in% gephi_red,
  In_Heatmap3 = all_markers %in% heatmap3
)
write.csv(summary_tbl, "comparing_markers_tbl.csv")

#### putting cluster 2 from UMAP into STRINGdb ####
stringdb_raw_umap_df <- raw_umap_df[, 3:4]
stringdb_raw_umap_df <-stringdb_raw_umap_df %>% tibble::rownames_to_column(var = "id")
stringdb_raw_umap_df <- stringdb_raw_umap_df[, c(1,3)]
uniprot_ids <- read.csv("STRINGdb stuff/TRAPMODS UNIPROT - Sheet1.csv", header = TRUE)
stringdb_raw_umap_df <- left_join(stringdb_raw_umap_df, uniprot_ids, by = "id")

# extract uniprot ids from each cluster and paste into STRINGdb
paste((stringdb_raw_umap_df %>% filter(Cluster == "1"))$uniprot, collapse = ",")
paste((stringdb_raw_umap_df %>% filter(Cluster == "2"))$uniprot, collapse = ",")
paste((stringdb_raw_umap_df %>% filter(Cluster == "3"))$uniprot, collapse = ",")
paste((stringdb_raw_umap_df %>% filter(Cluster == "4"))$uniprot, collapse = ",")
paste((stringdb_raw_umap_df %>% filter(Cluster == "5"))$uniprot, collapse = ",")


#### CLUSTERING BY PATIENTS (UMAP) USING 249 MARKERS ####
raw_umap_result_bypatient <- umap(t(marker_mat), metric = "correlation") # transpose back to make the rows be the patients

# UMAP coordinates
raw_umap_df_bypatient <- as.data.frame(raw_umap_result_bypatient) 
raw_umap_df_bypatient$Patient <- rownames(t(marker_mat))

# UMAP without colored clusters
ggplot(raw_umap_df_bypatient, aes(x = V1, y = V2, label = Patient)) +
  geom_point() +
  geom_text(size = 2, vjust = -0.5)

k <- 5
set.seed(123)
clusters <- kmeans(raw_umap_df_bypatient[, 1:2], centers = k)$cluster
raw_umap_df_bypatient$Cluster <- as.factor(clusters)
#### Fixed color scheme ####
fixed_colors <- c(
  "1" = "blue",
  "2" = "#ff785f",
  "3" = "green",
  "4" = "purple",
  "5" = "orange"
  # extend if needed
)

ggplot(raw_umap_df_bypatient, aes(x = V1, y = V2, color = Cluster)) +
  geom_point(size = 2) +
  scale_color_manual(values = fixed_colors, drop = FALSE) +  # <- keeps unused colors stable
  theme_minimal() + ggtitle("UMAP of patients with k-means clustering (k=5)")

patients_per_cluster <- split(raw_umap_df_bypatient$Patient, raw_umap_df$Cluster)

# Mapping clustered patients to latent classes
latent_classes <- read.csv("clinical data/trapmods_lca.csv")

# turn patients_per_cluster into a data frame
patients_per_cluster_df <- stack(patients_per_cluster)
colnames(patients_per_cluster_df) <- c("Patient", "Cluster")

# join with latent_classes by Patient
matched_patients_df <- patients_per_cluster_df %>%
  left_join(latent_classes, by = c("Patient" = "sample"))

umap_with_latent_classes <- raw_umap_df_bypatient %>%
  left_join(latent_classes, by = c("Patient" = "sample"))

# convert to factor so the color matching works
umap_with_latent_classes$lca <- factor(umap_with_latent_classes$lca)

latent_class_colors <- c(
  "1" = "blue",
  "2" = "#ff785f"
)

ggplot(umap_with_latent_classes, aes(x = V1, y = V2, color = lca)) +
  geom_point(size = 2, alpha = 0.8) +
  scale_color_manual(values = latent_class_colors) +
  theme_minimal() +
  labs(color = "Latent Class",
       title = "UMAP of Patients Colored By Latent Class")


#### CLUSTERING PATIENTS BASED ON 33 GEPHI MARKERS ####
markers_gephi_subset <- filtered_markers[, colnames(markers) %in% gephi_cluster0_genes]
patient_umap_gephi_subset <- umap(markers_gephi_subset, metric = "correlation")
patient_umap_gephi_subset <- as.data.frame(patient_umap_gephi_subset) 
patient_umap_gephi_subset$Patient <- rownames(markers_gephi_subset)

# UMAP without colored clusters
ggplot(patient_umap_gephi_subset, aes(x = V1, y = V2, label = Patient)) +
  geom_point() +
  geom_text(size = 2, vjust = -0.5)

# make elbow plot to determine k
set.seed(123)

wss <- sapply(1:10, function(k){
  kmeans(patient_umap_gephi_subset[, 1:2], centers = k, nstart = 25)$tot.withinss
})

plot(1:10, wss, type = "b", pch = 19, frame = FALSE,
     xlab = "Number of clusters k",
     ylab = "Total within-clusters sum of squares",
     main = "Elbow Plot for K-means clustering of patients with 33 markers")

# 3 clusters looks optimal

k <- 3
set.seed(123)
clusters <- kmeans(patient_umap_gephi_subset[, 1:2], centers = k)$cluster
patient_umap_gephi_subset$Cluster <- as.factor(clusters)
#### Fixed color scheme ####
fixed_colors <- c(
  "1" = "blue",
  "2" = "#ff785f",
  "3" = "green"
  # extend if needed
)

ggplot(patient_umap_gephi_subset, aes(x = V1, y = V2, color = Cluster)) +
  geom_point(size = 2) +
  scale_color_manual(values = fixed_colors, drop = FALSE) +  # <- keeps unused colors stable
  theme_minimal() + ggtitle("UMAP of patients based on expression profiles of Gephi cluster (k=3)")

# color UMAP by latent class

umap_with_latent_classes <- patient_umap_gephi_subset %>%
  left_join(latent_classes, by = c("Patient" = "sample"))

umap_with_latent_classes$lca <- factor(umap_with_latent_classes$lca)

latent_class_colors <- c(
  "1" = "blue",
  "2" = "#ff785f"
)

ggplot(umap_with_latent_classes, aes(x = V1, y = V2, color = lca)) +
  geom_point(size = 2, alpha = 0.8) +
  scale_color_manual(values = latent_class_colors) +
  theme_minimal() +
  labs(color = "Latent Class",
       title = "UMAP of Patients (Using Gephi Markers) Colored By Latent Class")


#### UPDATED LFC DATASET: PELOD trajectories and mortality outcomes ####
