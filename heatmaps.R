library(ComplexHeatmap)
library(circlize)

library(readxl)
sicol8_sictrl_counts <- read_excel("RNAseq siCOL8A1_for heatmap 2.xlsx")
sicol8_sictrl_counts_matrix <- as.matrix(sicol8_sictrl_counts[, 2:7])
sicol8_sictrl_counts_matrix <- apply(sicol8_sictrl_counts_matrix, 2, as.numeric)

# set row names from the first column
rownames(sicol8_sictrl_counts_matrix) <- sicol8_sictrl_counts[[1]]
colnames(sicol8_sictrl_counts_matrix) <- c("Ctrl siRNA 1", "Ctrl siRNA 2", "Ctrl siRNA 3", "COL8A1 siRNA 1", "COL8A1 siRNA 2", "COL8A1 siRNA 3")

pheatmap(sicol8_sictrl_counts_matrix,
         scale = "row",                  # Scale data by row (z-score)
         cluster_rows = FALSE,           # Disable row clustering if needed
         cluster_cols = FALSE,           # Disable column clustering if needed
         color = colorRampPalette(c("blue", "black", "red"))(50),
         show_rownames = TRUE,           # Show row names
         show_colnames = TRUE,           # Show column names
         name = "Z-score",               # Title for the color scale
         column_names_side = "top", 
         cellwidth = 20, cellheight = 10)

#gene annotations
color_palette <- c(
  "Cell adhesion molecule binding\nCytoskeletal protein binding\nCadherin binding" = "#1f77b4", # blue
  "Tubulin binding" = "#ff7f0e",                # orange
  "Cell cycle" = "#2ca02c",                     # green
  "Kinase activity" = "#27d6cb",         # teal
  "Steroid biosynthesis" = "#b373ee"            # purple
)

gene_annotation <- data.frame(
  Function = c(rep("Cell adhesion molecule binding\nCytoskeletal protein binding\nCadherin binding", 12),
               rep("Tubulin binding", 7),
               rep("Cell cycle", 7),
               rep("Kinase activity", 10),
               rep("Steroid biosynthesis", 4)),
  row.names = c("ANLN", "CTNNA1", "MARK2", "COLBLL1", "MICALL1", "NOP56", "STK38", "KDR", "TLN1", "USP8", "ARHGAP18", "PTPN11",
                "TPX2", "KIF4A", "CCDC88A", "CENPF", "DLGAP5", "APC", "RACGAP1", 
                "CDC27", "YWHAH", "CCNA1", "CCNA2", "BUB1B", "CDK1", "TFDP1",
                "BAZ1B", "PRKCH", "PIK3CB", "HIPK3", "EPHB2", "TBCK", "STK17A", "PKIA", "TOP1", "SGMS1",
                "MSMO1", "LIPA", "SC5D", "HMGCR"))

# deleted duplicates: MARK2 from Kinase activity, STK38 from Kinase activity, and BUB1B from Kinase activity

gene_annotation$Function <- factor(gene_annotation$Function,
                                   levels = c("Cell adhesion molecule binding\nCytoskeletal protein binding\nCadherin binding",
                                              "Tubulin binding",
                                              "Cell cycle",
                                              "Kinase activity", 
                                              "Steroid biosynthesis"
                                              ))
annotation_colors <- list(
  Function = color_palette[unique(gene_annotation$Function)]
)

pheatmap(sicol8_sictrl_counts_matrix,
         scale = "row",
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         color = colorRampPalette(c("blue", "black", "red"))(50),
         show_rownames = TRUE,
         show_colnames = TRUE,
         name = "Z-score",
         column_names_side = "top",
         cellwidth = 20, cellheight = 10,
         annotation_row = gene_annotation,
         annotation_colors = annotation_colors,
         annotation_names_row = FALSE)
