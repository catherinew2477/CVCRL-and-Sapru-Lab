# read in both datasets (significant genes from each)
library(dplyr)
set1 <- read.csv("sig_genes_siRNA.csv")
colnames(set1)[1] = "Gene"
set2 <- read.csv("sig_genes_recombinant.csv")
colnames(set2)[1] = "Gene"

#find genes that are in both datasets for fun idk
common_genes <- intersect(set1$Gene, set2$Gene)

#find genes that are oppositely regulated
opposite_regulated <- inner_join(set1, set2, by = "Gene", suffix = c(".set1", ".set2")) %>%
  filter((log2FoldChange.set1 > 0 & log2FoldChange.set2 < 0) | (log2FoldChange.set1 < 0 & log2FoldChange.set2 > 0))

# rank the 27 genes from biggest log2foldchanges
i <- 1
diff <- c()
while (i<28) {
  Log2FoldChange_diff <- abs((opposite_regulated$log2FoldChange.set1)[i]) + abs(opposite_regulated$log2FoldChange.set2[i])
  diff = c(diff, Log2FoldChange_diff)
  i <- i+1
}
opposite_regulated$Difference_log2FoldChange <- diff
opposite_regulated_sorted <- opposite_regulated %>% arrange(desc(Difference_log2FoldChange))
write.csv(opposite_regulated_sorted, "opposite_regulated_genes_sorted.csv")

# separate genes --> (1) positive log2foldchange in set 1 and negative in set2, (2) negative log2foldchange in set1 and positive in set 2
up1_down2 <- read.csv("opposite_regulated_genes_sorted.csv", header = TRUE, row.names = 1)
down1_up2 <- read.csv("opposite_regulated_genes_sorted.csv", header = TRUE, row.names = 1)

pos <- which(up1_down2$log2FoldChange.set1 > 0)
rows <- rownames(up1_down2)[pos]
up1_down2 <- up1_down2[rows,]
write.csv(up1_down2, "up1_down2.csv")

pos <- which(down1_up2$log2FoldChange.set1 < 0)
rows <- rownames(down1_up2)[pos]
down1_up2 <- down1_up2[rows,]
write.csv(down1_up2, "down1_up2.csv")


