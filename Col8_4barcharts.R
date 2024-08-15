library(ggplot2)

##setting up data frames
A1_data <- read.table("COL8A1_pTPM.tsv", sep = '\t', header = TRUE)
A1_transcripts <- A1_data[c("Tissue","nTPM")]
A1_proteomic <- A1_data[c("Tissue","pTPM")]

A2_data <- read.table("COL8A2_pTPM.tsv", sep = '\t', header = TRUE)
A2_transcripts <- A2_data[c("Tissue", "nTPM")]
A2_proteomic <- A2_data[c("Tissue","pTPM")]

tissues <- A1_data$Tissue

#for scaling purposes
max(A1_transcripts$nTPM) #36
max(A2_transcripts$nTPM) #19.3
max(A1_proteomic$pTPM) #60.9
max(A2_proteomic$pTPM) #21.8

## generate bar chart for COL8A1 transcript distribution
ggplot(A1_transcripts, aes(x = factor(tissues), y = A1_transcripts$nTPM)) + geom_bar(stat = "identity") + 
  xlab("Tissues") + ylab("Normalized transcripts per million") + ggtitle("COL8A1 gene expression levels across tissues") + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + ylim(0,38)

## generate bar chart for COL8A2 transcript distribution
ggplot(A2_transcripts, aes(x = factor(tissues), y = A2_transcripts$nTPM)) + geom_bar(stat = "identity") + 
  xlab("Tissues") + ylab("Normalized transcripts per million") + ggtitle("COL8A2 gene expression levels across tissues") + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + ylim(0,38)

## generate bar chart for COL8A1 transcript distribution
ggplot(A1_proteomic, aes(x = factor(tissues), y = A1_proteomic$pTPM)) + geom_bar(stat = "identity") + 
  xlab("Tissues") + ylab("Protein transcripts per million") + ggtitle("COL8A1 protein abundance across tissues") + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + ylim(0,65)

## generate bar chart for COL8A2 transcript distribution
ggplot(A2_proteomic, aes(x = factor(tissues), y = A2_proteomic$pTPM)) + geom_bar(stat = "identity") + 
  xlab("Tissues") + ylab("Protein transcripts per million") + ggtitle("COL8A2 protein abundance across tissues") + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + ylim(0,65)

## trying to plot two bars side by side for each tissue...

# comparing gene expression levels
nTPMdata <- read.table("bothgenes.tsv", sep = '\t', header = TRUE)
nTPMdata <- data[c("Gene.name","Tissue","nTPM")]

ggplot(nTPMdata, aes(x = Tissue, y = nTPM, fill = Gene.name)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  labs(title = "Comparing gene expression levels of COL8A1 and COL8A2 across tissues",
       x = "Tissue", y = "Normalized transcripts per million",
       fill = "Gene") +
  theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

# comparing protein transcripts
pTPMdata <- read.table("bothgenes.tsv", sep = '\t', header = TRUE)
pTPMdata <- pTPMdata[c("Gene.name","Tissue","pTPM")]

ggplot(pTPMdata, aes(x = Tissue, y = pTPM, fill = Gene.name)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  labs(title = "Comparing protein abundance of COL8A1 and COL8A2 across tissues",
       x = "Tissue", y = "Protein transcripts per million",
       fill = "Gene") +
  theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))


