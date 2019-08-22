

library("DESeq2")
#setwd("/home/thimmamp/Downloads/BESE_DE/AnalysingUUO")
filename <- "BetweenUUO_"
raw_data <- read.csv("fibrosis_counts.csv", header=TRUE, row.names = 1)

str(raw_data)
### raw RNA-Seq count distribution
### From the above file we are going to use only mutated samples ie with UUO

##inhouse processed
selected_data <- raw_data[, c(4:7, 11:14)] ### to get wt_uu0 and smoc_uuo
str(selected_data)
selected_data <- as.matrix(selected_data)

##Create the DESeq2 data object
genotype <- c(rep("WT", 4),rep("TRTED", 4))  ### for UUO
coldata <- data.frame(genotype)
rownames(coldata) <- colnames(selected_data)
coldata
### Create the DESeq dataset
dds <- DESeqDataSetFromMatrix(countData = selected_data, colData = coldata, design = ~ genotype)
summary(dds)

## Pre-filtering ; Filter genes that have 0 count in all samples ----
keep <- rowSums(counts(dds)) > 0 # 0 means no pre-filtering
dds <- dds[keep,]
summary(dds)
dds$genotype <- factor(dds$genotype, levels = c("WT", "TRTED"))
dds$genotype
dds <- estimateSizeFactors(dds)
normalized_data <- counts(dds, normalized=TRUE)
#Run DESeq2
dds <- DESeq(dds)
structure(dds)
res <- results(dds)
res_nona <- na.omit(res)
dim(res_nona)
write.table(res_nona, file="DE_genes_all_UUO.txt", sep="\t")
