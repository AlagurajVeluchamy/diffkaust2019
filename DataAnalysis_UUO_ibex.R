##### DE using DESeq2 for BESE workshop
# library(BiocManager)
#BiocManager::install("vsn")
#BiocManager::install('EnhancedVolcano')
#BiocManager::install('pheatmap')
#BiocManager::install('hexbin')
#BiocManager::install("apeglm")
#BiocManager::install("annotables")
#BiocManager::install('ggrepel')
#install.packages("devtools")
#devtools::install_github("stephenturner/annotables")
#BiocManager::install("gProfileR")

library("DESeq2")
library(ggplot2)
library("RColorBrewer")
library("pheatmap")
library(biobroom)
library("vsn")
library(apeglm)
library(annotables)
library(biobroom)
library(dplyr)
library(ggrepel)
library("EnhancedVolcano")
library(gProfileR)

#setwd("/home/thimmamp/Downloads/BESE_DE/AnalysingUUO")
filename <- "BetweenUUO_"
raw_data <- read.csv("fibrosis_counts.csv", header=TRUE, row.names = 1)

### Explore the structure of raw counts
str(raw_data)
### raw RNA-Seq count distribution
### From the above file we are going to use only mutated samples ie with UUO

##inhouse processed
selected_data <- raw_data[, c(4:7, 11:14)] ### to get wt_uu0 and smoc_uuo
str(selected_data)
selected_data <- as.matrix(selected_data)

###Before normalisation counts plot
dat_log2 <- stack(as.data.frame(log2(selected_data)))

pdf(file = paste0(filename, "raw_count_distribution_ofuuosamples.pdf"))
ggplot(data = dat_log2, mapping = aes(x = ind, y = values)) +
  geom_jitter(alpha = 0.3, color = "tomato") +
  geom_boxplot(alpha = 0)
dev.off()

##Experiment meta data ----
##Create the DESeq2 data object
genotype <- c(rep("WT", 4),rep("TRTED", 4))  ### for UUO
coldata <- data.frame(genotype)
rownames(coldata) <- colnames(selected_data)
coldata
### Create the DESeq dataset
dds <- DESeqDataSetFromMatrix(countData = selected_data,
                              colData = coldata,
                              design = ~ genotype)
dds
summary(dds)

###
## Pre-filtering ; Filter genes that have 0 count in all samples ----
keep <- rowSums(counts(dds)) > 0 # 0 means no pre-filtering
dds <- dds[keep,]
summary(dds)
dds$genotype <- factor(dds$genotype, levels = c("WT", "TRTED"))
dds$genotype

###now have a DESeq2 object storing our raw counts and metadata that we can use to explore the data with DESeq2 
# functions and to use for performing the differential expression analysis. 
###Normalization ----
### We have created the DESeq2 object and now wish to perform quality control on our samples. Therefore, we 
# need to generate the normalized counts (normalized for library size, which is the total number of gene counts 
# per sample, while accounting for library composition). To obtain the normalized counts, use the DESeq2 object 
# and generate the normalized counts matrix.

dds <- estimateSizeFactors(dds)
sizeFactors(dds)
View(counts(dds))  ### raw counts

##Normalize the data for library size as normalization factor
normalized_data <- counts(dds, normalized=TRUE)
normalized_data
View(normalized_data)


### visualize normalized data
dat_log2 <- stack(as.data.frame(log2(normalized_data)))
###Stacking vectors concatenates multiple vectors into a single vector along with a factor indicating where 
#each observation originated. Unstacking reverses this operation. (from utils package)
pdf(file = paste0(filename,"CountDistribution_after_normalisation.pdf"))
ggplot(data = dat_log2, mapping = aes(x = ind, y = values)) +
  geom_jitter(alpha = 0.3, color = "forestgreen") + 
  geom_boxplot(alpha = 0) 
  #geom_jitter(alpha = 0.3, color = "darkred") +  #color = gold1, tomato, forestgreeen
dev.off()

###we now have our normalized counts, which we can use to accurately compare gene expression between samples. 
#We will be using the normalized counts to explore similarities in gene expression between 
# each of our samples, with the expection that our biological replicates are more similar to each other and 
#the different conditions (wild type and fibrosis) are more different. 
  
#### Clustering of samples ----
#Result: Heatmap of all sample distance & PCA plot

####Unsupervised clustering analysis ####
#Quickly estimate dispersion trend and apply a variance stabilizing transformation
#blind option is to mention the VST transformation must be blind to sample information
#Hierarchical clustering with correlation heatmaps
#rld <- rlog(dds, blind=FALSE) # Extracting transformed values
# Dispersion in statistics is a way of describing how spread out a set of data is.
vsd <- vst(dds, blind=FALSE)
# Extract the vst matrix from the object
vsd_mat_wt <- assay(vsd)
# Compute pairwise correlation values
vsd_cor_wt <- cor(vsd_mat_wt) 
View(vsd_cor_wt)
# Plot heatmap
# Plot heatmap
pdf(file=paste0(filename, "Clustering_of_Samplesbyvst.pdf"))
pheatmap(vsd_cor_wt)
dev.off()

###dimensionality reduction of samples
#of vst values
pdf(file = paste0(filename,"PCA_ofSamplesby_vsd.pdf"))
plotPCA(vsd, intgroup=c("genotype"))
dev.off()

### RUN DE Analysis
#Run DESeq2
dds <- DESeq(dds)
structure(dds)
res <- results(dds)
res
# The base mean is the mean of normalized counts of all samples, normalizing for sequencing depth.


### MA plot
# What is MA plot?
# 
# 2-dimensional (2D) scatter plot to visualize gene expression datasets
# Visualize and identify gene expression changes from two different conditions (eg. normal vs. treated) in terms of log fold change (M) on Y-axis and log of mean of expression counts of normal and treated samples (A) on X-axis
# Genes with similar expression values in both normal and treated samples will cluster around M=0 value i.e genes expressed with no significant differences in between treatments
# Points away from M=0 line indicates genes with significant expression, For example, gene is upregulated and downregulated if point is above and below M=0 line respectively
# MA plot does not consider statistical measures (P-values or adjusted P-values) and therefore we can not tell genes with statistically significant differences between normal vs. treated from MA plot (Use Volcano plot if you want indicates genes with statistically significant differences)

#Differential expression across conditions.
#MA Plots before log fold change shrinkage
pdf(file = paste0(filename,"MA_plot_before_lfc_shrinkage.pdf"))
plotMA(res, ylim=c(-5,5))
dev.off()

###MA Plots after log fold change shrinkage
## Shrinkage of effect size (LFC estimates) is useful for visualization and ranking of genes. 
resLFC <- lfcShrink(dds, coef = "genotype_TRTED_vs_WT", type = "apeglm")
pdf(file = paste0(filename,"MA_plot_AFTER_lfc_shrinkage.pdf"))
plotMA(resLFC, ylim=c(-5, 5))
dev.off()

###Filter significant genes for alpha=0.05 and fold change 1.25
#Results: number of genes before filtering
summary(res)
class(res)
##Results: number of genes before filtering with > 1.25 fold change
res_fc_1.25 <- results(dds, lfcThreshold = log2(1.25)) # log2(1.25))
summary(res_fc_1.25)

#Results: number of genes after filtering: > 1.25 fold change & padj < 0.05
res_fc_1.25_alpha_0.05 <- results(dds, lfcThreshold=log2(1.25), alpha=0.05)
summary(res_fc_1.25_alpha_0.05)
resSig <- res_fc_1.25_alpha_0.05



### Heatmap of log FC of significant genes across samples
df <- as.data.frame(colData(dds)[,c("genotype")])
colnames(df) <- c("genotype")
rownames(df) <- colnames(vsd)

mat = assay(vsd)[ head(order(res$padj), 20), ] # select the top 20 genes with the lowest padj
mat = mat - rowMeans(mat)
pdf(file = paste0(filename, "Top20Genes_Samplewise_heatmap.pdf"))
pheatmap(mat, cluster_rows=TRUE, show_rownames=TRUE,
         cluster_cols=TRUE, annotation_col=df)
dev.off()

###plotting  Expression plots for top 20 genes across all samples
resSig
resSigOrdered <- resSig[order(resSig$padj),][1:20,]    #1:20 for larger sig genes
exp_plot_data <- as.data.frame(assay(vsd)[rownames(resSigOrdered), ])
colnames(exp_plot_data) <- coldata[, "genotype"]


class(exp_plot_data)


library(reshape2)

exp_plot_data_matrix <- cbind(ID=rownames(exp_plot_data),  exp_plot_data)
exp_plot_data_matrix <- melt(exp_plot_data_matrix) #melt reshaped broad df into narrow one
pdf(file = paste0(filename, "Top20_ExpressionPlot_geneSymbols.pdf"))
ggplot(data = exp_plot_data_matrix, mapping = aes(x = ID, y = value, color = variable)) +
  geom_jitter(alpha = 0.3, width = 0.25) +
  scale_y_log10() +
  theme(axis.text.x = element_text(colour = "grey20", size = 8, angle = 60, hjust = 1.0, vjust = 1.0),
        axis.text.y = element_text(colour = "grey20", size = 12),
        text = element_text(size = 16)) + 
  #scale_x_discrete(labels=x_label) + 
  xlab("Gene Symbol") +
  ylab("Normalized expression value")
dev.off()


res_nona <- na.omit(res)
dim(res_nona)
write.table(res_nona, file="DE_genes_all_UUO.txt", sep="\t")

