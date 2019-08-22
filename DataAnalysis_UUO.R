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

library("DESeq2")
library(ggplot2)
library("RColorBrewer")
library("pheatmap")
library("vsn")
library(apeglm)
library(annotables)
library(biobroom)
library(dplyr)
library(ggrepel)
library("EnhancedVolcano")
raw_data <- read.csv("fibrosis_counts.csv", header=TRUE, row.names = 1)
#raw_data <- read.table("AllSamples_Expression.txt", sep="\t", header = TRUE, row.names = 1)
class(raw_data)
dim(raw_data)

head(raw_data)

#rownames(raw_data) <- raw_data[, 1]

head(raw_data)
##inhouse processed
#selected_data <- raw_data[, c( 11:14,4:7)] ## to get wt_UUO and Smoc_UUO
#selected_data <- raw_data[, c( 8:10,1:3)] ## to get wt_norm and Smoc_norm
### count table from paper
selected_data <- raw_data[, c(1:3, 8:10)] ### to get wt_norm and smoc_norm
selected_data <- as.matrix(selected_data)
dim(selected_data)
class(selected_data)
head(selected_data)

###Before normalisation counts plot
dat_log2 <- stack(as.data.frame(log2(selected_data)))

ggplot(data = dat_log2, mapping = aes(x = ind, y = values)) +
  geom_jitter(alpha = 0.3, color = "tomato") +
  geom_boxplot(alpha = 0)
##Expt meta data ----
##Create the DESeq2 data object

#genotype <- c(rep("WT", 4),rep("TRTED", 4))  ### for UUO
genotype <- c(rep("WT", 3),rep("TRTED", 3))
coldata <- data.frame(genotype)
rownames(coldata) <- colnames(selected_data)
coldata

dds <- DESeqDataSetFromMatrix(countData = selected_data,
                              colData = coldata,
                              design = ~ genotype)
dds
summary(dds)

###
## 0.2 Pre-filtering ----
keep <- rowSums(counts(dds)) >= 0 # 0 means no pre-filtering
dds <- dds[keep,]

dds$genotype <- factor(dds$genotype, levels = c("WT", "TRTED"))
dds$genotype


###Normalization ----
dds <- estimateSizeFactors(dds)
normalized_data <- counts(dds, normalized=TRUE)
normalized_data

### visualize normalized data
dat_log2 <- stack(as.data.frame(log2(normalized_data)))
###Stacking vectors concatenates multiple vectors into a single vector along with a factor indicating where 
#each observation originated. Unstacking reverses this operation. (from utils package)

ggplot(data = dat_log2, mapping = aes(x = ind, y = values)) +
  geom_jitter(alpha = 0.3, color = "forestgreen") + 
  geom_boxplot(alpha = 0) 
  #geom_jitter(alpha = 0.3, color = "darkred") +  #color = gold1, tomato, forestgreeen
  
#### Clustering of samples ----
#Unsupervised clustering of samples
#Result: Heatmap of all sample distance & PCA plot

rld <- rlog(dds, blind=FALSE) # Extracting transformed values

vsd <- vst(dds, blind=FALSE)

sampleDists <- dist(t(assay(rld)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- rld$genotype
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

### Variance stabilised transformed
sampleDists <- dist(t(assay(vsd)))

sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- rld$genotype
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
###dimensionality reduction of samples
# of rlog values
plotPCA(rld, intgroup=c("genotype"))
#of vst values
plotPCA(vsd, intgroup=c("genotype"))

### RUN DE Analysis
#Run DESeq2
dds <- DESeq(dds)
structure(dds)
resultsNames(dds)
res <- results(dds)
res

### Mean - Variance relationship
#Mean-variance relationship of wild type and smoc samples
#Result: Mean Vs variance plot
ntd <- normTransform(dds)
meanSdPlot(assay(ntd))
meanSdPlot(assay(vsd))
meanSdPlot(assay(rld))
## Now mean-variance plot after DE using DESeq2
plotDispEsts(dds)

### MA plot
#Differential expression across conditions.
#MA Plots before log fold change shrinkage
plotMA(res, ylim=c(-5,5))


###MA Plots after log fold change shrinkage
resLFC <- lfcShrink(dds, coef = "genotype_TRTED_vs_WT", type = "apeglm")
plotMA(resLFC, ylim=c(-4, 4))

###Filter significant genes for alpha=0.05 and fold change 1.25
#Results: number of genes before filtering
summary(res)

##Results: number of genes before filtering with > 1.25 fold change
res_fc_1.25 <- results(dds, lfcThreshold = log2(1.25)) # log2(1.25))
summary(res_fc_1.25)

#Results: number of genes after filtering: > 1.25 fold change & padj < 0.05
res_fc_1.25_alpha_0.05 <- results(dds, lfcThreshold=log2(1.25), alpha=0.05)
summary(res_fc_1.25_alpha_0.05)


###Add ensembl annotation using annotables R package to filtered genes
#Result: table of filtered genes with annotation
resSig <- subset(res_fc_1.25_alpha_0.05, padj < 0.05)
summary(resSig)
#, lfcThreshold=log2(1.25), alpha=0.05)
resSig
####Annotation ## not needed if genes are symbols (NOT ensembleid)
library(annotables)


resSig_tidy <- tidy.DESeqResults(resSig)

resSig_tidy %>% 
  dplyr::arrange(p.adjusted) %>% 
  dplyr::inner_join(grcm38, by = c("gene" = "ensgene")) %>% 
  dplyr::select(gene, estimate, p.adjusted, symbol) %>% 
  knitr::kable(.)

# with description

resSig_tidy %>% 
  dplyr::arrange(p.adjusted) %>% 
  dplyr::inner_join(grcm38, by = c("gene" = "ensgene")) %>% 
  dplyr::select(gene, estimate, p.adjusted, symbol, description) %>% 
  knitr::kable(.)

####Get subset of normalized significant genes with padj < 0.05
#Result: table of those genes and heatmap of them across conditions

df <- as.data.frame(colData(dds)[,c("genotype")])
colnames(df) <- c("genotype")
rownames(df) <- colnames(vsd)

mat = assay(vsd)[ head(order(res$padj), 20), ] # select the top 20 genes with the lowest padj
mat = mat - rowMeans(mat)

pheatmap(mat, cluster_rows=TRUE, show_rownames=TRUE,
         cluster_cols=TRUE, annotation_col=df)

###Volcano plot of significant genes from step 10 Result: Volcano plot
#Result: Volcano plot on all genes
EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'pvalue',
                xlim = c(-10, 10),
                selectLab = c(""))
###Result: Volcano plot on significant genes
EnhancedVolcano(resSig,
                lab = rownames(resSig),
                x = 'log2FoldChange',
                y = 'pvalue',
                xlim = c(-10, 10),
                selectLab = c(""))

####Expression plot of top 20 genes (least padj) across samples
#Result: Expression plot
resSig_tidy_top_20 <- resSig_tidy %>% 
  dplyr::arrange(p.adjusted) %>% 
  head(20) %>% 
  dplyr::inner_join(grcm38, by = c("gene" = "ensgene")) %>% 
  dplyr::select(gene, estimate, p.adjusted, symbol)

x_label <- c(as.matrix(resSig_tidy_top_20[, "symbol"]))
x_label

# resSig_tidy_top_20 <- as.data.frame(res) %>% tibble::rownames_to_column("gene")%>%
#   dplyr::arrange(padj) %>% 
#   head(19) %>% dplyr::select(gene, log2FoldChange, padj)
#   #dplyr::inner_join(grcm38, by = c("gene" = "ensgene")) %>% 
#   #dplyr::select(gene, estimate, p.adjusted, symbol)
#x_label <- c(as.matrix(resSig_tidy_top_20[, "gene"]))
#x_label

###plotting top 20 genes across all samples
resSig
resSigOrdered <- resSig[order(resSig$padj),][1:20,]    #1:20 for larger sig genes
exp_plot_data <- as.data.frame(assay(vsd)[rownames(resSigOrdered), ])
colnames(exp_plot_data) <- coldata[, "genotype"]


class(exp_plot_data)


library(reshape2)

exp_plot_data_matrix <- cbind(ID=rownames(exp_plot_data),  exp_plot_data)
exp_plot_data_matrix <- melt(exp_plot_data_matrix)

ggplot(data = exp_plot_data_matrix, mapping = aes(x = ID, y = value, color = variable)) +
  geom_jitter(alpha = 0.3, width = 0.25) +
  scale_y_log10() +
  theme(axis.text.x = element_text(colour = "grey20", size = 8, angle = 60, hjust = 1.0, vjust = 1.0),
        axis.text.y = element_text(colour = "grey20", size = 12),
        text = element_text(size = 16)) + 
  scale_x_discrete(labels=x_label) + 
  xlab("Gene Symbol") +
  ylab("Normalized expression value")

ggplot(data = exp_plot_data_matrix, mapping = aes(x = ID, y = value, color = variable)) +
  geom_jitter(alpha = 0.3, width = 0.25) +
  scale_y_log10() +
  theme(axis.text.x = element_text(colour = "grey20", size = 8, angle = 60, hjust = 1.0, vjust = 1.0),
        axis.text.y = element_text(colour = "grey20", size = 12),
        text = element_text(size = 16)) +
  xlab("Ensembl ID") +
  ylab("Normalized expression value")


res_nona <- na.omit(res)
dim(res_nona)
write.table(res_nona, file="DE_genes_all.txt", sep="\t")
