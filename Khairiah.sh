for f in $(ls /scratch/dragon/intel/velucha/Heribert/Khairiah/1inputseq/*.gz | sed 's/[1-2]_001.fastq.gz//' | sort -u)
do
filename=$(basename $f)
echo "java -jar \$TRIMMOMATIC_JAR PE -phred33 -threads 20 -summary ${filename}1.summary ${f}1_001.fastq.gz ${f}2_001.fastq.gz  ${filename}1_paired.fq.gz ${filename}1_unpaired.fq.gz ${filename}2_paired.fq.gz ${filename}2_unpaired.fq.gz ILLUMINACLIP:\${adapters}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 " 
  done
###########tophat ##########

for f in $(ls /scratch/dragon/intel/velucha/Heribert/Khairiah/1inputseq/*.gz | sed 's/[1-2]_001.fastq.gz//' | sort -u)
do
filename=$(basename $f)
echo "tophat2 --output-dir ${filename}_dir --no-gtf-juncs -p 20 -g 1 -G ../../TAIR10_chrCM.gtf ../../Bowtie_build/TAIR10_chrCM  ../1Trim/${filename}1_paired.fq.gz ../1Trim/${filename}2_paired.fq.gz,../1Trim/${filename}1_unpaired.fq.gz,../1Trim/${filename}2_unpaired.fq.gz"
  done

########### star index ############



#!/bin/bash
#SBATCH --job-name=mapping
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=5:30:00
#SBATCH --mem=100gb
#SBATCH --mail-type=end
#SBATCH --cpus-per-task=40
##SBATCH --constraint=intel 
##SBATCH --qos=ibex-c2024
##SBATCH --constraint=[cpu_intel_gold_6148]

module load star
mkdir star_index

STAR --runThreadN 40 \
--runMode genomeGenerate \
--genomeDir star_index \
--genomeFastaFiles Arabidopsis_thaliana.TAIR10.dna.toplevel.fa \
--sjdbGTFfile Arabidopsis_thaliana.TAIR10.45.gtf \
--sjdbOverhang 150



################ star map ############
for f in $(ls /scratch/dragon/intel/velucha/Heribert/Khairiah/1inputseq/*.gz | sed 's/[1-2]_001.fastq.gz//' | sort -u)
do
filename=$(basename $f)
echo "STAR --runThreadN 40 --sjdbOverhang 150 --genomeDir ../star_index --sjdbGTFtagExonParentGene gene_id --sjdbGTFfile Arabidopsis_thaliana.TAIR10.45.gtf --outSAMtype BAM Unsorted SortedByCoordinate --outReadsUnmapped Fastx --quantMode GeneCounts --outFileNamePrefix ${f}_star_alignment_ --readFilesIn ${f}1_001.fastq.gz ${f}2_001.fastq.gz --readFilesCommand zcat"
  done

########################################## DESEQ2 ##################
ssh -X velucha@ilogin.ibex.kaust.edu.sa
#srun --time=5:00:00 --mem=16G --nodes=1 -c 16 --pty bash -l
git clone https://github.com/CBC-KAUST/diffkaust2019.git
cd diffkaust2019
module load R/3.6.0/gnu-6.4.0
module load RStudio_Desktop/1.1.383
rstudio &

setwd("/scratch/dragon/intel/velucha/Heribert/Khairiah/1inputseq")
raw_data <- read.table("All.tab", sep="\t", header = TRUE, row.names = 1)

library("DESeq2")

################ DE: pairwise comparison diff expression sample1
diff1 <- c('shootMock1', 'shootMock2', 'shootMock3', 'shootSA1901', 'shootSA1902', 'shootSA1903')
selected_data1 <- raw_data[, diff1]
selected_data1 <- as.matrix(selected_data1)
genotype1 <- c(rep("shootMock", 3),rep("shootSA190", 3))
coldata1 <- data.frame(genotype1)
rownames(coldata1) <- colnames(selected_data1)
dds1 <- DESeqDataSetFromMatrix(countData = selected_data1, colData = coldata1, design = ~genotype1)
keep <- rowSums(counts(dds1)) > 0
dds1 <- dds1[keep,]
dds1$genotype1 <- factor(dds1$genotype1, levels = c("shootMock", "shootSA190"))
dds1 <- estimateSizeFactors(dds1)
normalized_data1 <- counts(dds1, normalized=TRUE)
dds1 <- DESeq(dds1)
res1 <- results(dds1)
res1_nona <- na.omit(res1)
write.table(res1_nona, file="shootMockvsshootSA190.xls", sep="\t")

####################DE: pairwise comparison diff expression sample1
diff1 <- c('RootMock1', 'RootMock2', 'RootMock3', 'RootSA1901', 'RootSA1902', 'RootSA1903')
selected_data1 <- raw_data[, diff1]
selected_data1 <- as.matrix(selected_data1)
genotype1 <- c(rep("RootMock", 3),rep("RootSA190", 3))
coldata1 <- data.frame(genotype1)
rownames(coldata1) <- colnames(selected_data1)
dds1 <- DESeqDataSetFromMatrix(countData = selected_data1, colData = coldata1, design = ~genotype1)
keep <- rowSums(counts(dds1)) > 0
dds1 <- dds1[keep,]
dds1$genotype1 <- factor(dds1$genotype1, levels = c("RootMock", "RootSA190"))
dds1 <- estimateSizeFactors(dds1)
normalized_data1 <- counts(dds1, normalized=TRUE)
dds1 <- DESeq(dds1)
res1 <- results(dds1)
res1_nona <- na.omit(res1)
write.table(res1_nona, file="RootMock1vsRootSA190.xls", sep="\t")


#################DE: pairwise comparison diff expression sample1
diff1 <- c('RootMock1', 'RootMock2', 'RootMock3', 'RootPEG1', 'RootPEG2', 'RootPEG3')
selected_data1 <- raw_data[, diff1]
selected_data1 <- as.matrix(selected_data1)
genotype1 <- c(rep("RootMock", 3),rep("RootPEG", 3))
coldata1 <- data.frame(genotype1)
rownames(coldata1) <- colnames(selected_data1)
dds1 <- DESeqDataSetFromMatrix(countData = selected_data1, colData = coldata1, design = ~genotype1)
keep <- rowSums(counts(dds1)) > 0
dds1 <- dds1[keep,]
dds1$genotype1 <- factor(dds1$genotype1, levels = c("RootMock", "RootPEG"))
dds1 <- estimateSizeFactors(dds1)
normalized_data1 <- counts(dds1, normalized=TRUE)
dds1 <- DESeq(dds1)
res1 <- results(dds1)
res1_nona <- na.omit(res1)
write.table(res1_nona, file="RootMock1vsRootPEG.xls", sep="\t")

#################DE: pairwise comparison diff expression sample1

diff1 <- c('RootSA1901', 'RootSA1902', 'RootSA1903','RootSA190PEG1', 'RootSA190PEG2')
selected_data1 <- raw_data[, diff1]
selected_data1 <- as.matrix(selected_data1)
genotype1 <- c(rep("RootSA190", 3),rep("RootSA190PEG", 2))
coldata1 <- data.frame(genotype1)
rownames(coldata1) <- colnames(selected_data1)
dds1 <- DESeqDataSetFromMatrix(countData = selected_data1, colData = coldata1, design = ~genotype1)
keep <- rowSums(counts(dds1)) > 0
dds1 <- dds1[keep,]
dds1$genotype1 <- factor(dds1$genotype1, levels = c("RootSA190", "RootSA190PEG"))
dds1 <- estimateSizeFactors(dds1)
normalized_data1 <- counts(dds1, normalized=TRUE)
dds1 <- DESeq(dds1)
res1 <- results(dds1)
res1_nona <- na.omit(res1)
write.table(res1_nona, file="RootSA190vsRootSA190PEG.xls", sep="\t")

#############DE: pairwise comparison diff expression sample1
diff1 <- c('RootPEG1', 'RootPEG2', 'RootPEG3', 'RootSA190PEG1', 'RootSA190PEG2')
selected_data1 <- raw_data[, diff1]
selected_data1 <- as.matrix(selected_data1)
genotype1 <- c(rep("RootPEG", 3),rep("RootSA190PEG", 2))
coldata1 <- data.frame(genotype1)
rownames(coldata1) <- colnames(selected_data1)
dds1 <- DESeqDataSetFromMatrix(countData = selected_data1, colData = coldata1, design = ~genotype1)
keep <- rowSums(counts(dds1)) > 0
dds1 <- dds1[keep,]
dds1$genotype1 <- factor(dds1$genotype1, levels = c("RootPEG", "RootSA190PEG"))
dds1 <- estimateSizeFactors(dds1)
normalized_data1 <- counts(dds1, normalized=TRUE)
dds1 <- DESeq(dds1)
res1 <- results(dds1)
res1_nona <- na.omit(res1)
write.table(res1_nona, file="RootPEGvsRootSA190PEG.xls", sep="\t")

########### DE: pairwise comparison diff expression sample1
diff1 <- c('shootMock1', 'shootMock2', 'shootMock3', 'shootSA1901', 'shootSA1902', 'shootSA1903')
diff2 <- c('shootPEG1', 'shootPEG2', 'shootPEG3', 'shootSA190PEG1', 'shootSA190PEG2')
diff3 <- c('shootMock1', 'shootMock2', 'shootMock3', 'shootPEG1', 'shootPEG2', 'shootPEG3')
diff4 <- c('shootSA1901', 'shootSA1902', 'shootSA1903','shootSA190PEG1', 'shootSA190PEG2')

diff5 <- c('RootMock1', 'RootMock2', 'RootMock3', 'RootSA1901', 'RootSA1902', 'RootSA1903')
diff6 <- c('RootPEG1', 'RootPEG2', 'RootPEG3', 'RootSA190PEG1', 'RootSA190PEG2','RootSA190PEG3')
diff7 <- c('RootMock1', 'RootMock2', 'RootMock3', 'RootPEG1', 'RootPEG2', 'RootPEG3')
diff8 <- c('RootSA1901', 'RootSA1902', 'RootSA1903','RootSA190PEG1', 'RootSA190PEG2','RootSA190PEG3')

#####All normalization:

########## merge Annotation to table:function########

library(data.table)
b <- read.table("TAIR10_Ann.xls", sep="\t",header=TRUE)
dt1 <- data.table(b, key = "gid")
files <- list.files(path="/Users/velucha/Khairiah", pattern="*.resv", full.names=TRUE, recursive=FALSE)
lapply(files, function(x) {
    print(x)
    a <- read.table(x, sep="\t",header=TRUE)
    dt2 <- data.table(a, key = "gid")
    c <- dt1[dt2]
    fin <- paste(x, "fin", sep="_")
    write.table(c, fin, sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
})

#####################################

#######for plotting
library("DESeq2")
library(ggplot2)
library("RColorBrewer")

library(biobroom)
library("vsn")
library(apeglm)
library(annotables)
library(biobroom)
library(dplyr)
library(ggrepel)
library("EnhancedVolcano")
library(gProfileR)


diff1 <- c('shootMock1', 'shootMock2', 'shootMock3', 'shootSA1901', 'shootSA1902', 'shootSA1903','shootPEG1', 'shootPEG2', 'shootPEG3', 'shootSA190PEG1', 'shootSA190PEG2', 'RootMock1', 'RootMock2', 'RootMock3', 'RootSA1901', 'RootSA1902', 'RootSA1903','RootPEG1', 'RootPEG2', 'RootPEG3', 'RootSA190PEG1', 'RootSA190PEG2','RootSA190PEG3')
selected_data1 <- raw_data[, diff1]


########## plot1 before normalization

pdf(file = "raw_count_distribution_ofuuosamples.pdf")
#selected_data <- as.matrix(raw_data)
dat_log2 <- stack(as.data.frame(log2(selected_data1)))
ggplot(data = dat_log2, mapping = aes(x = ind, y = values)) +
  #geom_jitter(alpha = 0.3, color = "tomato") +
  geom_boxplot(alpha = 0) +
scale_x_discrete(name ="Conditions")+
scale_y_discrete(name ="Raw count")+
theme(axis.text.x = element_text(angle = 90))
dev.off()

########Plotting Individual replicates after normalization

diff2 <- c('shootMock1', 'shootMock2', 'shootMock3', 'shootSA1901', 'shootSA1902', 'shootSA1903','shootPEG1', 'shootPEG2', 'shootPEG3', 'shootSA190PEG1', 'shootSA190PEG2', 'RootMock1', 'RootMock2', 'RootMock3', 'RootSA1901', 'RootSA1902', 'RootSA1903','RootPEG1', 'RootPEG2', 'RootPEG3', 'RootSA190PEG1', 'RootSA190PEG2','RootSA190PEG3')
selected_data2 <- raw_data[, diff2]

genotype2 <- c(rep("shootMock1",1), rep("shootMock2",1),rep("shootMock3",1), rep("shootSA1901",1), rep("shootSA1902",1), rep("shootSA1903",1),rep("shootPEG1",1), rep("shootPEG2",1), rep("shootPEG3",1), rep("shootSA190PEG1",1), rep("shootSA190PEG2",1), rep("RootMock1",1), rep("RootMock2",1), rep("RootMock3",1), rep("RootSA1901",1), rep("RootSA1902",1), rep("RootSA1903",1),rep("RootPEG1",1), rep("RootPEG2",1), rep("RootPEG3",1), rep("RootSA190PEG1",1), rep("RootSA190PEG2",1),rep("RootSA190PEG3",1))

selected_data2 <- raw_data[, diff2]
coldata2 <- data.frame(genotype2)
rownames(coldata2) <- colnames(selected_data2)

dds2 <- DESeqDataSetFromMatrix(countData = selected_data2,
                              colData = coldata2,
                              design = ~ genotype2)
dds2 <- estimateSizeFactors(dds2)
sizeFactors(dds2)

normalized_data2 <- counts(dds2, normalized=TRUE)
dat_log2 <- stack(as.data.frame(log2(normalized_data2)))
pdf(file = "CountDistribution_after_normalisation4.pdf")
ggplot(data = dat_log2, mapping = aes(x = ind, y = values)) +
geom_boxplot(alpha = 0) +
scale_x_discrete(name ="Conditions")+
scale_y_discrete(name ="log2 Normalized count")+
theme(axis.text.x = element_text(angle = 90))
dev.off()

#############Get normalized table for heat map

diff2 <- c('shootMock1', 'shootMock2', 'shootMock3', 'shootSA1901', 'shootSA1902', 'shootSA1903','shootPEG1', 'shootPEG2', 'shootPEG3', 'shootSA190PEG1', 'shootSA190PEG2', 'RootMock1', 'RootMock2', 'RootMock3', 'RootSA1901', 'RootSA1902', 'RootSA1903','RootPEG1', 'RootPEG2', 'RootPEG3', 'RootSA190PEG1', 'RootSA190PEG2','RootSA190PEG3')
selected_data2 <- raw_data[, diff2]

genotype2 <- c(rep("shootMock", 3),rep("shootSA190", 3), rep("shootPEG", 3),rep("shootSA190PEG", 2),rep("RootMock", 3),rep("RootSA190", 3), rep("RootPEG", 3),rep("RootSA190PEG", 3))

coldata2 <- data.frame(genotype2)
rownames(coldata2) <- colnames(selected_data2)

dds2 <- DESeqDataSetFromMatrix(countData = selected_data2,
                              colData = coldata2,
                              design = ~ genotype2)

keep <- rowSums(counts(dds2)) > 0
dds2 <- dds2[keep,]
summary(dds2)
dds2$genotype2 <- factor(dds2$genotype2, levels = c("shootMock","shootSA190","shootPEG","shootSA190PEG","RootMock","RootSA190","RootPEG", "RootSA190PEG"))
dds2 <- estimateSizeFactors(dds2)
normalized_data2 <- counts(dds2, normalized=TRUE)

table<-counts(ddsColl)

write.table(normalized_data2, file="normalizedwithRep_counts.txt", sep="\t", quote=F)
write.table(table, file="normalized_counts.txt", sep="\t", quote=F)

################################### PCA plot after transformation of raw count #######

diff3 <- c('shootMock1', 'shootMock2', 'shootMock3', 'shootSA1901', 'shootSA1902', 'shootSA1903','shootPEG1', 'shootPEG2', 'shootPEG3', 'shootSA190PEG1', 'shootSA190PEG2', 'RootMock1', 'RootMock2', 'RootMock3', 'RootSA1901', 'RootSA1902', 'RootSA1903','RootPEG1', 'RootPEG2', 'RootPEG3', 'RootSA190PEG1', 'RootSA190PEG2','RootSA190PEG3')
selected_data3 <- raw_data[, diff3]
genotype3 <- c(rep("shootMock", 3),rep("shootSA190", 3), rep("shootPEG", 3),rep("shootSA190PEG", 2),rep("RootMock", 3),rep("RootSA190", 3), rep("RootPEG", 3),rep("RootSA190PEG", 3))
coldata3 <- data.frame(genotype3)
rownames(coldata3) <- colnames(selected_data3)
dds3 <- DESeqDataSetFromMatrix(countData = selected_data3,
                              colData = coldata3,
                              design = ~ genotype3)
vsd <- vst(dds3, blind=FALSE)

pdf(file="PCA_of_Samplesbyvstnew.pdf")
plotPCA(vsd, intgroup=c("genotype3"))
dev.off()


###############Distance plot

library("RColorBrewer")
library("pheatmap")
diff3 <- c('shootMock1', 'shootMock2', 'shootMock3', 'shootSA1901', 'shootSA1902', 'shootSA1903','shootPEG1', 'shootPEG2', 'shootPEG3', 'shootSA190PEG1', 'shootSA190PEG2', 'RootMock1', 'RootMock2', 'RootMock3', 'RootSA1901', 'RootSA1902', 'RootSA1903','RootPEG1', 'RootPEG2', 'RootPEG3', 'RootSA190PEG1', 'RootSA190PEG2','RootSA190PEG3')
selected_data3 <- raw_data[, diff3]
genotype3 <- c(rep("shootMock", 3),rep("shootSA190", 3), rep("shootPEG", 3),rep("shootSA190PEG", 2),rep("RootMock", 3),rep("RootSA190", 3), rep("RootPEG", 3),rep("RootSA190PEG", 3))
coldata3 <- data.frame(genotype3)
rownames(coldata3) <- colnames(selected_data3)
dds3 <- DESeqDataSetFromMatrix(countData = selected_data3,
                              colData = coldata3,
                              design = ~ genotype3)
vsd <- vst(dds3, blind=FALSE)

sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
colnames(sampleDistMatrix) <- colnames(vsd)
rownames(sampleDistMatrix) <- colnames(vsd)
pdf(file="Hierarchical_dist4.pdf")
pheatmap(sampleDistMatrix, clustering_distance_cols=sampleDists, clustering_distance_rows=sampleDists, col=colors)
dev.off()



######## plot2  After normalization

genotype <- c(rep("shootMock", 3),rep("shootSA190", 3), rep("shootPEG", 3),rep("shootSA190PEG", 2),rep("RootMock", 3),rep("RootSA190", 3), rep("RootPEG", 3),rep("RootSA190PEG", 3))
coldata <- data.frame(genotype)
rownames(coldata) <- colnames(selected_data1)
dds <- DESeqDataSetFromMatrix(countData = selected_data1,
                              colData = coldata,
                              design = ~ genotype)
keep <- rowSums(counts(dds)) > 0 
dds <- dds[keep,]


dds <- estimateSizeFactors(dds)
sizeFactors(dds)
View(counts(dds))

##Normalize the data for library size as normalization factor
normalized_data <- counts(dds, normalized=TRUE)
normalized_data
View(normalized_data)
pdf(file = paste0(filename,"CountDistribution_after_normalisation.pdf"))
ggplot(data = dat_log2, mapping = aes(x = ind, y = values)) +
  geom_jitter(alpha = 0.3, color = "forestgreen") + 
  geom_boxplot(alpha = 0) 
  #geom_jitter(alpha = 0.3, color = "darkred") +  #color = gold1, tomato, forestgreeen
dev.off()
#################################################################Unwanted Cuffdiff execution

cuffdiff -p 10 -b TAIR10_chrCM.fa -L  Wt_mock,wt3_flg22 --output-dir out ../1inputseq/M_19_2643_21_AD001_L001_R1_001.fastq.gz	../1inputseq/M_19_2645_29_AD005_L001_R1_001.fastq.gz	../1inputseq/M_19_2647_45_AD019_L001_R1_001.fastq.gz --library-norm-method quartile
cuffdiff -p 10 -b TAIR10_chrCM.fa -L  Wt_mock,wt3_flg22 --output-dir out ../1inputseq/M_19_2644_22_AD002_L001_R1_001.fastq.gz	../1inputseq/M_19_2646_30_AD006_L001_R1_001.fastq.gz	../1inputseq/M_19_2648_46_AD020_L001_R1_001.fastq.gz --library-norm-method quartile 
cuffdiff -p 10 -b TAIR10_chrCM.fa -L  Wt_mock,wt3_flg22 --output-dir out ../1inputseq/M_19_2649_55_AD003_L002_R1_001.fastq.gz	../1inputseq/M_19_2652_47_AD012_L002_R1_001.fastq.gz 
	
cuffdiff -p 10 -b TAIR10_chrCM.fa -L  Wt_mock,wt3_flg22 --output-dir out ../1inputseq/M_19_2650_56_AD004_L002_R1_001.fastq.gz	../1inputseq/M_19_2651_8_AD011_L002_R1_001.fastq.gz	../1inputseq/M_19_2653_16_AD025_L002_R1_001.fastq.gz
cuffdiff -p 10 -b TAIR10_chrCM.fa -L  Wt_mock,wt3_flg22 --output-dir out ../1inputseq/M_19_2654_65_AD009_L003_R1_001.fastq.gz	../1inputseq/M_19_2656_17_AD013_L003_R1_001.fastq.gz	../1inputseq/M_19_2658_41_AD015_L003_R1_001.fastq.gz
cuffdiff -p 10 -b TAIR10_chrCM.fa -L  Wt_mock,wt3_flg22 --output-dir out ../1inputseq/M_19_2655_66_AD010_L003_R1_001.fastq.gz	../1inputseq/M_19_2657_18_AD014_L003_R1_001.fastq.gz	../1inputseq/M_19_2659_42_AD016_L003_R1_001.fastq.gz
cuffdiff -p 10 -b TAIR10_chrCM.fa -L  Wt_mock,wt3_flg22 --output-dir out ../1inputseq/M_19_2660_3_AD021_L004_R1_001.fastq.gz	../1inputseq/M_19_2662_43_AD027_L004_R1_001.fastq.gz	../1inputseq/M_19_2664_61_AD007_L004_R1_001.fastq.gz
cuffdiff -p 10 -b TAIR10_chrCM.fa -L  Wt_mock,wt3_flg22 --output-dir out ../1inputseq/M_19_2661_4_AD022_L004_R1_001.fastq.gz	../1inputseq/M_19_2663_44_AD018_L004_R1_001.fastq.gz	../1inputseq/M_19_2665_62_AD008_L004_R1_001.fastq.gz
