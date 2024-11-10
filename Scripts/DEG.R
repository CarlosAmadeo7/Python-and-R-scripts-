## ----setup, include=FALSE------------------------------------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE,warning = FALSE, message=FALSE,
                      fig.height = 4, fig.width = 8)


## ----Obtaining the counts from the bam files-----------------------------------------------------------------------------
## Loading libraries

#if (!requireNamespace("BiocManager", quietly = TRUE))
    #install.packages("BiocManager")
#BiocManager::install("DESeq2")

library(systemPipeR)
library(GenomicAlignments)
library(GenomicFeatures)
library(BiocParallel)
library(ggplot2)
library("limma")
library(DESeq2)
library("tidyverse")
library("apeglm")
#library("ggfortify")
library("ggplot2")

# Loading the sqlite file 
txdb <- loadDb("C:/Users/amade/Documents/UofSC/FALL 2024/Genomic data science/STAT718_Homework5/data/tair10.sqlite")
# Loading the features 
eByg <- exonsBy(txdb, by = c("gene"))

# Importing the BAM files 
outpaths <- list.files('C:/Users/amade/Documents/UofSC/FALL 2024/Genomic data science/STAT718_Homework5/results/hisat2_mapping/', pattern='sorted.bam$',full.names=TRUE)
bfl <- BamFileList(outpaths, yieldSize = 50000, index = character())


## ----Obtaining the counts with summarizeoverlaps and the loading metadata from the counts--------------------------------
unstranded <- summarizeOverlaps(eByg, bfl, mode="Union",
                                ignore.strand=TRUE,inter.feature=FALSE, 
                                singleEnd=FALSE)  # the single end is false since this data is paired-end data
unstranded <- assays(unstranded)$counts
unstranded[1:5,]


## ----Uploading the sampleData/metadata-----------------------------------------------------------------------------------
url <- "https://people.stat.sc.edu/hoyen/STAT718/Data/targets.txt"
sampleData <- read.delim(file = url, comment = "#")
sampleData

## ------------------------------------------------------------------------------------------------------------------------
## Filtering low expressed genes 
dim(unstranded)  ## 145 genes and 18 samples
keep <- rowSums(unstranded) > 5  # arbitrary picking of genes 
filtCounts <- unstranded[keep,]
filtCounts[1:3,]   ## Filter counts, new dim is 121 genes and 18 columns 

# Deleting the .sorted bam names 
colnames(filtCounts)<-sub(".sorted.bam"," ",colnames(filtCounts))

## Removing unkwanted spaces so we can match with counts.
sampleData$SampleName<-trimws(tolower(sampleData$SampleName))
# Removing spaces in the colnames of the counts 
colnames(filtCounts) <- trimws(tolower(colnames(filtCounts)))
## Matching the sampleData with the filtCounts
sampleData<-sampleData[match(colnames(filtCounts),sampleData$SampleName),]



## ----PCA, Differentail expression analysis-------------------------------------------------------------------------------

### Adding new groups of treatment and time, to perform DE
## Treatment column
sampleData$Treatment <- sub("^([A-Za-z]+)\\..*", "\\1", sampleData$SampleLong)
### Time, extracting them from the SampleLong
sampleData$Time <- sub(".*\\.([0-9]+h)\\..*", "\\1", sampleData$SampleLong)

# Re doing the sample features to match factors and not characters 
sampleData$Factor <- as.factor(sampleData$Factor)
sampleData$SampleLong <- as.factor(sampleData$SampleLong)
sampleData$SampleName<-as.factor(sampleData$SampleName)
sampleData$Treatment<-as.factor(sampleData$Treatment)
sampleData$Time <-as.factor(sampleData$Time )

# Creating the Deseq2 object with the Factor columns
dds<- DESeqDataSetFromMatrix(countData = filtCounts,
                             colData = sampleData, design = ~  Treatment + Time  )
# Performing DESEq2 analysis
dds<- DESeq(dds)

### Performing PCA analysis of the dds
### Normalizing the counts with vst
vst<-varianceStabilizingTransformation(dds, blind = TRUE)
vst_counts<-assay(vst)
## PCA based on the Treatment effect
pcaData <- plotPCA(vst, intgroup = "Treatment", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color = Treatment, shape = Treatment)) +
    geom_point(size = 4) +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    ggtitle("PCA of normalized counts based on treatment effect") +
    theme_minimal()

## PCA based on the time effect
pcaData <- plotPCA(vst, intgroup = "Time", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color = Time, shape = Time)) +
    geom_point(size = 4) +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    ggtitle("PCA of normalized counts based on time effect") +
    theme_minimal()



## ----Differential expression genes---------------------------------------------------------------------------------------
## Identify DE genes between treatment conditions: Mock and Vir.
## Results of the dds
#resultsNames(dds)

## Performing Differential expression analysis and Wald Test
### Between Mock and Avr
contrast.trt <- c("Treatment", "Mock", "Avr")
res.trt <- DESeq2::results(dds, contrast = contrast.trt)
head(res.trt)
dim(res.trt)   ##3 121 genes and 6 columns 

### plotMA shows no interesting data since the number of genes is really low 
#plotMA(res.trt, main="DESeq2 trt vs untrt")
 

## ------------------------------------------------------------------------------------------------------------------------
### Perform shrunken foldchange. Are the results similar? 
# Shrink log2 fold changes
resLFC <- lfcShrink(dds, coef = 2, type = "apeglm")
head(resLFC)
dim(resLFC)


## ------------------------------------------------------------------------------------------------------------------------
# Use Adjust p-value=0.05 as the cutoff and identify significant DE genes.
## i am using the p value since the padj is not that relevant because the dataset is really small 
DEG<- res.trt[which(res.trt$pvalue < 0.05),]
DEG_sh<- resLFC[which(resLFC$pvalue < 0.05), ]
dim(DEG)  ## 10 genes were significant 
dim(DEG_sh)  ## 10 genes were significant

## ------------------------------------------------------------------------------------------------------------------------
### Create volcano plots using the unshrunken and shrunken estimates and the corresponding p values.
#Label significant DE genes labeled in different color then non-DE genes. 
#jpeg("Volcano.png", width=1000, height=500)
padj2<-res.trt[,6]
padj3<-resLFC[,5]

par(mfrow=c(1,2))
plot(res.trt[,2], -log(padj2,10), pch=16, main="Unshrunk data", xlab="log2FoldChange", ylab="-log10(p)")
plot(resLFC[,2], -log(padj3,10), pch=16, main="shrunk data ", xlab="log2FoldChange", ylab="-log10(p)")
abline(h=3, lty=2, col=2)
abline(v=2, lty=2, col=2)
abline(v=-2, lty=2, col=2)



## ------------------------------------------------------------------------------------------------------------------------
## Installing libraries 
#if (!require("BiocManager", quietly = TRUE))
   # install.packages("BiocManager")

#BiocManager::install("AnnotationDbi")
## Generate a heatmap with the shrunk data 
library("pheatmap")
#library("colorRamps")
#library(RColorBrewer)
library(AnnotationDbi)
library(org.Hs.eg.db)

### Normalizing the data
ndata<-counts(dds, normalized=T)
# Identify significant genes (adjusted p-value < 0.05)
whsel <- which(resLFC[,4] < 0.05)
sel <- rownames(resLFC)[whsel]
# Match significant genes with normalized data
matchsel <- match(sel, rownames(ndata))
selcnt <- ndata[matchsel, ]
# Create a heatmap of the significant genes
pheatmap(selcnt, 
         scale = "row", 
         clustering_distance_rows = "euclidean", 
         clustering_distance_cols = "euclidean", 
         clustering_method = "complete", 
         show_rownames = TRUE, 
         show_colnames = TRUE, 
         main = "Heatmap of Significant DE Genes")


