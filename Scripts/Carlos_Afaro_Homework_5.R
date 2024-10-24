## ----setup, include=FALSE--------------------------------------------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE,warning = FALSE, message=FALSE,
                      fig.height = 4, fig.width = 8)


## ----Loading libraries, dataset and obtianing the sense, antisense and unstranded counts-----------------------------------------
## Loading the libraries 
library(systemPipeR)
library(GenomicAlignments)
library(GenomicFeatures)
library(BiocParallel)
library(ggplot2)

# Loading the sqlite file 
txdb <- loadDb("C:/Users/amade/Documents/UofSC/FALL 2024/Genomic data science/STAT718_Homework5/data/tair10.sqlite")

# Loading the features 
eByg <- exonsBy(txdb, by = c("gene"))

# Importing the BAM files 
outpaths <- list.files('C:/Users/amade/Documents/UofSC/FALL 2024/Genomic data science/STAT718_Homework5/results/hisat2_mapping/', pattern='sorted.bam$',full.names=TRUE)

bfl <- BamFileList(outpaths, yieldSize = 50000, index = character())

# Generate strand-specific for positive (sense) strand and Strand-specific for negative (antisense) strand.
## Generating the positive strand sense counts 
### Using the singleEND=FALSE (for paired end data ) and strandMode=1 (for postive sense strand)
pos_sense<- summarizeOverlaps(eByg, bfl, mode="Union",
                                ignore.strand=FALSE,inter.feature=FALSE,
                                singleEnd=FALSE,strandMode=1)

## Generating the negative strand sense counts 
### Using the singleEND=FALSE (for paired end data) and strandMode=2 (for anti sense strand)
anti_sense <-summarizeOverlaps(eByg, bfl, mode="Union",
                               ignore.strand=FALSE, inter.feature=FALSE,
                               singleEnd=FALSE, strandMode=2)

## Generating the unstranded counts 
### Using the singleEND=FALSE (for paired end data)
unstranded <- summarizeOverlaps(eByg, bfl, mode="Union",
                                ignore.strand=TRUE,inter.feature=FALSE,
                                singleEnd=FALSE)  #### THIS IS FOR PAIRED END.

### Summarizing data
pos_sense<-assays(pos_sense)$counts
pos_sense[1:3,]

anti_sense<-assays(anti_sense)$counts
anti_sense[1:3,]

unstranded <- assays(unstranded)$counts
unstranded[1:3,]


## ----Assesing similarity or differences------------------------------------------------------------------------------------------
#install.packages("pheatmap")
library(pheatmap)
## Summing the two strand-specific read count table and compare to the unstranded count table. Assesing similarity
summing_counts<-pos_sense + anti_sense
## Plotting the absolute diff bettwen the 
difference <- abs(summing_counts - unstranded)
## Finding a correlation of bettwen the sense+antisense AND the unstranded
correlation <- cor(summing_counts, unstranded)
## Visualizing the correlation as a cluster heatmap
pheatmap(correlation,
         main="Correlation bettween the unstranded and totalcounts")
## Plotting a scatter plot to assess similarity bettwen the unstranded and the sum of the sense and antisense data 
plot(summing_counts[,1], unstranded[,1], 
     xlab = "Total Stranded Counts (sense + antisense)", 
     ylab = "Unstranded Counts", 
     main = "Is this dataset unstranded or stranded?")
abline(0, 1, col = "blue")


