######################################################################################################### 
## Libraries
library(readxl)
library(DESeq2)
library(ggplot2)
###########################################################################################################
# Uploading datacounts
file_path <- "C:/Users/amade/Documents/UofSC/Bioinformatics/MPP_analysis/RNA_seq_data_Cabezas/data_cabezas.csv" 
data <- read.csv(file_path, header = T)
str(data)
#### Prepare metadata
col_data <- data.frame(
  condition = factor(c("HSC", "HSC","HSC","HSC",
                       "MPP1", "MPP1", "MPP1",
                       "MPP2", "MPP2", "MPP2",
                       "MPP3", "MPP3", "MPP3",
                       "MPP4", "MPP4","MPP4"
  ))
)
# rownames 
rownames(col_data) <- colnames(data)
### Creating a Deseq object 
dds <- DESeqDataSetFromMatrix(countData = as.matrix(data), colData = col_data, design = ~ condition)
dds
### Filtering lowering expression ( lower expressed genes , less than 10)
dds <- dds[rowSums(counts(dds)) > 10, ]
dds
### Normalization of counts, varianceStabilizingTransformation
vsd <- vst(dds, blind = FALSE)
#### Performing PCA, returns two principal components 
pcaData <- plotPCA(vsd, intgroup = c("condition"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
head(pcaData)
##### Visualization plot 
ggplot(pcaData, aes(PC1, PC2, color = condition)) +
  geom_point(size = 4, alpha = 0.8) +  # Bigger points with slight transparency
  scale_color_brewer(palette = "Set1") +  # Color-blind friendly palette
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_minimal() +  # Clean theme
  theme(
    axis.line = element_blank(),  # Remove all axis lines initially
    axis.line.x.bottom = element_line(color = "black", size = 1),  # Highlight bottom axis
    axis.line.y.left = element_line(color = "black", size = 1),  # Highlight left axis
    axis.title = element_text(size = 14),  # Bigger axis labels
    axis.text = element_text(size = 12),  # Bigger axis tick labels
    legend.title =element_blank(),
    legend.text = element_text(size = 10),
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank()   # Remove minor grid lines
  )
