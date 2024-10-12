# Libraries 
library(rtracklayer)
library(readr)
library(tximport)
library(rtracklayer)
library(readr)
library(dplyr)

####################### Analysis of the transcripts counts predicted by the paper

# Uploading the transcript dataset from pangenomes'graph results.
## Change the filepath once you have the file to yours.
file_path_1= "C:/Users/amade/Documents/UofSC/Bioinformatics/Pangenome_project/Analysis/trasncripts_rpvg.txt"
transcript_counts<- read.table(file_path_1,sep="\t",header = T)
head(transcript_counts)
## Getting just the Names, Read_COunts, and TPMs for mapping the transcripts
transcript_counts_1<-transcript_counts[,c("Name","ReadCount","TPM")]
head(transcript_counts_1)


# Uploading the human gtf file
gtf <- import("C:/Users/amade/Documents/UofSC/Bioinformatics/Pangenome_project/gencode.v46.annotation.gtf")
head(gtf)
##Extracting the columns to map transcripts IDs (ENST) to gene IDs(ENGS), seqnames, start, end,strand, transcript_id, gene_id
transcript_entries <- gtf[gtf$type == "transcript"]
transcript_gene_mapping <- data.frame(
  ensembl_transcript_id = transcript_entries$transcript_id,
  ensembl_gene_id = transcript_entries$gene_id,
  gene_name = transcript_entries$gene_name
)
head(transcript_gene_mapping)


# Merging transcript-level data (ENST) to gene (ENSG) file with the GTF mapping
merged_data <- merge(transcript_counts_1, transcript_gene_mapping, 
                     by.x = "Name", by.y = "ensembl_transcript_id")
str(merged_data)
head(merged_data)

# Summing up the ReadCounts by gene (gene-level counts), Obtaing the gene level data from rpvg method 
gene_counts <- aggregate(ReadCount ~ ensembl_gene_id + gene_name, data = merged_data, sum)
str(gene_counts)


### Visualization 
## Boxplot
boxplot(gene_counts$ReadCount,
        main = "Predicted Gene Counts from pangenome graph",
        ylab = "Read Counts",
        col = "lightgreen")
## top 10 most expressed genes 
top_genes <- gene_counts[order(gene_counts$ReadCount, decreasing = TRUE), ][1:10, ]
barplot(top_genes$ReadCount, 
        names.arg = top_genes$gene_name,
        las = 2,  # Makes the gene names vertical
        col = "skyblue", 
        main = "Top 10 Expressed Genes")





################## Analysis of the counts shown by the paper
#### Change the path 
file_path_2= "C:/Users/amade/Documents/UofSC/Bioinformatics/Pangenome_project/Analysis/GSE117823_All_Counts.txt"
paper_counts<- read.table(file_path_2,sep="\t",header = T)
head(paper_counts)
### Getting the first column 
paper_counts_first<-paper_counts[,c("GeneID","G.1")]

## Visualization for paper counts
### Box plot 
boxplot(paper_counts_first$G.1,
        main = "Predicted paper gene counts",
        ylab = "Read Counts",
        col = "lightgreen")
### Top 10 most expressed genes 
top_genes_paper <- paper_counts_first[order(paper_counts_first$G.1, decreasing = TRUE), ][1:10, ]
barplot(top_genes_paper$G.1, 
        names.arg = top_genes_paper$GeneID,
        las = 2,  # Makes the gene names vertical
        col = "skyblue", 
        main = "Top 10 Paper Expressed Genes")
