## ----setup, include=FALSE-----------------------------------------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE,warning = FALSE, message=FALSE,
                      fig.height = 4, fig.width = 8)


## ----Uploading the file-------------------------------------------------------------------------------------------------------
# Uploading the file to RStudio

library(GenomeInfoDb)
human_genome <- read.delim("C:/Users/amade/Downloads/hhh.txt", header=TRUE, sep = '\t')


## ----Finding the lenght of genes in the dataset-------------------------------------------------------------------------------
L<- length(human_genome$name)
print(paste("The number of genes in the dataset is", L))


## ----Number of genes per chromosome-------------------------------------------------------------------------------------------
# Creation of a character vector with all the chromosomes names 
chromosomes<- c(paste0("chr", 1:22), "chrX", "chrY","chrM")
# Looping each one of them 
for (i in chromosomes){
  chr_count <- sum(human_genome$chrom ==i)
  print(paste(i, ":",chr_count))
}


## ----Positive and negative strand dataset-------------------------------------------------------------------------------------
# Creating a dataset" positive strand "
positive_strand <- human_genome[human_genome$ strand =="+",]
# Creating a dataset" negative strand "
negative_strand <- human_genome[human_genome$ strand =="-",]
# Printing the number of rows in each dataset
x <-nrow(positive_strand)
y <-nrow(negative_strand)
print(paste("The number of genes in the + strand is" , x))
print(paste("The number of genes in the - strand is" , y))



## ----Gene Length--------------------------------------------------------------------------------------------------------------
# Creating another column with all the gene length. 
human_genome$gene_lenght<- human_genome$txEnd- human_genome$txStart


## ----Plotting the positive and negative strand--------------------------------------------------------------------------------
# Creating positive and negative strand datasets 
p_strand<- human_genome[human_genome$ strand =="+",]
n_strand<- human_genome[human_genome$ strand =="-",]
# Plotting the distribution difference as histograms 
par(mfrow = c(1, 2))
hist(p_strand$gene_lenght,
     main="Gene Lenght distribution on the positive strand",
     cex.main=0.9,
     xlab="Gene Lenght",
     ylab= "Frequency",
     ylim=c(0,80000),
     xlim=c(0,700000),
     col="red",
     freq=TRUE,
     breaks=50,
     ps=30)
hist(n_strand$gene_lenght,
     main="Gene Lenght distribution on the negative strand",
     cex.main=0.9,
     xlab="Gene Lenght",
     ylab= "Frequency",
     ylim=c(0,80000),
     xlim=c(0,700000),
     col="blue",
     freq=TRUE,
     breaks=50,
     ps=30)



## ----Longest and Shortest gene------------------------------------------------------------------------------------------------
# Finding the shortest and longest gene 
longest_gene <- human_genome[which.max(human_genome$gene_lenght),]
shortest_gene <- human_genome[which.min(human_genome$gene_lenght),]
l_gene <- longest_gene$name
s_gene <- shortest_gene$name
print(paste("The longest gene in the dataset is", l_gene))
print(paste("The shortest gene in the dataset is", s_gene))



## ----Exon counts--------------------------------------------------------------------------------------------------------------
# Getting the positive and negative exon counts 
Positive_exons<- p_strand$exonCount
Negative_exons<- n_strand$exonCount
# Plotting the differences as histograms
par(mfrow = c(1, 2))
hist(Positive_exons,
     main="Exon count distribution on the positive strand",
     cex.main=0.9,
     xlab="Exon Count",
     ylab= "Frequency",
     ylim=c(0,80000),
     xlim = c(0,150),
     col="green",
     freq=TRUE,
     breaks=30,
     )
hist(Negative_exons,
     main="Exon count distribution on the negative strand",
     cex.main=0.9,
     xlab="Exon Count",
     ylab= "Frequency",
     ylim=c(0,80000),
     xlim = c(0,150),
     col="gray",
     freq=TRUE,
     breaks=30,
     )



## ----Correlation in the dataset-----------------------------------------------------------------------------------------------
# Exon and gene length in the whole dataset 
exonCount<-human_genome$exonCount
geneLenght<-human_genome$gene_lenght
plot(exonCount,geneLenght,
     col="red",
     main="Distribution of gene lenght and exon Count ",
     type="p")


## ----Pearson correlation------------------------------------------------------------------------------------------------------
# Pearson correlation of geneLenght and Exoncounts in the dataset 
correlation <- cor(exonCount, geneLenght, use = "complete.obs", method = "pearson")
print(paste("The correlation between exonCount and gene_length in the whole dataset is :", correlation))
#### Pearson correlation of geneLenght and Exoncounts in the positive strand 
ps<-cor(p_strand$exonCount, p_strand$gene_lenght,use = "complete.obs", method = "pearson")
###  Pearson correlation of geneLenght and Exoncounts in the negative strand 
ns<- cor(n_strand$exonCount, n_strand$gene_lenght,use = "complete.obs", method = "pearson")
print(paste("The correlation between exonCount and gene_length in the positive strand is :", ps))
print(paste("The correlation between exonCount and gene_length in the negative strand is : ", ns))



## ----Gene density-------------------------------------------------------------------------------------------------------------
# Installing the package 
#if (!require("BiocManager", quietly = TRUE))
  #install.packages("BiocManager")
# Find out the length of the chromosomes
# Loading the information and getting the sizes of the chromosomes 
hg38_chrom_length<-getChromInfoFromUCSC("hg38")
#print(hg38_chrom_length)
hg38_length<-hg38_chrom_length$size[1:25]
# Find the length of the genes in each chromosomes and store it in a vector
chromosomes<- c(paste0("chr", 1:22), "chrX", "chrY","chrM")
ans <-vector(length=length(chromosomes))
k=1
for (i in chromosomes){
  chrom_data <- human_genome[human_genome$chrom == i,]
  ans[k]= nrow(chrom_data)
  k=k+1
}
# Then convert the length of each chromosome in MB by dividing it by 1e6
ratio_mb <- hg38_length/1e6
# Calculating the gene density by dividing the number of genes per chromosome and the ratio_mb
gene_density <- ans / ratio_mb
# Putting everything in a dataframe for visualization purposes 
gene_density_data <- data.frame(
  chromosome <-chromosomes,
  Number_of_Genes = ans,
  Chromosome_Length_Mb = ratio_mb,
  Gene_Density = gene_density
)
print(gene_density_data)


