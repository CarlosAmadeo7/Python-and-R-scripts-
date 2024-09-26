## ----setup, include=FALSE----------------------------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE,warning = FALSE, message=FALSE,
                      fig.height = 4, fig.width = 8)
## ----Function to create the 10000 length sequence----------------------------------------------------------------
## Change the length of the sequence according to convenience
string_vector<- function(length=10000){
  # The nucleotides in a single vector 
  n_vector<-c("A","T","G","C")
  # Sample function with the specified length
  simulation <-sample(n_vector,size = length,replace=T)
  return(paste(simulation, collapse = ""))
}
# Seed for reproducibility 
set.seed(42)
# Testing the function
random_sequence<-string_vector(10000)
# Transforming the random sequence to create a table with the function strsplit
dna_vector <- strsplit(random_sequence, "")[[1]]
# Table to see the number of nucleotides 
table<-table(dna_vector)
print(table)
print(dna_vector)
