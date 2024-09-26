## ----setup, include=FALSE----------------------------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE,warning = FALSE, message=FALSE,
                      fig.height = 4, fig.width = 8)

## ----Function for contingency table k-tuples-----------------------------------------------------------------------
#Defining the function with the sequence and k-tuple input 
def<-function(sequence,k){
  #Sub stringing the sequence
  x<-substring(sequence,1:nchar(sequence),1:nchar(sequence))
  #creating empty vector to store the codons
  vector<-c()
  #Initializing function 
  for (i in seq(1,length(x)-k+1)){
    ans<-paste(x[i:(i+k-1)],collapse = "")
    vector<-c(vector,ans)
  }
  # Returning a table to count the number of codons 
  table_s<-table(vector)
  return(table_s)
}

# testing sequence called sequence
sequence <- "CAGACAAAACGAGAGAGAGAGGAAAGATAGATG"  
# Determining the length of k
k <- 3 
#Storing and printing the table 
codon_table <- def(sequence, k)
print(codon_table)