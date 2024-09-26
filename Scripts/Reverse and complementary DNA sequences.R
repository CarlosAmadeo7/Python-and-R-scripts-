## ----setup, include=FALSE----------------------------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE,warning = FALSE, message=FALSE,
                      fig.height = 4, fig.width = 8)


## ----Proving that is a DNA sequence to proceed-------------------------------------------------------------------
# Testing sequence
x<-c("ATATATATATATAGAGGAGAGAGAT")
#################################

isdna<-function(sequence){
  # Sub stringing the sequence 
  sequence<-substring(sequence,nchar(sequence),nchar(sequence))
  # Testing if the sequence belongs to the vector of nucleotides
  isdna<-all(sequence%in%c("A","T","C","G"))
  # if and else 
  if (isdna){
    print("This is valid sequence, please proceed")
  } else{
    print("give me a valid sequence")
  }
}
# Testing the testing sequence for the function
isdna(x)
## ----Writing a function for complement and reverse sequence------------------------------------------------------
# Writing a function to get the reverse and complement sequence
RevComp<-function(sequence, option ="Complementary_Reverse"){
  # Defining the complementary sequence of the testing sequence in a vector 
  complementary_part = c(A="T",C="G",G="C",T="A")
  # Splitting the sequence using substring, using from 1 to the respective length
  splitted_Seq = substring(sequence, 1:nchar(sequence), 1:nchar(sequence))
  # Reversing the split sequence using the function rev 
  reverse_Seq <-rev(splitted_Seq)
  # Finding the complementary sequence going through the complementary part and the original sequence
  compl_Seq <-complementary_part[splitted_Seq]
  # Adding both the reverse and complementary 
  both <- complementary_part[reverse_Seq]
  
  # If statement for multiple options
  if (option == "Complementary"){
    return(paste(compl_Seq,collapse = ""))
  } else if (option =="Reverse"){
    return(paste( reverse_Seq, collapse = ""))
  } else if (option == "Complementary_Reverse"){
    return(paste(both ,collapse = "" ))
  } else{
    print("There is not option like that, please chack again ")
  }
}
# There are three options for displaying:
# 1. Complementary <- it will display the complementary sequence.
# 2. Reverse <- it will display the Reverse sequence.
# 3. Complementary_Reverse <- it will display the Complementary_Reverse sequence.

# Testing the function :
Complementary<-RevComp(x,option="Complementary")
Reverse<-RevComp(x,option="Reverse")
Compl_and_rev<- RevComp(x,option="Complementary_Reverse")
print(paste("This is the complementary sequence :",Complementary))
print(paste("This is the reverse sequence:",Reverse))
print(paste("This is the complementary and reverse sequence:",Compl_and_rev))
## ----Writing a function for multiple sequences and export function in a text file--------------------------------
# Defining the Mult_seq function that gives you the reverse and complementary sequence for different sequence
Mult_seq<-function(mult_seq, option="Complementary_Reverse"){
  #Applying sapply for all the sequences
  rev_mult_sequences<-sapply(mult_seq,RevComp,option=option)
  # Returning the rev_mul_sequences
  return(rev_mult_sequences)
}
# Testing the function with a vector of 5 small 4 bp sequences
mult_seq <- c("ATGG", "CGTA", "GAGA","TCAT","CAGA","TACG","AAAA","TAGC")
# Three options to specify depends of individual purposes.
results<-Mult_seq(mult_seq, option="Reverse")
print(results)

## ----Writing an export function----------------------------------------------------------------------------------
# Writing an export function to a text file
############### Setting up your path file where you want to save the text file 
file_path <- "put_your_file_in_here"
###############
#Creating the Exporting function 
Exporting<- function(sequences,file_path=file_path){
  write.table(sequences,file = file_path,col.names = F,row.names = F,quote = F)
  cat("The sequences were saved into", file_path, "\n")
}
# Testing the function, you need to give a name to the text file 
Exporting(results,file_path="sequences.txt")

