## ----setup, include=FALSE-------------------------------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE,warning = FALSE, message=FALSE,
                      fig.height = 4, fig.width = 8)


## -------------------------------------------------------------------------------------------------------------------
# Installing libraries to get the package ALL
#if (!require("BiocManager", quietly = TRUE))
    #install.packages("BiocManager")

#BiocManager::install("ALL")

# Loading libraries and data 
library("ALL")
# Acceding and processing data
data("ALL")
bcell = grep("^B", as.character(ALL$BT))
moltyp = which(as.character(ALL$mol.biol) 
    %in% c("NEG", "BCR/ABL"))
ALL_bcrneg = ALL[, intersect(bcell, moltyp)]
ALL_bcrneg$mol.biol = factor(ALL_bcrneg$mol.biol)
data<-exprs(ALL_bcrneg)
whp<-which(rownames(data)=="1636_g_at")
tempdat=data.frame(data[whp,],ALL_bcrneg$mol.biol)
colnames(tempdat)=c("y","mol_biol")

## Finding the means/group 
mean_bcr_abl <- mean(tempdat$y[tempdat$mol_biol == "BCR/ABL"])
mean_neg <- mean(tempdat$y[tempdat$mol_biol == "NEG"])

#Plotting
stripchart(y ~ mol_biol, data = tempdat,
           method = "jitter",
           jitter=0.2,
           vertical = TRUE,   
           pch = 1,                  
           xlab = "Molecular types",  
           ylab = "1636_g_at",          
           main = "Distribution of 1636_g at probe by cancer molecular subtypes",
            ) # Title of the plot

x_bcr_abl <- 1  # x-axis position for BCR/ABL
x_neg <- 2      

# Add blue lines for the means of each group
lines(c(x_bcr_abl - 0.1, x_bcr_abl + 0.1), c(mean_bcr_abl, mean_bcr_abl), col = 4, lwd = 2) 
lines(c(x_neg - 0.1, x_neg + 0.1), c(mean_neg, mean_neg), col = 4, lwd = 2)   



## ----Drawing and plotting uniform distributions---------------------------------------------------------------------
# Determining sample size
n <- 5
# Number of experiments 
n_experiments <- 200
# Empty vector to store the means 
sample_means <- numeric(n_experiments)
# Reproducibility
set.seed(123456)
# for loop function
for (i in 1:n_experiments) {
  # Give a uniform distribution for the sample size
  samples <- runif(n)
  # Obtaining the means 
  sample_means[i] <- mean(samples)
}
# Plotting the distribution using density function
plot(density(sample_means), main="Sample's mean distribution", xlab="Sample average", ylab="Frequency/Density")
rug(sample_means)
# Printing the mean and standard deviation of the 
Mean<-mean(sample_means)
STD <-sd(sample_means)
print(paste("The average of the samples is:",Mean))
print(paste("The standard deviation of the sample is:",STD))



## ----100 sample size------------------------------------------------------------------------------------------------
# New sample size
n <- 100
# creating a new empty vector 
sample_means <- numeric(n_experiments)
set.seed(123456)
# Foor loop 
for (i in 1:n_experiments) {
  samples <- runif(n)
  sample_means[i] <- mean(samples)
}
# Plotting distribution
plot(density(sample_means), main="100 Sample's mean distribution", xlab="Sample average", ylab="Frequency/Density")
rug(sample_means)
### Sample_means
Mean<-mean(sample_means)
STD <-sd(sample_means)
print(paste("The average of the samples is:",Mean))
print(paste("The standard deviation of the sample is:",STD))


## ----Simulating normal samples (a,b,c)------------------------------------------------------------------------------
# Seed for reproducibility
set.seed(123456)  
# creating Normal distribution sample
samples <- rnorm(30, mean = 0, sd = 1)
# Shuffling the indices using samples and then picking the treatment(second 15s) and control groups (first 15s)
group_indices <- sample(1:30, 30)  
## Control group 
control_group <- samples[group_indices[1:15]] 
## Treatment group
treatment_group <- samples[group_indices[16:30]]  
# Performing t-test to compare control vs treatment
t_test_1 <- t.test(control_group, treatment_group)
print(paste("The p-value is :",t_test_1$p.value))



## -------------------------------------------------------------------------------------------------------------------
# Reproducibility
set.seed(123456)
# Generating 1000 samples
un_samples <- runif(1000)
# Plotting the histogram 
hist(un_samples, main="Histogram of 1000 uniform samples distribution", xlab="Values", ylab="Distribution/Frequency", col="lightblue", border="black")


## -------------------------------------------------------------------------------------------------------------------
# Creating a empty vector with 1000 entries
p_values <- numeric(1000)
# Repeat the process 1000 times using for loop
for (i in 1:1000) {
  # (a) Simulate 30 samples from a normal distribution with mean 0 and sd 1
  samples <- rnorm(30, mean = 0, sd = 1)
  # (b) Randomly assign 15 samples to control and 15 to treatment group
  group_indices <- sample(1:30, 30)
  control_group <- samples[group_indices[1:15]]
  treatment_group <- samples[group_indices[16:30]]
  # (c) Perform two-sample t-test and store the p-value
  t_test_result <- t.test(control_group, treatment_group)
  p_values[i] <- t_test_result$p.value
}
# (d) Plot histogram of the 1000 p-values
hist(p_values, main="Histogram of p-values from 1000 sample T-tests", xlab="p-value", ylab="Distribution/Frequency", col="lightgreen", border="black")
lines(density(p_values), col="red", lwd=2)

