#rarefy CCGP samples to minimum number of samples
#Jan. 22nd 2024
#jnadams

library(vegan)
data(BCI)
head(BCI)
dim(BCI)
# Your data should be organised like the BCI data - species in columns, sites in rows.

# Use rowSums to find the minimum number found at a site:
raremax <- min(rowSums(BCI))
raremax <- max(rowSums(BCI))

# Rarefy to that number
rarefied_data <- rrarefy(BCI, raremax)
max(rowSums(rarefied_data))
dim(rarefied_data)

# Loading library  
library(dplyr) 

# Create a data frame 
d <- data.frame(sample_number = c(00, 01, 02, 03))

# Printing three rows  
sample_n(d, 3) 
