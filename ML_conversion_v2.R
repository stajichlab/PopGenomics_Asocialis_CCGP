#conversion for ML tree
#Jan. 20th 2024
#jnadams

#set working directory
setwd("~/shared/projects/Population_Genomics/Asocialis_CCGP")
#load required libs
library(vcfR)
install.packages("poppr")
library(poppr)
library(ape)
library(RColorBrewer)
sale.VCF <- read.vcfR("vcf/CCGP_output_filtered_pruned_v20.vcf")
head(sale.VCF)

#VCF file to fasta for ML tree

gt <- extract.gt(sale.VCF, element = "GT")
gt[gt == "0/0"] <- 0
gt[gt == "0/1"] <- 1
gt[gt == "1/1"] <- 2
gt <- t(gt)
dim(gt)
gt.pyl <- apply(gt, 1, paste0, collapse="")
write.table(gt.pyl, sep = "\t", col.names = F, file = "142genomes_Jan20.phy")
