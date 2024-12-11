#convervsion for ML tree
#Jan. 19th 2024
#jnadams

#load required libs
library(vcfR)
library(poppr)
library(ape)
library(RColorBrewer)
sale.VCF <- read.vcfR(LD_88genomes_AUG29.vcf.gz)
head(sale.VCF)

#VCF file to fasta for ML tree

gt <- extract.gt(sale.VCF, element = GT)
gt[gt == 0/0] <- 0
gt[gt == 0/1] <- 1
gt[gt == 1/1] <- 2
gt <- t(gt)
dim(gt)
gt.pyl <- apply(gt, 1, paste0, collapse=“”)
write.table(gt.pyl, sep = “\t”, col.names = F, file = “88genomes_AUG29.phy”)