#script for calculating Fst with popgenome
#Jan 26th 2024
#jnadams

#load required libs
install.packages("devtools")
library(devtools)

devtools::install_github("pievos101/PopGenome")
library(PopGenome)

#read in fasta file
# Reading data, read fasta from folder 
GENOME.class <- readData("fasta") 
#read in meta data
meta2 <- read_csv("CCGP_Metadata_142isolates_v1.csv")
get.sum.data(GENOME.class)
get.individuals(GENOME.class)

# Available statistics and examples 
show.slots(GENOME.class) 
# Run necessary module 
GENOME.class <- F_ST.stats(GENOME.class)
GENOME.class <- F_ST.stats(GENOME.class, list(Socal_ind, Norcal_ind))
View(GENOME.class)
GENOME.class <- neutrality.stats(GENOME.class) 
GENOME.class@n.sites 
GENOME.class@Pi
GENOME.class@Tajima.D

# Available region data and statistics 
GENOME.class@region.data
GENOME.class@region.stats 
# Examples
GENOME.class@region.data@biallelic.sites[[1]][1:10]
GENOME.class@region.data@transitions[[1]][1:10]

# Without defining populations 
get.individuals(GENOME.class)
GENOME.class <- neutrality.stats(GENOME.class,detail=TRUE)
get.neutrality(GENOME.class)[[1]] 
# Define populations with lists
Socal_ind <- meta2 %>% filter(NorthSouth == "Southern CA") %>% dplyr::select(label)
Norcal_ind <- meta2 %>% filter(NorthSouth == "Northern CA") %>% dplyr::select(label)

GENOME.class <- set.populations(GENOME.class,list(c("Socal_ind", "Norcal_ind")))
View(GENOME.class)
# Check whether grouping is set correctly 
GENOME.class@region.data@populations
GENOME.class@region.data@populations2 
GENOME.class@region.data@outgroup
# Recalculate statistics for populations 
GENOME.class <-neutrality.stats(GENOME.class,detail=TRUE) 
GENOME.class@Tajima.D 
# Each population 
get.neutrality(GENOME.class)[[1]]
get.neutrality(GENOME.class)[[2]] 
# Set an outgroup 
GENOME.class <-set.outgroup(GENOME.class,c("Alyrata"))
GENOME.class@region.data@outgroup 
GENOME.class <- neutrality.stats(GENOME.class,detail=TRUE)
get.neutrality(GENOME.class)[[1]] 
get.neutrality(GENOME.class)[[2]]

#analyzing vcf files for whole genome sequence data
# What parameters need to be defined 
GENOME.class.vcf <- readData("vcf_popgenome")
GENOME.class.vcf@region.names 
GENOME.class.vcf <- neutrality.stats(GENOME.class.vcf, FAST=TRUE) 
get.sum.data(GENOME.class.vcf) 
GENOME.class.vcf@region.data
# FST 
GENOME.class.vcf <- F_ST.stats(GENOME.class.vcf)
get.F_ST(GENOME.class.vcf)[[1]] 
GENOME.class.vcf@nucleotide.F_ST

#loading vcf files with annotation
GENOME2.class <- readData("great_tit/vcf2",format="VCF", gffpath="great_tit/gff") 
get.sum.data(GENOME2.class)
GENOME2.class@region.data 
GENOME2.class <- set.synnonsyn(GENOME2.class, ref.chr="great_tit/fasta/LGE22.fasta")
GENOME2.class@region.data@synonymous
GENOME2.class@region.data@CodingSNPS 
GENOME2.class.syn <- neutrality.stats(GENOME2.class,subsites="syn")
GENOME2.class.syn@Tajima.D 
GENOME2.class.syn@theta_Watterson

#plot a site frequency spectrum for each population
# Concatenate loci
CON <- concatenate.regions(GENOME.class) 
CON <- detail.stats(CON,site.spectrum=TRUE,site.FST=TRUE) 
results <-get.detail(CON) 
allele_Freqs <- CON@region.stats@minor.allele.freqs[[1]] 
freq.table <- list()
freq.table[[1]] <- table(allele_Freqs) 
sfs <- data.frame(freq.table)

library(ggplot2) 
ggplot(sfs, aes(x=allele_Freqs, y=Freq)) + geom_bar(stat = ’identity’)
