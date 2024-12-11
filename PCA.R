#PCA script using adegenet
#Jan. 14th 2024
#jnadams

#load required libs
install.packages("devtools")
library("devtools")
install_github("thibautjombart/adegenet")
library("adegenet")
install.packages("ape")
library("ape")
install.packages("pegas")
library("pegas")

#read in the data
#read in vcf file
CCGP_vcf_pruned <- read.vcfR("vcf/CCGP_output_filtered_pruned_v20.vcf.gz", verbose = FALSE)
head(CCGP_vcf_pruned, n = 20)
tail(CCGP_vcf_pruned, n = 20)

#convert vcfR file to genind object
my_genind <- vcfR2genind(CCGP_vcf_pruned)
head(my_genind)

x.lichen <- tab(my_genind, freq=TRUE, NA.method="mean")
pca.lichen <- dudi.pca(x.lichen, center=TRUE, scale=FALSE, scannf=FALSE, nf=10)
pca.lichen

s.label(pca.lichen$li)
s.class(pca.lichen$li, fac=pop(my_genind), col=funky(15))

s.class(pca.lichen$li, fac=pop(my_genind),
        col=transp(funky(15),.6),
        axesel=FALSE, cstar=0, cpoint=3)
add.scatter.eig(pca.lichen$eig[1:50],3,1,2, ratio=.3)
sum(is.na(my_genind$tab))

X <- tab(my_genind, freq = TRUE, NA.method = "mean")
class(X)
dim(X)
X[1:5,1:5]
pca1 <- dudi.pca(X, scale = FALSE, scannf = FALSE, nf = 3)
barplot(pca1$eig[1:50], main = "PCA eigenvalues", col = heat.colors(50))
pca1
pca1$li <- pca1$li %>% 
  as_tibble(rownames = "individuals") %>%
  mutate(individuals = str_extract(individuals, "(\\d+$)"), individuals = case_when(is.na(individuals) ~ "55", TRUE ~ individuals)) %>%
  column_to_rownames("individuals")

s.label(pca1$li, boxes = FALSE, clabel = 0.7, grid = TRUE, addaxes = TRUE)
title("PCA of lichen dataset\naxes 1-2")
add.scatter.eig(pca1$eig[1:20], 3,1,2)

s.class(pca1$li, pop(X))
title("PCA of lichen dataset\naxes 1-2")
add.scatter.eig(pca1$eig[1:20], 3,1,2)

#make a PCA based on RGB scale with 2 PCs
colorplot(pca1$li, pca1$li, transp=TRUE, cex=3, xlab="PC 1", ylab="PC 2")
title("PCA of lichen dataset\naxes 1-2")
abline(v=0,h=0,col="grey", lty=2)
par(new = TRUE)
s.label(pca1$li, boxes = FALSE, clabel = 0.7, grid = FALSE, addaxes = FALSE)
#make a PCA based on RGB scale with PC1 and PC3
colorplot(pca1$li[c(1,3)], pca1$li, transp=TRUE, cex=3, xlab="PC 1", ylab="PC 3")
title("PCA of lichen dataset\naxes 1-3")
abline(v=0,h=0,col="grey", lty=2)

