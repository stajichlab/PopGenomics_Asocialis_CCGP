#pophelper script
#Jan. 13th 2024
#jnadams

# install dependencies
#install.packages(c("ggplot2","gridExtra","label.switching","tidyr","remotes"),repos="https://cloud.r-project.org")
library(pacman)
p_load(ggplot2, gridExtra, label.switching, tidyr, 
       remotes, reshape2, dplyr, stringr, cowplot)
# install pophelper package from GitHub; one time install
remotes::install_github('royfrancis/pophelper')
# load library for use
library(pophelper)
# check version
packageDescription("pophelper", fields="Version")

#read in files
# use ADMIXTURE files (do not use this command to read local files)
afiles<- list.files(path="/bigdata/stajichlab/shared/projects/Population_Genomics/Asocialis_CCGP/vcf/admixture_runs", full.names=T)
alist <- readQ(files=afiles, indlabfromfile=T)
head(alist[[1]])

## tabulateQ
tr1 <- tabulateQ(qlist=alist)
head(tabulateQ(alist))

## summarizeQ
sr1 <- summariseQ(tr1, writetable=T, exportpath=getwd())
head(summariseQ(tabulateQ(alist)))

# choose plots to make
names(alist)

# pull group assignments for coloring pca
## merge K = 2
k2<-mergeQ(alignK(alist[c(10)]))

## merge K = 3
k3<-mergeQ(alignK(alist[c(12)]))

## merge K = 4
k4<-mergeQ(alignK(alist[c(14)]))

## merge K = 5
k5<-mergeQ(alignK(alist[c(16)]))

# align clusters b/t K
pop1 <- plotQ(alignK(c(k2, k3, k4, k5)),imgoutput="join", showindlab=T, returnplot=T,
            exportplot=F,basesize=11, sortind="all", sharedindlab = F)
grid.arrange(pop1$plot[[1]])
