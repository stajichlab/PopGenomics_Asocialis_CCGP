#run mmrr tutorial Jan. 3rd 2024
# Install packages
mmrr_packages()
install.packages("tseries")
library(tseries)
#run mmrr do everything
#to run mmrr we need a genetic distance matrix, which can be generated from a vcf file
#Now let's calculate PC distances
pc_dists_Acarospora_full <- gen_dist(CCGP_vcf_pruned, dist_type = "pc", npc_selection = "auto", criticalpoint = 2.0234)
set.seed(01)
#make a matrix of gendist
Y_full <- as.matrix(pc_dists_Acarospora_full)
head(Y_full)
dim(Y_full)

# Extract enviro vars
env_full <- raster::extract(CA_env, CCGP_coords_reordered_v4)
# Calculate environmental distances
X_full <- env_dist(env_full)
head(X_full)
dim(X_full)
# Add geographic distance to X
Z_full <- geo_dist(CCGP_coords_reordered_v4)
head(Z_full)
dim(Z_full)

# MMRR performs Multiple Matrix Regression with Randomization analysis
# Y is a dependent distance matrix
# X is a list of independent distance matrices (with optional names)

MMRR<-function(Y,X,nperm=999){
  #compute regression coefficients and test statistics
  nrowsY<-nrow(Y)
  y<-unfold(Y)
  if(is.null(names(X)))names(X)<-paste("X",1:length(X),sep="")
  Xmats<-sapply(X,unfold)
  fit<-lm(y~Xmats)
  coeffs<-fit$coefficients
  summ<-summary(fit)
  r.squared<-summ$r.squared
  tstat<-summ$coefficients[,"t value"]
  Fstat<-summ$fstatistic[1]
  tprob<-rep(1,length(tstat))
  Fprob<-1
  
  #perform permutations
  for(i in 1:nperm){
    rand<-sample(1:nrowsY)
    Yperm<-Y[rand,rand]
    yperm<-unfold(Yperm)
    fit<-lm(yperm~Xmats)
    summ<-summary(fit)
    Fprob<-Fprob+as.numeric(summ$fstatistic[1]>=Fstat)
    tprob<-tprob+as.numeric(abs(summ$coefficients[,"t value"])>=abs(tstat))
  }
  
  #return values
  tp<-tprob/(nperm+1)
  Fp<-Fprob/(nperm+1)
  names(r.squared)<-"r.squared"
  names(coeffs)<-c("Intercept",names(X))
  names(tstat)<-paste(c("Intercept",names(X)),"(t)",sep="")
  names(tp)<-paste(c("Intercept",names(X)),"(p)",sep="")
  names(Fstat)<-"F-statistic"
  names(Fp)<-"F p-value"
  return(list(r.squared=r.squared,
              coefficients=coeffs,
              tstatistic=tstat,
              tpvalue=tp,
              Fstatistic=Fstat,
              Fpvalue=Fp))
}

# unfold converts the lower diagonal elements of a matrix into a vector
# unfold is called by MMRR

unfold<-function(X){
  x<-vector()
  for(i in 2:nrow(X)) x<-c(x,X[i,1:i-1])
  x<-scale(x, center=TRUE, scale=TRUE)  # Comment this line out if you wish to perform the analysis without standardizing the distance matrices! 
  return(x)
}


# Tutorial for data files gendist.txt, geodist.txt, and ecodist.txt

# Read the matrices from files.
# The read.matrix function requires {tseries} package to be installed and loaded.
# If the files have a row as a header (e.g. column names), then specify 'header=TRUE', default is 'header=FALSE'.
#install tseries package

genMat <- Y_full
head(genMat)
geoMat <- Z_full
ecoMat_PCA1 <- X_full$CA_rPCA1
ecoMat_PCA2 <- X_full$CA_rPCA2
ecoMat_PCA3 <- X_full$CA_rPCA3
# Make a list of the explanatory (X) matrices.
# Names are optional.  Order doesn't matter.
# Can include more than two matrices, if desired.
Xmats <- list(geography=geoMat,ecology=ecoMat)
Xmats_v2 <- list(geography=geoMat,eco1=ecoMat_PCA1, eco2=ecoMat_PCA2, eco3=ecoMat_PCA3)
head(Xmats)
# Run MMRR function using genMat as the response variable and Xmats as the explanatory variables.
# nperm does not need to be specified, default is nperm=999)
MMRR(genMat,Xmats,nperm=999)
MMRR(genMat,Xmats_v2,nperm=999)
# These data should generate results of approximately:
# Coefficient of geography = 0.778 (p = 0.001)
# Coefficient of ecology = 0.167 (p = 0.063)
# Model r.squared = 0.727 (p = 0.001)
# Note that significance values may change slightly due to the permutation procedure.

