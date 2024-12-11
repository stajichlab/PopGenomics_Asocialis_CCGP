#!/usr/bin/env Rscript
#algatr test 
#Dec. 1st 2023 
#jnadams

library(algatr)
library(here)
library(raster)
library(sp)
library(vcfR)
library(adegenet)
library(tess3r)

pdf('algatr_plots.pdf')

#load example files
load_algatr_example()

#set working directory
getwd()
#read in vcf file
CCGP_vcf <- read.vcfR("vcf/ASCCGP_v3.All.SNP.combined_selected.vcf", verbose = FALSE)
head(CCGP_vcf, n = 20)
#read in filtered vcf file
CCGP_vcf_filtered <- read.vcfR("vcf/filtered.vcf", verbose = FALSE)
head(CCGP_vcf_filtered, n = 20)
tail(CCGP_vcf_filtered)
View(CCGP_vcf_filtered)
#CCGP_vcf_filtered_v2 <- is.na(CCGP_vcf_filtered)
#head(CCGP_vcf_filtered_v2)

#convert vcf to genotype matrix
Acarospora_dosage <- vcf_to_dosage(CCGP_vcf)
head(Acarospora_dosage)
#read in environmental data; we can use an environmental layer (PC1)
krig_raster <- raster::aggregate(CA_env[[1]], fact = 6)
terra::plot(CA_env[[1]], col = mako(100), axes = FALSE)

#run with compressed file
#read in vcf file
CCGP_vcf_gz <- read.vcfR("vcf/ASCCGP_v3.TestFresh.SNP.combined_selected.vcf.gz", verbose = FALSE)
head(CCGP_vcf_gz)
#run with compressed file from a year ago
CCGP_vcf_gz_v2 <- read.vcfR("vcf/ASCCGP_v1.All.INDEL.combined_selected.vcf.gz", verbose = FALSE)
head(CCGP_vcf_gz_v2)
#run with compressed file from a year ago with goodmapper 
CCGP_vcf_gz_v3 <- read.vcfR("vcf/ASCCGP_v1.GoodMapper.SNP.combined_selected.vcf.gz", verbose = FALSE)
head(CCGP_vcf_gz_v3)
CCGP_vcf_gz_v4 <- read.vcfR("vcf/ASCCGP_v1.All.SNP.combined_selected.vcf.gz", verbose = FALSE)
head(CCGP_vcf_gz_v4)
tail(CCGP_vcf_gz_v4)

CCGP_vcf_gz_v5 <- read.vcfR("vcf/ASCCGP_v3.All.SNP.selected_pruned.vcf.gz", verbose = FALSE)
head(CCGP_vcf_gz_v5)
tail(CCGP_vcf_gz_v5)

CCGP_vcf_gz_v6 <- read.vcfR("vcf/CCGP_output_filtered_pruned_v19.vcf.gz", verbose = FALSE)
head(CCGP_vcf_gz_v6)
tail(CCGP_vcf_gz_v6)

CCGP_vcf_gz_v7 <- read.vcfR("vcf/plink_v2.vcf.gz", verbose = FALSE)
head(CCGP_vcf_gz_v7)
tail(CCGP_vcf_gz_v7)

#convert vcf to genotype matrix
Acarospora_dosage_gz <- vcf_to_dosage(CCGP_vcf_gz)
head(Acarospora_dosage_gz)
Acarospora_dosage_gz <- vcf_to_dosage("vcf/ASCCGP_v3.TestFresh.SNP.combined_selected.vcf.gz")

#read in coordinates file
CCGP_coords <- read.csv("/bigdata/stajichlab/shared/projects/Population_Genomics/Asocialis_CCGP/CCGP_coords_v3.csv")
head(CCGP_coords)
tail(CCGP_coords)
#read in coordinates file for 11 samples
CCGP_coords_filtered <- read.csv("/bigdata/stajichlab/shared/projects/Population_Genomics/Asocialis_CCGP/CCGP_coords_filtered.csv")
head(CCGP_coords_filtered, n = 20)
tail(CCGP_coords_filtered)
#read in coordinate file that is reordered 
CCGP_coords_reordered <- read.csv("/bigdata/stajichlab/shared/projects/Population_Genomics/Asocialis_CCGP/CCGP_coords_reordered.csv")
head(CCGP_coords_reordered)
tail(CCGP_coords_reordered)
CCGP_coords_reordered_v2 <- read.csv("/bigdata/stajichlab/shared/projects/Population_Genomics/Asocialis_CCGP/CCGP_coords_reordered_v2.csv")
head(CCGP_coords_reordered_v2)
tail(CCGP_coords_reordered_v2)
CCGP_coords_reordered_v3 <- read.csv("/bigdata/stajichlab/shared/projects/Population_Genomics/Asocialis_CCGP/CCGP_coords_reordered_v3.csv")
head(CCGP_coords_reordered_v3)
tail(CCGP_coords_reordered_v3)
CCGP_coords_reordered_v4 <- read.csv("/bigdata/stajichlab/shared/projects/Population_Genomics/Asocialis_CCGP/CCGP_coords_reordered_v4.csv")
head(CCGP_coords_reordered_v4)
tail(CCGP_coords_reordered_v4)
#read in coordinates file
CCGP_coords_v2 <- read.csv("/bigdata/stajichlab/shared/projects/Population_Genomics/Asocialis_CCGP/CCGP_coords_v4.csv")
head(CCGP_coords_v2)
tail(CCGP_coords_v2)
#get rid of NA in CCGP_coords_v2
CCGP_coords_v3 <- na.omit(CCGP_coords_v2)
head(CCGP_coords_v3)
tail(CCGP_coords_v3)

tot <- rowSums(CCGP_coords_filtered)
print(tot)
CCGP_coords_filtered_v2 <- is.na(CCGP_coords_filtered)
head(CCGP_coords_filtered_v2)
tail(CCGP_coords_filtered_v2)

qmat <- qmatrix(tess3_obj, K = bestK)
# First, create a grid for kriging
# We can use one environmental layer (PC1), aggregated (i.e., increased cell size) to increase computational speed
krig_raster <- raster::aggregate(CA_env[[1]], fact = 6)
head(krig_raster)
krig_admix <- tess_krig(tess_CCGP_results_v6$qmat, tess_CCGP_results_v6$coords, tess_CCGP_results_v6$grid)

#Run tess do everything 
tess_CCGP_results_v2 <- tess_do_everything(CCGP_vcf_gz_v3, CCGP_coords_v3, grid = CA_env[[1]], Kvals = 1:10, K_selection = "auto")
tess_CCGP_results_v2 <- tess_do_everything(CCGP_vcf_gz_v2, CCGP_coords_v3, grid = CA_env[[1]], Kvals = 1:10, K_selection = "auto")
tess_CCGP_results_v3 <- tess_do_everything(CCGP_vcf_gz_v2, CCGP_coords_reordered, grid = CA_env[[1]], Kvals = 1:10, K_selection = "auto")
tess_CCGP_results_v4 <- tess_do_everything(CCGP_vcf_gz_v2, CCGP_coords_reordered_v2, grid = CA_env[[1]], Kvals = 1:10, K_selection = "auto")

tess_CCGP_results_v5 <- tess_do_everything(CCGP_vcf_gz_v2, CCGP_coords_reordered_v3, grid = CA_env[[1]], Kvals = 1:10, K_selection = "auto")
tess_CCGP_results_v6 <- tess_do_everything(CCGP_vcf_gz_v4, CCGP_coords_reordered_v3, grid = CA_env[[1]], Kvals = 1:10, K_selection = "auto")
tess_CCGP_results_v7 <- tess_do_everything(CCGP_vcf_gz_v5, CCGP_coords_reordered_v3, grid = CA_env[[1]], Kvals = 1:10, K_selection = "auto")
tess_CCGP_results_v8 <- tess_do_everything(CCGP_vcf_gz_v6, CCGP_coords_reordered_v3, grid = CA_env[[1]], Kvals = 1:10, K_selection = "auto")
tess_CCGP_results_v9 <- tess_do_everything(CCGP_vcf_gz_v7, CCGP_coords_reordered_v4, grid = CA_env[[1]], Kvals = 1:10, K_selection = "auto")

tess_CCGP_results_filtered <- tess_do_everything(CCGP_vcf_filtered, CCGP_coords_filtered, grid = CA_env[[1]], Kvals = 1:10, K_selection = "auto")
tess_CCGP_results_all <- tess_do_everything(CCGP_vcf, CCGP_coords_v3, grid = CA_env[[1]], Kvals = 1:10, K_selection = "auto")

#tess_CCGP_results_filtered <- tess_do_everything(CCGP_vcf_filtered_v2, CCGP_coords_filtered, grid = CA_env[[1]], Kvals = 1:10, K_selection = "auto")

#tess_CCGP_results_filtered_v2 <- tess_do_everything(CCGP_vcf_filtered, CCGP_coords_filtered_v2, grid = CA_env[[1]], Kvals = 1:10, K_selection = "auto")
class(Coord) #check the class of Coord
any(is.na(Coord))  # Check for missing values
all(is.numeric(Coord))  # Check if all values are numeric



par(mfrow = c(2, 2), pty = "s", mar = rep(0, 4))

tess_ggplot(tess_CCGP_results_v6$krig_admix, plot_method = "maxQ")



tess_CCGP_results_v6$krig_admix@cpp

#Run tess with k=3
#convert vcf to genotype matrix
Acarospora_dosage <- vcf_to_dosage(CCGP_vcf_gz_v7)
head(Acarospora_dosage)
tess3_obj_noK <- tess3(Acarospora_dosage, coord = as.matrix(CCGP_coords_reordered_v4), K = 3, method = "projected.ls", ploidy = 2)
# Get TESS object and best K from results
tess3_obj <- tess3_result$tess3_obj_noK
bestK <- tess3_result[["K"]]
summary(tess3_obj_noK)

# Get Qmatrix with ancestry coefficients
qmat <- qmatrix(tess3_obj_noK, K = 3)

# qmat contains ancestry coefficient values for each individual (row) and each K value (column)
head(qmat)

krig_admix <- tess_krig(qmat, CCGP_coords_reordered_v4, krig_raster)

#visualize tess results
tess_barplot(qmat)

tess_Acarospora <- tess_ggbarplot(qmat, legend = FALSE)
tess_Acarospora
head(tess_Acarospora$data)
head(tess_Acarospora$labels)

tess_barplot(qmat, legend = FALSE)

par(mfrow = c(2, 2), pty = "s", mar = rep(0, 4))

tess_ggplot(krig_admix, plot_method = "maxQ")

#run tess3r with labeled structure plot
# retrieve tess3 Q matrix for K = 3 clusters 
q.matrix <- qmatrix(tess3_obj_noK, K = 3)
# STRUCTURE-like barplot for the Q-matrix 
barplot(q.matrix, border = NA, space = 0, 
        xlab = "SampleID", ylab = "Ancestry proportions", 
        main = "Ancestry matrix") -> bp
axis(1, at = 1:nrow(q.matrix), labels = bp$order, las = 3, cex.axis = .4)

#run wingen do everything: 
CCGP_lyr <- coords_to_raster(CCGP_coords_reordered_v4, res = 0.5, buffer = 5)
envlayer <- rast(CA_env$CA_rPCA1)
CCGP_wingen_results <- wingen_do_everything(
  gen = CCGP_vcf_gz_v7,
  lyr = CCGP_lyr,
  coords = CCGP_coords_reordered_v4,
  wdim = 3,
  fact = 0,
  sample_count = TRUE,
  preview = FALSE,
  min_n = 2,
  stat = "pi",
  rarify = FALSE,
  kriged = TRUE,
  grd = CCGP_lyr,
  index = 1:2,
  agg_grd = NULL, disagg_grd = 4,
  agg_r = NULL, disagg_r = NULL,
  masked = TRUE, mask = envlayer,
  bkg = envlayer, plot_count = TRUE
)

#run mmrr do everything
#to run mmrr we need a genetic distance matrix, which can be generated from a vcf file
#Now let's calculate PC distances
pc_dists_Acarospora <- gen_dist(CCGP_vcf_gz_v7, dist_type = "pc", npc_selection = "auto", criticalpoint = 2.0234)
set.seed(01)
mmrr_full_Acarospora_everything <- mmrr_do_everything(pc_dists_Acarospora, CCGP_coords_reordered_v4, env = CA_env, geo = TRUE, model = "full")
#run the mmrr vignette
# Install packages
mmrr_packages()
#load required libs
library(algatr)
library(here)
library(raster)
#load example data
#load_algatr_example()
#make a matrix of gendist
Y_v2 <- as.matrix(pc_dists_Acarospora)
# Extract enviro vars
env_v3 <- raster::extract(CA_env, CCGP_coords_reordered_v4)
# Calculate environmental distances
X_v2 <- env_dist(env_v3)
# Add geographic distance to X
X_v2[["geodist"]] <- geo_dist(CCGP_coords_reordered_v4)
#Run mmrr_run
set.seed(10)
results_full_Acarospora <- mmrr_run(Y_v2, X_v2, nperm = 99, stdz = TRUE, model = "full")
# Run MMRR with all variables
set.seed(01)
results_best_Acarospora <- mmrr_run(Y_v2, X_v2, nperm = 99, stdz = TRUE, model = "best")
# Single variable plot
mmrr_plot(Y_v2, X_v2, mod = results_full_Acarospora$mod, plot_type = "vars", stdz = TRUE)
# Fitted variable plot
mmrr_plot(Y_v2, X_v2, mod = results_full_Acarospora$mod, plot_type = "fitted", stdz = TRUE)
# Covariance plot
mmrr_plot(Y_v2, X_v2, mod = results_full_Acarospora$mod, plot_type = "cov", stdz = TRUE)
#mmrr_plot best model
mmrr_plot(Y_v2, X_v2, mod = results_best_Acarospora$mod, plot_type = "all", stdz = TRUE)
#mmrr summary statistics table
mmrr_table(results_full_Acarospora, digits = 2, summary_stats = TRUE)
# here we aggregate the layer for computational speed
lyr_v2 <- aggregate(CA_env$CA_rPCA1, 50)
plot(lyr_v2)
points(CCGP_coords_reordered_v4)
# Recreate MMRR input with resistance distances
# Calculate environmental distances
X_v2 <- env_dist(env_v3)
# Add geographic distance to X
X_v2[["resistdist"]] <- geo_dist(CCGP_coords_reordered_v4, type = "resistance", lyr = lyr_v2)
# Run MMRR with resistance distances
results_resist_Acarospora <- mmrr_run(Y_v2, X_v2, nperm = 99, stdz = TRUE, model = "full")
mmrr_plot(Y_v2, X_v2, mod = results_resist_Acarospora$mod, plot_type = "all", stdz = TRUE)
#table of stats for results_resist
mmrr_table(results_resist_Acarospora)

#Run GDM do everything
gdm_full_Acarospora_everything <- gdm_do_everything(pc_dists_Acarospora,
                                         CCGP_coords_reordered_v4,
                                         envlayers = CA_env,
                                         model = "full",
                                         scale_gendist = TRUE
)

#run GDM vignette
# Install packages
gdm_packages()
library(algatr)
library(raster)
library(terra)
#load example datasets
load_algatr_example()
env_Acarospora <- raster::extract(CA_env, CCGP_coords_reordered_v4)
gdm_full_Acarospora <- gdm_run(
  gendist = pc_dists_Acarospora,
  coords = CCGP_coords_reordered_v4,
  env = env_Acarospora,
  model = "full",
  scale_gendist = TRUE
)

gdm_full_Acarospora$varimp
gdm_full_Acarospora$model
gdm_full_Acarospora$pvalues
CA_env$CA_rPCA1
CA_env$CA_rPCA2
CA_env$CA_rPCA3

#model selection with gdm_run; this code is not run within the vignette, the code does not work
#gdm_best <- gdm_run(gendist = liz_gendist, 
#                    coords = liz_coords, 
#                    env = env_liz, 
#                    model = "best", 
#                    scale_gendist = TRUE,
#                    nperm = 500, 
#                    sig = 0.05)

# Look at p-values
#gdm_best$pvalues
#gdm_best$varimp

summary(gdm_full_Acarospora$model)
#plot the gdm
gdm_plot_diss(gdm_full_Acarospora$model)
#plotting fitted I-splines for variables
par(mfrow = c(2, 2))
gdm_plot_isplines(gdm_full_Acarospora$model)
#generate a table of GDM statistics
gdm_table(gdm_full_Acarospora)
# Extract the GDM map from the GDM model object
map_Acarospora <- gdm_map(gdm_full_Acarospora$model, CA_env, CCGP_coords_reordered_v4, plot_vars = FALSE)

maprgb_Acarospora <- map_Acarospora$pcaRastRGB

# Now, use `extrap_mask()` to do buffer-based masking
map_mask_Acarospora <- extrap_mask(CCGP_coords_reordered_v4, maprgb_Acarospora, method = "buffer", buffer_width = 1.25)

# Plot the resulting masked map
p <- plot_extrap_mask(maprgb_Acarospora, map_mask_Acarospora, RGB = TRUE, mask_col = rgb(1, 1, 1, alpha = 0.6))
ggsave("extrap_mask.pdf",p,width=12,height=12)
#RDA vignette
# Install packages
rda_packages()

