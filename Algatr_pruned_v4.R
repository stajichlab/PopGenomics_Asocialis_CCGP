#algatr analysis on 50kb pruned dataset without Rhizocarpon for a total of 142 samples
#rarify wingen analysis
#April 23rd 2024
#jnadams

#load required libs
library(algatr)
library(here)
library(raster)
library(sp)
library(vcfR)
library(adegenet)
library(tess3r)
library(wingen)
library(ggplot2)
library(terra)
library(fields)
library(rworldmap)
library(automap)
library(cowplot)

#set working directory
setwd("bigdata/stajichlab/shared/projects/Population_Genomics/Asocialis_CCGP/")
getwd()
#read in vcf file
CCGP_vcf_pruned <- read.vcfR("vcf/CCGP_output_filtered_pruned_v20.vcf", verbose = FALSE)
CCGP_vcf_str_dos <- read.vcfR(str_dos, verbose = FALSE)
head(CCGP_vcf_pruned, n = 20)
tail(CCGP_vcf_pruned, n = 20)
dim(CCGP_vcf_pruned)
#read in coordiantes file
CCGP_coords_reordered_v4 <- read.csv("/bigdata/stajichlab/shared/projects/Population_Genomics/Asocialis_CCGP/CCGP_coords_reordered_v4.csv")
head(CCGP_coords_reordered_v4)
tail(CCGP_coords_reordered_v4)

# Install packages for TESS
tess_packages()

# Convert vcf to genotype matrix
Acarospora_pruned_dosage <- vcf_to_dosage(CCGP_vcf_pruned)

# First, create a grid for kriging
# We can use one environmental layer (PC1), aggregated (i.e., increased cell size) to increase computational speed
krig_raster_Acarospora <- raster::aggregate(CA_env[[1]], fact = 6)

# If you want to see the difference between the non-aggregated (original) and aggregated rasters:
terra::plot(CA_env[[1]], col = mako(100), axes = FALSE)

#this code doesn't work and returns an error: plot rendering error
terra::plot(krig_raster_Acarospora, col = mako(100), axes = FALSE)

#run automatic k-selection with ploidy = 2
# Best K is 2; this provides a more reasonable estimate for the "best" K compared to manual selection above
tess3_result_Acarospora <- tess_ktest(Acarospora_pruned_dosage, CCGP_coords_reordered_v4, Kvals = 1:10, ploidy = 2, K_selection = "auto")

#run automatic k-selection with ploidy = 1; this does not work error in Check(X, ploidy)
# Best K is 2; this provides a more reasonable estimate for the "best" K compared to manual selection above
#tess3_result_Acarospora_v2 <- tess_ktest(Acarospora_pruned_dosage, CCGP_coords_reordered_v4, Kvals = 1:10, ploidy = 1, K_selection = "auto")

#run TESS with k = 2
tess3_obj_K2_Acarospora <- tess3(Acarospora_pruned_dosage, coord = as.matrix(CCGP_coords_reordered_v4), K = 2, method = "projected.ls", ploidy = 2)
# Get TESS object and best K from results
tess3_obj_K2_Acarospora <- tess3_result_Acarospora$tess3_obj
bestK_K2 <- tess3_result_Acarospora[["K"]]
# Get Qmatrix with ancestry coefficients
qmat_K2 <- qmatrix(tess3_obj_K2_Acarospora, K = bestK_K2)
# qmat contains ancestry coefficient values for each individual (row) and each K value (column)
head(qmat_K2)
#krige q-values
krig_admix_Acarospora_K2 <- tess_krig(qmat_K2, CCGP_coords_reordered_v4, krig_raster_Acarospora)
#visualize TESS results with K = 2
tess_barplot(qmat_K2)
tess_ggbarplot(qmat_K2, legend = FALSE)

#building maps in tess3 using tess_ggplot
par(mfrow = c(2, 2), pty = "s", mar = rep(0, 4))
tess_ggplot(krig_admix_Acarospora_K2, plot_method = "maxQ")

#run TESS with k = 3
tess3_obj_K3_Acarospora <- tess3(Acarospora_pruned_dosage, coord = as.matrix(CCGP_coords_reordered_v4), K = 3, method = "projected.ls", ploidy = 2)
#tess3_obj_noK <- tess3(liz_dosage, coord = as.matrix(liz_coords), K = 3, method = "projected.ls", ploidy = 2)
#run tess with k = 3
#tess3_obj_K4 <- tess3(Acarospora_dosage_v2, coord = as.matrix(CCGP_coords_reordered_v4), K = 4, method = "projected.ls", ploidy = 2)
# Get TESS object and best K from results
tess3_obj_K3_v2 <- tess3_obj_K3_Acarospora$tess3_obj
bestK_K3_v2 <- tess3_obj_K3_Acarospora[["K"]]
head(tess3_obj_K3_Acarospora)
tail(tess3_obj_K3_Acarospora)
# Get Qmatrix with ancestry coefficients
qmat_Acarospora_K3 <- qmatrix(tess3_obj_K3_Acarospora, K = bestK_K3_v2)
# qmat contains ancestry coefficient values for each individual (row) and each K value (column)
head(qmat_Acarospora_K3)
tail(qmat_Acarospora_K3)
#visualize results for k = 3
tess_barplot(qmat_Acarospora_K3)
tess_ggbarplot(qmat_Acarospora_K3, legend = FALSE)
#krige q-values
krig_admix_Acarospora_K3 <- tess_krig(qmat_Acarospora_K3, CCGP_coords_reordered_v4, krig_raster_Acarospora)
#building maps in tess3 using tess_ggplot
par(mfrow = c(2, 2), pty = "s", mar = rep(0, 4))
tess_ggplot(krig_admix_Acarospora_K3, plot_method = "maxQ")

#run TESS with k = 4
tess3_obj_K4_Acarospora <- tess3(Acarospora_pruned_dosage, coord = as.matrix(CCGP_coords_reordered_v4), K = 4, method = "projected.ls", ploidy = 2)
# Get TESS object and best K from results
tess3_obj_K4_v2 <- tess3_obj_K4_Acarospora$tess3_obj
bestK_K4_v2 <- tess3_obj_K4_Acarospora[["K"]]
head(tess3_obj_K4_Acarospora)
tail(tess3_obj_K4_Acarospora)
# Get Qmatrix with ancestry coefficients
qmat_Acarospora_K4 <- qmatrix(tess3_obj_K4_Acarospora, K = bestK_K4_v2)
# qmat contains ancestry coefficient values for each individual (row) and each K value (column)
head(qmat_Acarospora_K4)
tail(qmat_Acarospora_K4)
#visualize results for k = 4
tess_barplot(qmat_Acarospora_K4)
tess_ggbarplot(qmat_Acarospora_K4, legend = FALSE)
#sort by q-value
tess_ggbarplot(qmat_Acarospora_K4, sort_by_Q = TRUE, legend = FALSE)
#tess_ggbarplot(qmat_Acarospora_K4, sort_by_Q = FALSE, legend = FALSE)
#q-values for rows 70-80
qmat_Acarospora_K4[70:80,]
krig_admix_Acarospora_K4 <- tess_krig(qmat_Acarospora_K4, CCGP_coords_reordered_v4, krig_raster_Acarospora)
#building maps in tess3 using tess_ggplot
par(mfrow = c(2, 2), pty = "s", mar = rep(0, 4))
tess_ggplot(krig_admix_Acarospora_K4, plot_method = "maxQ")

#run wingen vignette
# Install packages
wingen_packages()
library(sf)

# First, we reformat our dataframe of coordinates into sf coordinates
coords_longlat <- st_as_sf(CCGP_coords_reordered_v4, coords = c("x", "y"), crs = "+proj=longlat")

# Next, the coordinates and raster can be projected to an equal area projection, in this case NAD83 / California Albers (EPSG 3310)
coords_proj <- st_transform(coords_longlat, crs = 3310)
#save one of the environmental PC layers as a SpatRaster object
envlayer <- rast(CA_env$CA_rPCA1)
envlayer <- rast(CA_env[[1]])

# Aggregate the layer so plotting is a bit faster
envlayer <- aggregate(envlayer, 5)

# Reproject to same crs as the projected coordinates
envlayer <- project(envlayer, crs(coords_proj))
#generate raster layer for sliding window
#Acarospora_lyr <- coords_to_raster(coords_proj, res = 10.0, buffer = 0)
Acarospora_lyr <- coords_to_raster(coords_proj, res = 50000, buffer = 5, plot = TRUE)
#preview the moving window
par(mfrow = c(1, 2), oma = rep(1, 4), mar = rep(2, 4))
preview_gd(Acarospora_lyr,
           CCGP_coords_reordered_v4,
           wdim = 3,
           fact = 0
)
#preview the moving window
sample_count <- preview_gd(Acarospora_lyr, coords_proj, wdim = 3, fact = 0)

# Visualize the sample count layer
ggplot_count(sample_count)
#Run the moving window
wgd_Acarospora <- window_gd(CCGP_vcf_pruned,
                            CCGP_coords_reordered_v4,
                            Acarospora_lyr,
                            stat = "pi",
                            wdim = 3,
                            fact = 0,
                            rarify = TRUE
)
wgd_Acarospora <- window_gd(CCGP_vcf_pruned,
                 coords_proj,
                 Acarospora_lyr,
                 stat = "pi",
                 wdim = 3,
                 fact = 0,
                 rarify = TRUE
)
#visualize wingen results using plot_gd() and plot_count()
par(mfrow = c(1, 2), oma = rep(1, 4), mar = rep(2, 4))
# Plot map of pi
plot_gd(wgd_Acarospora, main = "Moving window pi", bkg = envlayer)
# Plot map of pi
ggplot_gd(wgd_Acarospora, bkg = envlayer) + ggtitle("Moving window pi")
# Plot sample count map
plot_count(wgd_Acarospora, main = "Sample counts")

# Plot sample count map
ggplot_count(wgd_Acarospora) + ggtitle("Sample count")

#krige our moving window map using the krig_gd() function
kgd_Acarospora <- krig_gd(wgd_Acarospora, index = 1:2, Acarospora_lyr, disagg_grd = 5)
summary(kgd_Acarospora)

kgd_Acarospora <- krig_gd(wgd_Acarospora, index = 1:2, Acarospora_lyr, disagg_grd = 5)
summary(kgd_Acarospora)

#visualize kriged results
par(mfrow = c(1, 2), oma = rep(1, 4), mar = rep(2, 4))
# Plot kriged map of pi
plot_gd(kgd_Acarospora, main = "Kriged pi")
# Plot kriged map of pi
ggplot_gd(kgd_Acarospora) + ggtitle("Kriged pi")
# Plot kriged sample count map
plot_count(kgd_Acarospora, main = "Kriged sample counts")

# Plot kriged sample count map
ggplot_count(kgd_Acarospora) + ggtitle("Kriged sample counts")

#mask using sample counts
mgd_1 <- mask_gd(kgd_Acarospora, kgd_Acarospora[["sample_count"]], minval = 1)
mgd_2 <- mask_gd(kgd_Acarospora, kgd_Acarospora[["sample_count"]], minval = 2)

mgd_1 <- mask_gd(kgd_Acarospora, kgd_Acarospora[["sample_count"]], minval = 1)
mgd_2 <- mask_gd(kgd_Acarospora, kgd_Acarospora[["sample_count"]], minval = 2)

ggplot_gd(mgd_1, bkg = envlayer) + ggtitle("Kriged & masked pi")
#ggplot_gd(mgd_1, bkg = envlayer) + ggtitle("Kriged & masked pi")
ggplot_gd(mgd_2, bkg = envlayer) + ggtitle("Kriged & masked pi")

par(mfrow = c(1, 2), oma = rep(1, 4), mar = rep(2, 4))
plot_gd(mgd_1, main = "Kriged & masked pi", bkg = envlayer)
plot_gd(mgd_2, main = "Kriged & masked pi", bkg = envlayer)

#mask using state boundaries
# Resample envlayer based on masked layer
r <- terra::resample(envlayer, mgd_1)

# Perform masking
mgd_Acarospora <- mask_gd(mgd_1, r)

# Plot masked map
plot_gd(mgd_Acarospora, main = "Kriged & masked pi", bkg = r)

# Resample envlayer based on masked layer
r <- terra::resample(envlayer, mgd_1)

# Perform masking
mgd_Acarospora <- mask_gd(mgd_1, r)

# Plot masked map
ggplot_gd(mgd_Acarospora, bkg = envlayer) + ggtitle("Kriged & masked pi")

#run wingen_do_everything
results_Acarospora <- wingen_do_everything(
  gen = CCGP_vcf_pruned,
  lyr = Acarospora_lyr,
  coords = CCGP_coords_reordered_v4,
  wdim = 3,
  fact = 0,
  sample_count = TRUE,
  preview = FALSE,
  min_n = 2,
  stat = "pi",
  rarify = FALSE,
  kriged = TRUE,
  grd = Acarospora_lyr,
  index = 1:2,
  agg_grd = NULL, disagg_grd = 4,
  agg_r = NULL, disagg_r = NULL,
  masked = TRUE, mask = envlayer,
  bkg = envlayer, plot_count = TRUE
)

results_Acarospora <- wingen_do_everything(
  gen = CCGP_vcf_pruned,
  lyr = Acarospora_lyr,
  coords = coords_proj,
  wdim = 3,
  fact = 0,
  sample_count = TRUE,
  preview = FALSE,
  min_n = 2,
  stat = "pi",
  rarify = TRUE,
  kriged = TRUE,
  grd = Acarospora_lyr,
  index = 1:2,
  agg_grd = NULL, disagg_grd = 4,
  agg_r = NULL, disagg_r = NULL,
  masked = TRUE, mask = envlayer,
  bkg = envlayer, plot_count = TRUE
)

#run mmrr vigenette Dec. 18th 2023
# Install packages
mmrr_packages()
#run mmrr do everything
#to run mmrr we need a genetic distance matrix, which can be generated from a vcf file
#Now let's calculate PC distances
# Install packages 
gen_dist_packages()
pc_dists_Acarospora_full <- gen_dist(CCGP_vcf_pruned, dist_type = "pc", npc_selection = "auto", criticalpoint = 2.0234)
pc_dists_Acarospora_full_v2 <- gen_dist(str_dos_v2, dist_type = "pc", npc_selection = "auto", criticalpoint = 2.0234)
set.seed(01)
#make a matrix of gendist
Y_full <- as.matrix(pc_dists_Acarospora_full)

# Extract enviro vars
env_full <- raster::extract(CA_env, CCGP_coords_reordered_v4)
# Calculate environmental distances
X_full <- env_dist(env_full)
# Add geographic distance to X
X_full[["geodist"]] <- geo_dist(CCGP_coords_reordered_v4)

#run mmrr with mmrr_run() and all variables
set.seed(10)
results_full_Acarospora_v2 <- mmrr_run(Y_full, X_full, nperm = 99, stdz = TRUE, model = "full")

#run mmrr with variable selection
set.seed(01)
results_best_Acarospora_v2 <- mmrr_run(Y_full, X_full, nperm = 99, stdz = TRUE, model = "best")

#visualize mmrr results
# Single variable plot
mmrr_plot(Y_full, X_full, mod = results_full_Acarospora_v2$mod, plot_type = "vars", stdz = TRUE)
#run mmrr 
mmrr_plot_vars(Y_full, X_full, stdz = TRUE)
# Fitted variable plot
mmrr_plot(Y_full, X_full, mod = results_full_Acarospora_v2$mod, plot_type = "fitted", stdz = TRUE)
# Covariance plot
mmrr_plot(Y_full, X_full, mod = results_full_Acarospora_v2$mod, plot_type = "cov", stdz = TRUE)
#how do results compare to the best model?
mmrr_plot(Y_full, X_full, mod = results_best_Acarospora_v2$mod, plot_type = "all", stdz = TRUE)

#generate statistics with mmrr_table()
mmrr_table(results_full_Acarospora_v2, digits = 2, summary_stats = TRUE)

# here we aggregate the layer for computational speed
lyr_v2 <- aggregate(CA_env$CA_rPCA1, 50)
plot(lyr_v2)
points(CCGP_coords_reordered_v4)
# Recreate MMRR input with resistance distances
# Calculate environmental distances
X_full_v2 <- env_dist(env_full)
# Add geographic distance to X
X_full_v2[["resistdist"]] <- geo_dist(CCGP_coords_reordered_v4, type = "resistance", lyr = lyr_v2)
# Run MMRR with resistance distances
results_resist_Acarospora_v2 <- mmrr_run(Y_full, X_full_v2, nperm = 99, stdz = TRUE, model = "full")
mmrr_plot(Y_full, X_full_v2, mod = results_resist_Acarospora_v2$mod, plot_type = "all", stdz = TRUE)
#table of stats for results_resist
mmrr_table(results_resist_Acarospora_v2)

#run mmrr do everything
mmrr_full_Acarospora_everything_v2 <- mmrr_do_everything(pc_dists_Acarospora_full, CCGP_coords_reordered_v4, env = CA_env, geo = TRUE, model = "full")

#run gdm vignette
# Install packages
gdm_packages()
#extract environmental variables
env_Acarospora <- raster::extract(CA_env, CCGP_coords_reordered_v4)
#run gdm with gdm_run()
gdm_full_Acarospora_v2 <- gdm_run(
  gendist = pc_dists_Acarospora_full,
  coords = CCGP_coords_reordered_v4,
  env = env_Acarospora,
  model = "full",
  scale_gendist = TRUE
)
#interpreting gdm results
summary(gdm_full_Acarospora_v2$model)
#We can look at the relationships between both predicted ecological distance (the raw predictor) and predicted compositional dissimilarity against observed compositional dissimilarity using the gdm_plot_diss() function. 
gdm_plot_diss(gdm_full_Acarospora_v2$model)

#plotting fitted I-splines for variables
par(mfrow = c(2, 2))
gdm_plot_isplines(gdm_full_Acarospora_v2$model)

#table of gdm results
gdm_table(gdm_full_Acarospora_v2)

#visualizing gdm results
gdm_map(gdm_full_Acarospora_v2$model, CA_env, CCGP_coords_reordered_v4)

#masking out relevant areas
# Extract the GDM map from the GDM model object
map_Acarospora_v2 <- gdm_map(gdm_full_Acarospora_v2$model, CA_env, CCGP_coords_reordered_v4, plot_vars = FALSE)
maprgb_Acarospora_v2 <- map_Acarospora_v2$pcaRastRGB

# Now, use `extrap_mask()` to do buffer-based masking
map_mask_Acarospora_v2 <- extrap_mask(CCGP_coords_reordered_v4, maprgb_Acarospora_v2, method = "buffer", buffer_width = 1.25)

# Plot the resulting masked map
plot_extrap_mask(maprgb_Acarospora_v2, map_mask_Acarospora_v2, RGB = TRUE, mask_col = rgb(1, 1, 1, alpha = 0.6))

#run rda vignette
# Install packages
rda_packages()

#load required libs
library(algatr)
library(dplyr)
library(raster)
library(vegan)

load_algatr_example()
# Convert from vcf to dosage matrix
gen_dosage <- vcf_to_dosage(CCGP_vcf_pruned)
gen_dosage[1:100, 1:5]
# Are there NAs in the data?
gen_dosage[1:5, 1:5]
#impute missing values
str_dos <- str_impute(gen = gen_dosage, K = 1:5, entropy = TRUE, repetitions = 10, quiet = FALSE, save_output = FALSE)
str_dos_v2 <- str_impute(gen = CCGP_vcf_pruned, K = 1:5, entropy = TRUE, repetitions = 10, quiet = FALSE, save_output = FALSE)
# Check that NAs are gone
str_dos[1:5, 1:5]
str_dos[1:15, 1:15]
str_dos[1:100, 1:5]
head(str_dos)

# Extract environmental vars
env <- raster::extract(CA_env, CCGP_coords_reordered_v4)

# Standardize environmental variables and make into dataframe
env <- scale(env, center = TRUE, scale = TRUE)
env <- data.frame(env)
head(env)
tail(env)
View(env)

#replace NA with 0
env$CA_rPCA1[is.na(env$CA_rPCA1)] <- 0
env$CA_rPCA2[is.na(env$CA_rPCA2)] <- 0
env$CA_rPCA3[is.na(env$CA_rPCA3)] <- 0
is.na(env$CA_rPCA1)
is.na(env$CA_rPCA2)
is.na(env$CA_rPCA3)

#run a sipmle rda with no variable selection
mod_full <- rda_run(str_dos, env, model = "full")
mod_full$call
head(summary(mod_full))
RsquareAdj(mod_full)

#run a simple RDA with variable selection
mod_best <- rda_run(str_dos, env,
                    model = "best",
                    Pin = 0.05,
                    R2permutations = 1000,
                    R2scope = T
)


mod_best$call
mod_best$anova
RsquareAdj(mod_best)

#run a partial RDA with geography as a covariable
mod_pRDA_geo <- rda_run(str_dos, env, CCGP_coords_reordered_v4,
                        model = "full",
                        correctGEO = TRUE,
                        correctPC = FALSE
)

anova(mod_pRDA_geo)
RsquareAdj(mod_pRDA_geo) # 0.0305
head(summary(mod_pRDA_geo))

#run partial RDA with population genetic structure and geography as covariables
mod_pRDA_gs <- rda_run(str_dos, env, CCGP_coords_reordered_v4,
                       model = "full",
                       correctGEO = TRUE,
                       correctPC = TRUE,
                       nPC = 4
)

head(summary(mod_pRDA_gs))

#variance partitioning with partial RDA
varpart <- rda_varpart(str_dos, env, CCGP_coords_reordered_v4,
                       Pin = 0.05, R2permutations = 1000,
                       R2scope = T, nPC = 3
)
varpart
rda_varpart_table(varpart, call_col = TRUE)
#identifying candidate loci
mod_pRDA <- rda_run(str_dos, env, model = "best", correctPC = TRUE, nPC = 2)
mod_pRDA_full <- rda_run(str_dos, env, model = "full", correctPC = TRUE, nPC = 2)
View(mod_pRDA_full)
mod_pRDA_full$anova
rda_plot(mod_pRDA_full, axes = "all", binwidth = 20)

rda_sig_z <- rda_getoutliers(mod_pRDA_full, naxes = "all", outlier_method = "z", z = 3, plot = TRUE)

# How many outlier SNPs were detected?
length(rda_sig_z$rda_snps)

rda_sig_p <- rda_getoutliers(mod_best, naxes = "all", outlier_method = "p", p_adj = "fdr", sig = 0.01, plot = TRUE)

# How many outlier SNPs were detected?
length(rda_sig_p$rda_snps)

# Extract SNP names; choices is number of axes
snp_names <- rownames(scores(mod_best, choices = 2, display = "species"))

# Identify outliers that have q-values < 0.1
q_sig <-
  rda_sig_p$rdadapt %>%
  mutate(snp_names = snp_names) %>%
  filter(q.values <= 0.1)

# How many outlier SNPs were detected?
nrow(q_sig)

Reduce(intersect, list(
  q_sig$snp_names,
  rda_sig_p$rda_snps,
  rda_sig_z$rda_snps
))

rda_plot(mod_best, rda_sig_p$rda_snps, biplot_axes = c(1, 2), rdaplot = TRUE, manhattan = FALSE)

# In the case of our partial RDA, there was only one RDA axis, so a histogram is generated
rda_plot(mod_pRDA_full, rda_sig_z$rda_snps, rdaplot = TRUE, manhattan = FALSE, binwidth = 0.01)
#plot manhattan plot
rda_plot(mod_best, rda_sig_p$rda_snps, rda_sig_p$pvalues, rdaplot = FALSE, manhattan = TRUE)
#calculate significance threshold for RDA
signif_v2 <- 0.05 / length(rda_sig_p$rda_snps)
#plot manhattan plot with corrected threshold
rda_plot(mod_best, rda_sig_p$rda_snps, rda_sig_p$pvalues, sig = signif_v2, rdaplot = FALSE, manhattan = TRUE)

# Extract genotypes for outlier SNPs
rda_snps <- rda_sig_p$rda_snps
rda_gen <- str_dos[, rda_snps]

# Run correlation test
cor_df <- rda_cor(rda_gen, env)

# Make a table from these results (displaying only the first 5 rows):
rda_table(cor_df, nrow = 5)

#order by strength of the correlation

# Order by the strength of the correlation
rda_table(cor_df, order = TRUE, nrow = 5)

# Only retain the top variable for each SNP based on the strength of the correlation
rda_table(cor_df, top = TRUE, nrow = 5)

# Display results for only one environmental variable
rda_table(cor_df, var = "CA_rPCA2", nrow = 5)

#run lfmm vignette
# Install packages
lfmm_packages()

#load required libs
library(algatr)
library(here)
library(raster)

load_algatr_example()
# Convert from vcf to dosage matrix
gen_dosage <- vcf_to_dosage(CCGP_vcf_pruned)

# Are there NAs in the data?
gen_dosage[1:5, 1:5]
#impute missing values
str_dos <- str_impute(gen = gen_dosage, K = 1:5, entropy = TRUE, repetitions = 10, quiet = FALSE, save_output = FALSE)
# Check that NAs are gone
str_dos[1:5, 1:5]
env <- raster::extract(CA_env, CCGP_coords_reordered_v4)

#replace NA with 0
env$CA_rPCA1[is.na(env$CA_rPCA1)] <- 0
env$CA_rPCA2[is.na(env$CA_rPCA2)] <- 0
env$CA_rPCA3[is.na(env$CA_rPCA3)] <- 0
is.na(env$CA_rPCA1)
is.na(env$CA_rPCA2)
is.na(env$CA_rPCA3)

# Keep relevant params but retain default values for them
options(expressions = 5e5)
Cstack_info()
select_K(str_dos, K_selection = "tracy_widom", criticalpoint = 2.0234) # this gives me an error: protect(): protection stack overflow

select_K(str_dos, K_selection = "quick_elbow", low = 0.08, max.pc = 0.90) # 4

select_K(str_dos, K_selection = "tess", coords = CCGP_coords_reordered_v4, Kvals = 1:10) # 3

select_K(str_dos, K_selection = "find_clusters", perc.pca = 90, max.n.clust = 10) # 5

ridge_results <- lfmm_run(str_dos, env, K = 5, lfmm_method = "ridge")
lasso_results <- lfmm_run(str_dos, env, K = 5, lfmm_method = "lasso")


# Build tables for each of our LFMM runs, displaying only significant SNPs and ordering according to effect size (B)
lfmm_table(lasso_results$df, order = TRUE)
lfmm_table(ridge_results$df, order = TRUE)

lfmm_table(lasso_results$df, order = TRUE, var = "CA_rPCA3")

#plot a qqplot
lfmm_qqplot(lasso_results$df)

#plot a manhattan plot with the lasso method
lfmm_manhattanplot(lasso_results$df, sig = 0.05)
#plot a manhattan plot with corrected significance threshold

signif <- 0.05 / dim(CCGP_vcf_pruned)[1]

lfmm_manhattanplot(lasso_results$df, sig = signif)

#plot a manhattan plot with the ridge method
lfmm_manhattanplot(ridge_results$df, sig = 0.05)

#plot a manhattan plot with the ridge method with corrected signficant threshold 
lfmm_manhattanplot(ridge_results$df, sig = signif)