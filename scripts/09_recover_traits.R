# 09_recover_traits.R
# -----------------------------------
# Recover original trait values from a focal topological gap and examine species nearby

library(dplyr)

source("scripts/functions_gap_computation.R")

# ---------------------------------------------
# Extract focal gap geometry
# ---------------------------------------------
tda <- readRDS("data/empirical/tda_0mya.rds")
fg_coords <- get_gap_coordinates(tda)[[468]] # the focal gap coordinates
centroid_size <- compute_gap_centroid_size(fg_coords)
centroid <- centroid_size$centroid
size <- centroid_size$size

# ---------------------------------------------
# Recover original trait values from centroid
# ---------------------------------------------
pca <- readRDS("data/processed/pca.rds")
passeroid_scaled <- readRDS("data/processed/trait_passeroid.rds")
non_passeroid <- readRDS("data/processed/trait_non_passeroid_raw.rds")

# Inverse PCA and scaling to recover original (log-transformed) trait values
trait_org <- centroid %*% t(pca$rotation[, 1:4])
trait_org <- trait_org * attr(passeroid_scaled, "scaled:scale") + attr(passeroid_scaled, "scaled:center")

# Reverse log-transformation (except HWI)
trait_org[c(1:7, 9, 10)] <- exp(trait_org[c(1:7, 9, 10)])
trait_org 

# ---------------------------------------------
# Project non-passeroid species into PC space
# ---------------------------------------------
passeroid_pca <- read.csv("data/processed/trait_passeroid_pca.csv", row.name = 1)

# Apply the same scaling and projection to non-passeroids
non_passeroid_scaled <- scale(non_passeroid,
                              center = attr(passeroid_scaled, "scaled:center"),
                              scale = attr(passeroid_scaled, "scaled:scale"))
# non_passeroid_scaled <- t(apply(non_passeroid, 1, function(x) (x - attr(passeroid_scaled, "scaled:center")) / attr(passeroid_scaled, "scaled:scale")))
non_passeroid_pca <- as.data.frame(as.matrix(non_passeroid_scaled) %*% pca$rotation[, 1:4])

# --------------------------------------------------------
# Find species within or near the focal gap
# --------------------------------------------------------
# Non-passeroids within gap (distance < size/2)
dist_centroid_nonpas <- sapply(1:nrow(non_passeroid_pca), function(x) sqrt(sum((non_passeroid_pca[x, ] - centroid)^2)))
nonpas_within_gap <- non_passeroid_pca[which(dist_centroid_nonpas < size/2), ]
rownames(nonpas_within_gap)

# Passeroids near gap (distance < size)
dist_centroid_pas <- sapply(1:nrow(passeroid_pca), function(x) sqrt(sum((passeroid_pca[x, ] - centroid)^2)))
pas_near_gap <- passeroid_pca[which(dist_centroid_pas < size), ]
rownames(pas_near_gap)


