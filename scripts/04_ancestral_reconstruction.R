# 04_ancestral_reconstruction.R
# -----------------------------
# This script performs ancestral state reconstruction, simulates null trait evolution,
# interpolates traits over evolutionary time, and visualizes evolutionary trajectories.

library(dplyr)
library(phytools)
library(phylolm)
library(ape)
library(viridis)

# Load custom functions
source("scripts/functions_ancestral_reconstruction.R")

# --------------------------------------------------------------
# Step 1: Load trait data and consensus tree
# --------------------------------------------------------------
trait <- read.csv("data/processed/trait_passeroid_pca.csv", row.names = 1)
tree <- readRDS("data/processed/contree_pas.rds")

# --------------------------------------------------------------
# Step 2: Ancestral state reconstruction (empirical)
# --------------------------------------------------------------
# estimate node values for each PC (~2 minutes)
anc_emp <- lapply(1:4, function(i) {
  fastAnc(tree, setNames(trait[, i], rownames(trait)))
})

# put all tip and node values into one data frame
anc_emp <- do.call(cbind, anc_emp)
colnames(anc_emp) <- colnames(trait)
anc_emp <- rbind(trait[tree$tip.label, ], anc_emp)
rownames(anc_emp) <- 1:nrow(anc_emp)

saveRDS(anc_emp, "data/processed/anc_emp.rds")

# --------------------------------------------------------------
# Step 3: Simulate trait evolution under Brownian motion
# --------------------------------------------------------------
# for each PC, simulate 10 sets of trait evolutionary data using the consensus tree
sim_list <- lapply(1:ncol(trait), function(pc) {
  data <- setNames(trait[, pc], rownames(trait))
  model <- phylolm(data ~ 1, phy = tree, model = "BM")
  set.seed(pc)
  # simulate traits (output: tip & node values)
  fastBM(tree, a = model$coef, sig2 = model$sigma2, nsim = 10, internal = TRUE)
})

# merge 10 sets of trait data into an object
anc_sim <- lapply(1:10, function(i) {
  do.call(cbind, lapply(1:4, function(pc) sim_list[[pc]][, i]))
})

saveRDS(anc_sim, "data/processed/anc_sim.rds")

# ---------------------------------------------------------------------
# Step 4: Interpolate trait values through time (empirical + simulated)
# ---------------------------------------------------------------------
tmax <- 10  # set maximum Mya to be processed

## Empirical data
anc_emp <- readRDS("data/processed/anc_emp.rds")
dir.create("data/empirical", showWarnings = FALSE)

# for every million year, interpolate trait values and save as a separate file
for (mya in 1:tmax) {
  interpolated <- get_traits_at_timepoint(tree, anc_emp, mya)
  write.csv(interpolated, file = paste0("data/empirical/trait_", mya, "mya.csv"), row.names = FALSE)
}
write.csv(trait, "data/empirical/trait_0mya.csv", row.names = FALSE) # save modern time trait value

## Simulated data
# anc_sim <- readRDS("data/processed/anc_sim.rds")
dir.create("data/simulate", showWarnings = FALSE)

for (i in 1:10) {
  sim_dir <- file.path("data/simulate", paste0("sim", i))
  dir.create(sim_dir, showWarnings = FALSE)
  for (mya in 0:tmax) {
    interpolated <- get_traits_at_timepoint(tree, anc_sim[[i]], mya)
    write.csv(interpolated, file = file.path(sim_dir, paste0("trait_", mya, "mya.csv")), row.names = FALSE)
  }
}

# --------------------------------------------------------------
# Step 5: Visualize evolutionary trajectories
# --------------------------------------------------------------
# Build named vectors for traits (tip + internal nodes)
names <- c(tree$tip.label, (nrow(trait) + 1):nrow(anc_emp))
traits_emp <- lapply(1:4, function(i) setNames(anc_emp[, i], names))

## Plot 1D evolutionary trajectories for each PC
for (i in 1:4) {
  plot_1D_trajectory(tree, traits_emp[[i]], xbreaks = seq(0, 40, 5))
}

## Plot 2D trajectory with color gradient by time (PC1 vs PC2)
H <- nodeHeights(tree)
contmap <- contMap(tree, setNames(trait[, 1], rownames(trait)))
colors <- viridis(1500)[c(1:800, seq(801, 1500, length.out = 201))]
contmap$cols[] <- colors

# Name branch segments by relative time
for (i in 1:nrow(H)) {
  segment_times <- round((H[i, 1] + cumsum(contmap$tree$maps[[i]])) / max(H) * 1000)
  names(contmap$tree$maps[[i]]) <- segment_times
}

# Plot morphospace trajectory (PC1 vs. PC2)
phylomorphospace(contmap$tree,
                 trait[, 1:2],
                 anc_emp[(nrow(trait) + 1):nrow(anc_emp), 1:2],
                 colors = contmap$cols,
                 lwd = 1.5, node.size = c(0, 0), label = "off")
