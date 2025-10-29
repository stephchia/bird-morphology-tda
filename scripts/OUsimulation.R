# simulation of OU

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
# Step 3: Simulate trait evolution under Brownian motion
# --------------------------------------------------------------
# for each PC, simulate 10 sets of trait evolutionary data using the consensus tree
sim_list <- lapply(1:ncol(trait), function(pc) {
  data <- setNames(trait[, pc], rownames(trait))
  model <- phylolm(data ~ 1, phy = tree, model = "OUfixedRoot")
  set.seed(pc)
  # simulate traits (output: tip & node values)
  fastBM(tree, a = model$coef, sig2 = model$sigma2, alpha = model$optpar, theta = model$coef, nsim = 10, internal = TRUE)
})

# merge 10 sets of trait data into an object
anc_sim <- lapply(1:10, function(i) {
  do.call(cbind, lapply(1:4, function(pc) sim_list[[pc]][, i]))
})

saveRDS(anc_sim, "data/processed/anc_sim_ou.rds")

# ---------------------------------------------------------------------
# Step 4: Interpolate trait values through time (empirical + simulated)
# ---------------------------------------------------------------------
tmax <- 10  # set maximum Mya to be processed

## Simulated data
# anc_sim <- readRDS("data/processed/anc_sim.rds")
dir.create("data/simulate_ou", showWarnings = FALSE)

for (i in 1:10) {
  sim_dir <- file.path("data/simulate_ou", paste0("sim", i))
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
dataset <- anc_sim[[2]]
names <- c(tree$tip.label, (nrow(trait) + 1):nrow(dataset))
traits_sim <- lapply(1:4, function(i) setNames(dataset[, i], names))

## Plot 1D evolutionary trajectories for each PC
for (i in 1:4) {
  plot_1D_trajectory(tree, traits_sim[[i]], xbreaks = seq(0, 40, 5))
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
