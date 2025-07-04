# 05_persistent_homology.R
# -----------------------------
# This script computes persistent homology (H0 and H1) across time slices
# for empirical and simulated datasets using the TDA package.
# Intended to be run on an HPC cluster (e.g., SLURM).
# Recommend: assign ≥500GB RAM for mya = 0

library(TDA)
library(parallel)

# ----------------------------------------
# Function to compute persistence diagram
# ----------------------------------------
compute_rips <- function(data, max_dim = 1, max_scale = 10) {
  ripsDiag(
    X = data,
    maxdimension = max_dim,
    maxscale = max_scale,
    dist = "euclidean",
    library = "Dionysus",
    location = TRUE,
    printProgress = TRUE
  )
}

# ----------------------------------------
# Empirical dataset
# ----------------------------------------
for (mya in 20:0) {
  data <- read.csv(paste0("data/empirical/trait_", mya, "mya.csv"))
  tic <- proc.time()
  diag <- compute_rips(data)
  saveRDS(diag, file = paste0("data/empirical/tda_", mya, "mya.rds"))
  
  elapsed <- (proc.time() - tic)[3]
  cat(paste0("Empirical | ", mya, " Mya | Time: ", elapsed, " sec\n"))
}

# ----------------------------------------
# Simulated datasets
# ----------------------------------------
compute_sim_tda <- function(sim_index) {
  for (mya in 20:0) {
    data <- read.csv(paste0("data/simulate/sim", sim_index, "/trait_", mya, "mya.csv"))
    tic <- proc.time()
    diag <- compute_rips(data)
    saveRDS(diag, file = paste0("data/simulate/sim", sim_index, "/tda_", mya, "mya.rds"))
    
    elapsed <- (proc.time() - tic)[3]
    cat(paste0("Sim", sim_index, " | ", mya, " Mya | Time: ", elapsed, " sec\n"))
  }
}
mclapply(1:10, compute_sim_tda, mc.cores = 10)
