# 03_build_consensus_tree.R
# --------------------------
# This script builds a consensus phylogenetic tree for Passeroidea 
# from 1,000 sampled trees from the Hackett backbone.

library(phytools)

# ------------------------------------------------------------------
# STEP 1: Randomly select 1,000 trees from the 10,000 Hackett trees
# ------------------------------------------------------------------
# NOTE: Original tree files were downloaded from https://birdtree.org/
#       and are not included in the repository due to size.

# Sample 100 trees from every 1k tree file (Warning: long runtime, ~30min)
tree1k <- list()
for (i in 1:10) {
  if (i == 1) {
    trees <- read.tree("AllBirdsHackett1.tre")
  } else {
    trees <- read.tree(paste0("BirdzillaHackett", i, ".tre"))
  }
  tree1k <- c(tree1k, trees[sample(1000, 100)])
}

saveRDS(tree1k, "data/tree1k.rds")

# ------------------------------------------------------------------------
# STEP 2: Trim trees and build a consensus tree
# ------------------------------------------------------------------------
# load trait and tree data
tree1k <- readRDS("data/raw/tree1k.rds")
trait <- read.csv("data/processed/trait_passeroid_pca.csv", row.names = 1)

# Drop non-overlapping species from each tree
passeroidea_tree_set <- lapply(tree1k, function(tree) {
  drop.tip(tree, setdiff(tree$tip.label, rownames(trait)))
})
saveRDS(passeroidea_tree_set, "data/processed/tree1k_pas.rds")

# Convert to multiPhylo object
class(passeroidea_tree_set) <- "multiPhylo"

# Compute consensus tree using least-squares branch lengths (2 mins)
consensus_tree <- consensus.edges(
  passeroidea_tree_set,
  method = "least.squares",
  consensus.tree = consensus(passeroidea_tree_set, p = 0.5, rooted = TRUE)
)
saveRDS(consensus_tree, "data/processed/contree_pas.rds")
