
# Functions: Ancestral state reconstruction and simulation
# -------------------------------------------------------------

# Interpolate reconstructed ancestral traits at specific time point
get_traits_at_timepoint <- function(tree, tipnode, mya) {
  if (mya == 0) {
    return(tipnode[1:length(tree$tip.label), , drop = FALSE]) # Return tip values
  } 
  
  tree_depth <- max(nodeHeights(tree))
  time <- tree_depth - mya
  X <- nodeHeights(tree)
  
  # Determine which branches are crossed at this time slice
  idx <- which(X[, 1] < time & X[, 2] > time)
  
  # Interpolate traits for each branch
  interpolated <- t(sapply(idx, function(i) {
    edge <- tree$edge[i, ]
    trait_start <- tipnode[edge[1], ]
    trait_end   <- tipnode[edge[2], ]
    ratio <- (time - X[i, 1]) / (X[i, 2] - X[i, 1])
    trait_start + (trait_end - trait_start) * ratio
  }))
  
  colnames(interpolated) <- c("pc1","pc2","pc3","pc4")
  return(interpolated)
}

# Plot 1D trait ecolutionary trajectory (trait vs time)
plot_1D_trajectory <- function(tree, trait_vector, xbreaks = NULL, yspace = 1) {
  ylims <- range(trait_vector, na.rm = TRUE)
  par(mar = c(4, 4, 2, 0))
  
  phenogram(tree, trait_vector, fsize = 0, lwd = 0.7,
            spread.labels = FALSE, ylim = ylims,
            xlab = "Time since root (million years)", ylab = "Trait value")
  
  par(xaxt = "s", yaxt = "s")
  axis(side = 1, at = xbreaks)
  axis(side = 2, at = seq(floor(min(trait_vector)), ceiling(max(trait_vector)), yspace))
}
