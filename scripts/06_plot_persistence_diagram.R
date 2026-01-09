# This script generates persistence diagrams and combined
# trait space + persistence diagram plots for empirical and simulated data.

library(ggplot2)
library(dplyr)
library(cowplot)

source("scripts/function_persistence_diagram.R")

# ------------------------------------------------------------
# Individual persistence diagrams for empirical data (at 0 and 5 Mya)
# ------------------------------------------------------------
cols <- c("gray50", "#4fbd9d", "#cc4949")
lim <- c(-.2, 2.45)

persistence_diag("empirical", mya = 0, cols, lim, idx = 468) # excluding one point with large y value (Figure 3D)
persistence_diag("empirical", mya = 5, cols, lim, idx = 283) # excluding one point with large y value (Figure 3E)

# ------------------------------------------------------------
# Trait space + persistent diagram combo plot (Figure S6)
# ------------------------------------------------------------
# empirical dataset
save_trait_pd_plot("empirical", mya = 0, "output/img_trait_pd/pd_emp.pdf")

# simulated datasets
for (i in 1:10) {
  dataset <- paste0("simulate/sim", i)
  filename <- paste0("output/img_trait_pd/pd_sim", i, ".pdf")
  save_trait_pd_plot(dataset, mya = 0, filename)
}
