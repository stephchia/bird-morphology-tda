# Visualize topological gaps traced across evolutionary time
# for empirical and simulated trait datasets.

library(dplyr)
library(phytools)
library(ggplot2)
library(cowplot)

# Load custom functions
source("scripts/functions_gap_computation.R")

# Parameters
THRES_PERSIST <- 0.4
THRES_DIST <- 1
THRES_SIZE <- 0
MAX_MYA <- 10

# plots gap history and print gap summary for empirical dataset (Figure 3)
for (MYA in 0:10) {
  plot_gap_history_print_summary("empirical", MYA, MAX_MYA, THRES_PERSIST, THRES_DIST, THRES_SIZE, color = TRUE)
  ggsave(paste0("output/gap_history/gap_history_", MYA, "mya.pdf"), width = 10, height = 5)
}
