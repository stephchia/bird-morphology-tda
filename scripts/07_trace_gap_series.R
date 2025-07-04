# 07_trace_gap_series.R
# ----------------------------------------------------
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
MYA <- 0

# plots gap history and print gap summary for empirical dataset
plot_gap_history_print_summary("empirical", MYA, MAX_MYA, THRES_PERSIST, THRES_DIST, THRES_SIZE, color = TRUE)
ggsave("output/gap_history.pdf", width = 10, height = 5)
