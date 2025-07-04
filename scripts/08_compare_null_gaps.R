# 08_compare_null_gaps.R
# ------------------------------
# This script compares the distribution of topological gaps
# between empirical and simulated trait evolution data.

library(dplyr)
library(ggplot2)
library(ggdensity)
library(cowplot)

# Load custom functions
source("scripts/functions_gap_computation.R")
source("scripts/functions_plot_null_comparison.R")

# Set parameters
THRES_PERSIST <- 0
THRES_DIST <- 1
THRES_SIZE <- 0
MAX_MYA <- 10

# -----------------------------
# Load and summarize gap data
# -----------------------------
# Empirical gaps
gap_stats_emp <- get_gap_summary(
  "data/empirical", MAX_MYA, THRES_PERSIST, THRES_DIST, THRES_SIZE
  )$df_gaps

# Simulated gaps (takes a minute)
gap_stats_sim <- bind_rows(lapply(1:10, function(sim_id) { # a few minutes
  get_gap_summary(
    paste0("data/simulate/sim", sim_id), MAX_MYA, THRES_PERSIST, THRES_DIST, THRES_SIZE
    )$df_gaps %>%
    mutate(sim = sim_id)
}))

# ----------------------------------------------------
# Visualize gap distributions at modern time (0 Mya)
# ----------------------------------------------------
gap_emp_now <- gap_stats_emp %>% filter(mya == 0)
gap_sim_now <- gap_stats_sim %>% filter(mya == 0)

# Highlight the focal gap for plotting (optional)
focal_gap <- gap_emp_now %>% filter(idx == "0_468")

# 1D density plots (Figure S5)
plot_gap_stats_1D(gap_sim_now, gap_emp_now, var = "persistence",xlab = "Topological persistence", ymax = 7)
ggsave("output/img_gap_stats/1d_persistence.pdf", width = 6, height = 2.5)

plot_gap_stats_1D(gap_sim_now, gap_emp_now, var = "size", xlab = "Gap size", ymax = 2.8)
ggsave("output/img_gap_stats/1d_size.pdf", width = 6, height = 2.5)

plot_gap_stats_1D(gap_sim_now, gap_emp_now, var = "sparsity", xlab = "Sparsity", ymax = 2.5)
ggsave("output/img_gap_stats/1d_sparsity.pdf", width = 6, height = 2.5)

# 2D plot: persistence vs. sparsity (Figure 4A)
plot_gap_stats_2D(gap_sim_now, gap_emp_now,
                  highlight = focal_gap, # Highlight the focal gap
                  x = "sparsity", y = "persistence",
                  xlab = "Sparsity (PC unit)", ylab = "Topological persistence (PC unit)")
ggsave("output/img_gap_stats/2d_persistence_sparsity.pdf", width = 5, height = 4.5)

# --------------------------------------------
# Analyze gap series lifespan vs. sparsity
# --------------------------------------------
THRES_PERSIST <- 0.4
THRES_DIST <- 1
THRES_SIZE <- 0
MAX_MYA <- 10

# Summarize each gap's trajectory over time
summarize_gap_series <- function(path_base, MAX_MYA, THRES_PERSIST, THRES_DIST, THRES_SIZE, dataset_label) {
  out <- get_gap_summary(path_base, MAX_MYA, THRES_PERSIST, THRES_DIST, THRES_SIZE)$df_gaps %>% 
    group_by(group) %>%
    summarise(n_species = mean(n_species),
              persistence = mean(persistence),
              size = mean(size),
              sparsity = mean(sparsity),
              evo_lifespan = max(mya) - min(mya)) %>%
    mutate(dataset = dataset_label)
  return(out)
}

# Summarize empirical gap series
gap_series_emp <- summarize_gap_series("data/empirical", MAX_MYA, THRES_PERSIST, THRES_DIST, THRES_SIZE, "empirical")

# Summarize simulated gap series
gap_series_sim <- bind_rows(lapply(1:10, function(sim_id) {
  summarize_gap_series(paste0("data/simulate/sim", sim_id), 
                       MAX_MYA, THRES_PERSIST, THRES_DIST, THRES_SIZE, paste0("sim", sim_id))
}))

# Plot: evolutionary lifespan vs. sparsity (Figure 4B)
plot_gap_series_2d(data_sim = gap_series_sim, data_emp = gap_series_emp,
                   x = "evo_lifespan", y = "sparsity",
                   xlab = "Evolutionary lifespan (million years)", ylab = "Sparsity (PC unit)")
ggsave("output/img_gap_stats/2d_sparsity_lifespan.pdf", width = 5, height = 4.5)
