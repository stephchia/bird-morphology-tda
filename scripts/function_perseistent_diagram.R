# Functions: Plot persistence diagram and trait distribution
# ------------------------------------------------------------

plot_trait_pc <- function(dataset, mya) {
  data <- read.csv(paste0("data/", dataset, "/trait_", mya, "mya.csv"))
  
  ggplot(data) +
    geom_point(aes(x = pc1, y = pc2), color = "gray20", alpha = 0.5) +
    scale_x_continuous(limits = c(-15.5, 15.5), expand = c(0, 0)) +
    scale_y_continuous(limits = c(-5.5, 5.5), expand = c(0, 0)) +
    xlab("PC1") +
    ylab("PC2") +
    theme_bw() +
    theme(panel.grid = element_blank())
}

persistence_diag <- function(dataset, mya, 
                             cols = c("cornflowerblue", "#fc9e3a", "brown3"), 
                             xlim = c(-.2, 3.5), ylim = c(-.2, 3.5), 
                             THRES_PERSISTENCE = 0.4) {
  
  tda <- readRDS(paste0("data/", dataset, "/tda_", mya, "mya.rds"))
  diag <- tda$diagram
  
  H0 <- which(diag[, 1] == 0)
  H1 <- which(diag[, 1] == 1)
  
  dt_H1 <- data.frame(birth = diag[H1, 2], death = diag[H1, 3], persistence = diag[H1, 3] - diag[H1, 2])
  dt_H0 <- data.frame(birth = diag[H0, 2], death = diag[H0, 3], persistence = diag[H0, 3] - diag[H0, 2])
  
  ggplot(dt_H1) +
    geom_point(dt_H0, mapping = aes(x = birth, y = death), 
               shape = 0, size = 1.8, stroke = 0.5, color = cols[1]) +
    geom_point(subset(dt.H1, persistence < THRES_PERSISTENCE), mapping = aes(x = birth, y = death), 
               shape = 21, size = 1.8, stroke = 0.5, color = cols[2]) +
    geom_point(subset(dt.H1, persistence > THRES_PERSISTENCE), mapping = aes(x = birth, y = death), 
               shape = 21, size = 1.8, stroke = 1, color = cols[3]) +
    geom_abline(intercept = 0, slope = 1, linewidth = 0.3) +
    scale_x_continuous(limits = xlim, expand = c(0, 0)) +
    scale_y_continuous(limits = ylim, expand = c(0, 0)) +
    xlab("Gap birth") + 
    ylab("Gap death") +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.text = element_text(size = 10))
}

save_trait_pd_plot <- function(dataset, mya, filename) {
  p1 <- plot_trait_pc(dataset, mya)
  p2 <- persistence_diag(dataset, mya)
  plot_grid(p1, p2, rel_widths = c(2, 1)) %>%
    save_plot(filename, ., base_height = 3, base_width = 9)
}
