# Functions: Empirical and simulated gap comparison
# ------------------------------------------------------

plot_gap_stats_1D <- function(dt_sim, dt_emp, highlight = NULL, var, xlab, ymax) {
  p <- ggplot(dt_sim, aes_string(x = var)) +
    geom_density(color = "gray60", fill = "gray85", linewidth = 0.5) +
    geom_density(dt_emp, mapping = aes_string(x = var), 
                 color = "#3eb08f", fill = "#3eb08f", linewidth = 0.5, alpha = 0.3) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, ymax), expand = c(0, 0)) + 
    xlab(xlab) +
    ylab("Density") +
    theme_classic() +
    theme(panel.grid = element_blank())
  
  # add highlight arrow
  if (!is.null(highlight)) {
    ggdata <- ggplot_build(p)$data[[2]]
    dt.arrow <- data.frame(x = highlight[, var], 
                           yend = ggdata[which.min(abs(ggdata$x - highlight[, var])), "ymax"])
    p <- p + geom_segment(data = dt.arrow, aes(x = x, xend = x, y = yend + ymax/10, yend = yend + ymax/80), 
                          lineend = "round", color = "#cc4949", linewidth = .7,
                          arrow = arrow(length = unit(2, "mm"), type = "closed"))
  }
  return(p)
}

plot_gap_stats_2D <- function(dt_sim, dt_emp, highlight, x, y, xlab, ylab) {
  p <- ggplot(dt_sim, aes(x = .data[[x]], y = .data[[y]])) + 
    geom_hdr(probs = c(0.99, 0.95, 0.9, 0.75, 0.5)) + 
    geom_point(alpha = 0.15, stroke = 0) + 
    geom_point(dt_emp, mapping = aes(x = .data[[x]], y = .data[[y]]), 
               color = "#3eb08f", stroke = 0, alpha = 0.5, size = 1.5) + 
    geom_point(highlight, mapping = aes(x = .data[[x]], y = .data[[y]]), 
               color = "#cc4949", size = 2.5) + 
    xlab(xlab) + 
    ylab(ylab) + 
    theme_bw() + 
    theme(panel.grid = element_blank(), 
          legend.position = "none")
  return(p)
}

# plot gap series stats
plot_gap_series_2d <- function(data_sim, data_emp, x, y, xlab, ylab) {
  data_sim <- data_sim %>% group_by(evo_lifespan) %>%
    mutate(n = n(), show_vio = n > 5) %>% ungroup()
  
  p <- ggplot(filter(data_sim, show_vio == TRUE), aes(x = factor(.data[[x]]), y = .data[[y]])) +
    geom_violin(fill = "gray90", color = "gray50", width = 0.8, size = 0.2) +
    geom_point(data_sim, mapping = aes(x = factor(.data[[x]]), y = .data[[y]]), 
               color = "gray60", stroke = 0, size = 2.5, alpha = 0.8) +
    geom_point(data_emp, mapping = aes(x = factor(.data[[x]]), y = .data[[y]]), 
               shape = 21, fill = "#3eb08f", color = "#266e59", size = 3) +
    geom_point(filter(data_emp, evo_lifespan == 7), mapping = aes(x = factor(.data[[x]]), y = .data[[y]]), 
               shape = 21, fill = "#cc4949", color = "#9c3636", size = 4) +
    xlab(xlab) +
    ylab(ylab) +
    theme_bw() +
    theme(panel.grid = element_blank())
  return(p)
}
