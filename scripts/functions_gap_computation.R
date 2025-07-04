# Functions: Gap computation
# ----------------------------------------------------

# Compute topological persistence (death - birth) of H1 features
get_topological_persistence <- function(tda) {
  h1 <- which(tda$diagram[, 1] == 1) # select only H1 features
  return(tda$diagram[h1, "Death"] - tda$diagram[h1, "Birth"])
}

get_gap_coordinates <- function(tda) {
  diag <- tda$diagram
  H1 <- which(diag[, 1] == 1)
  out <- lapply(H1, function(i) as.data.frame(tda$cycleLocation[[i]]))
  return(out)
}

compute_gap_centroid_size <- function(dt_coord) {
  vertices <- distinct(as.data.frame(rbind(
    as.matrix(dt_coord[, c(1, 3, 5, 7)]),
    as.matrix(dt_coord[, c(2, 4, 6, 8)]))))
  
  # compute centroid (mean coordinates)
  centroid <- colMeans(vertices)
  
  # compute gap size (mean Euclidean distance of vertices to centroid)
  dist_matrix <- matrix(centroid, ncol = 4, nrow = nrow(vertices), byrow = TRUE)
  size <- mean(sqrt(rowSums((vertices - dist_matrix)^2)))
  
  out <- list(centroid = centroid, size = size)
  return(out)
}

get_sparsity <- function (centroid, data) {
  distances <- sqrt(rowSums((data[, 1:4] - centroid[1:4])^2))
  sparsity <- sort(distances)[1:(ceiling(nrow(data)*0.05))] %>% mean
  return(sparsity)
}

get_gap_geometry_at_timepoint <- function(path_base, mya, THRES_PERSIST, THRES_SIZE) {
  trait <- read.csv(file.path(path_base, paste0("trait_", mya, "mya.csv")))
  tda <- readRDS(file.path(path_base, paste0("tda_", mya, "mya.rds")))
  persistence <- get_topological_persistence(tda)
  idx <- which(persistence > THRES_PERSIST) # filter gaps by persistence
  if (length(idx) == 0) return(NULL)
  
  coords <- get_gap_coordinates(tda)[idx] # get gap coordinates
  
  out <- bind_rows(lapply(seq_along(idx), function(i) { # for each hole
    temp <- compute_gap_centroid_size(coords[[i]])
    sparsity <- get_sparsity(temp$centroid, trait)
    data.frame(mya = mya, 
               n_species = nrow(coords[[i]]), # number of vertices (species)
               persistence = persistence[idx[i]], 
               size = temp$size, 
               sparsity = sparsity, 
               idx = paste0(mya, "_", idx[i]),
               pc1 = temp$centroid[1], pc2 = temp$centroid[2], pc3 = temp$centroid[3], pc4 = temp$centroid[4])
  })) %>% filter(size > THRES_SIZE)
  return(out)
}

link_gaps_across_time <- function(gaps, THRES_DIST) {
  # Step 1: Compute distances
  # compute pairwise distance (centroids) of all gaps across time
  D <- as.matrix(dist(gaps[, c("pc1", "pc2", "pc3", "pc4")]))
  D_bin <- D < THRES_DIST # binary pairwise matrix for whether gaps are nearby each other
  
  # Step 2: Initialize groupings
  gap_series <- as.list(seq_len(nrow(gaps)))
  
  # iterate over time points to link nearby gaps (skip first time point)
  time_points <- sort(unique(gaps$mya))[-1]
  
  # Step 3: Link gaps through time
  for (t in time_points) {
    prev <- which(gaps$mya == t - 1)  # gaps from previous time step
    curr <- which(gaps$mya == t)      # gaps from current time step
    if (length(prev) == 0 || length(curr) == 0) next  # skip if no gaps
    
    # link holes between consecutive time points
    for (j in prev) {
      connected <- curr[D_bin[j, curr]]
      gap_series[[j]] <- unique(c(gap_series[[j]], connected))
    }
  }
  
  # Step 4: Merge overlapping groups
  i <- 1
  while (i < length(gap_series)) {
    j <- i + 1
    while (j <= length(gap_series)) {
      if (length(intersect(gap_series[[i]], gap_series[[j]])) > 0) {
        gap_series[[i]] <- union(gap_series[[i]], gap_series[[j]])
        gap_series <- gap_series[-j]
      } else {
        j <- j + 1
      }
    }
    i <- i + 1
  }
  
  # Step 5: Sort by lifespan
  out <- gap_series[order(sapply(gap_series, length), decreasing = TRUE)]
  return(out)
}

get_gap_summary <- function(path_base, MAX_MYA, THRES_PERSIST, THRES_DIST, THRES_SIZE) {
  # compute hole centroid and size for all ancestors
  df <- bind_rows(lapply(0:MAX_MYA, function(mya) {
    get_gap_geometry_at_timepoint(path_base, mya, THRES_PERSIST, THRES_SIZE)
  }))
  if (nrow(df) == 0) return("No gaps found")
  
  # link gaps across time so gaps at similar locations in consecutive time point belong to the same gap series
  gap_series <- link_gaps_across_time(df, THRES_DIST)
  
  # assign index of gap series to each gap
  df$group <- factor(sapply(1:nrow(df), function(i) {
    match(TRUE, sapply(gap_series, function(g) i %in% g))
  }))
  
  # reorder group numbers to ensure 0 Mya gaps have the smallest numbers
  df$group <- factor(match(df$group, unique(df$group)))
  
  # calculate the persistence time of the longest persisting hole
  lifespan <- sapply(gap_series, function(g) df$mya[max(g)] - df$mya[min(g)])
  
  return(list(df_gaps = df, gap_lifespan = lifespan))
}

# ----------------------------------------------------
# Gap plotting
# ----------------------------------------------------
# Generates 3-panel visualization of gap distribution across time
plot_gap_history <- function(trait_data, mya, df_segments, gap_df, 
                             lim_pc1, lim_pc2, THRES_PERSIST, color = FALSE) {
  expand <- c(0.06, 0.06)
  gap_df <- arrange(gap_df, desc(group)) # adjust the plotting order (which points are on top of others)
  
  p_trait <- ggplot(trait_data, aes(x = pc1, y = pc2)) +
    geom_point(alpha = 0.4, size = 2, stroke = 0, color = "gray60") +
    geom_segment(df_segments,
                 mapping = aes(x = x1, xend = x2, y = y1, yend = y2, group = group, color = group),
                 linewidth = 1, lineend = "round") +
    scale_x_continuous(limits = lim_pc1, expand = expand) +
    scale_y_continuous(limits = lim_pc2, expand = expand) +
    theme_classic() +
    theme(panel.grid = element_blank(),
          axis.ticks.length = unit(-0.1, "cm"),
          axis.text.x = element_text(margin = margin(t = -12), size = 10),
          axis.text.y = element_text(margin = margin(r = -16), size = 10),
          axis.title = element_blank(),
          legend.position = "none",
          plot.margin = margin(0,0,0,0))
  
  p_x <- ggplot(gap_df, aes(x = pc1, y = mya, size = size, color = group, 
                            stroke = 3 / (1 - THRES_PERSIST) * (persistence - THRES_PERSIST) + 1.5)) +
    geom_point(shape = 1) +
    scale_size(range = c(1, 4)) +
    scale_x_continuous(limits = lim_pc1, expand = expand) +
    scale_y_reverse(limits = c(10, 0), expand = expand, breaks = c(10, 5, 0)) +
    theme_minimal() +
    theme(panel.grid = element_blank(),
          axis.title = element_blank(),
          axis.line.y = element_line(linewidth = 0.5),
          axis.ticks.y = element_line(linewidth = 0.5),
          axis.ticks.length.y = unit(-0.1, "cm"),
          axis.text.y = element_text(margin = margin(r = -16), hjust = 0, size = 10),
          axis.text.x = element_blank(),
          legend.position = "none",
          plot.margin = margin(0, 0, 0, 0))
  
  p_y <- ggplot(gap_df, aes(x = pc2, y = mya, size = size, color = group, 
                           stroke = 3 / (1 - THRES_PERSIST) * (persistence - THRES_PERSIST) + 1.5)) +
    geom_point(shape = 1) +
    scale_size(range = c(1, 4)) +
    scale_x_continuous(limits = lim_pc2, expand = expand) +
    scale_y_reverse(limits = c(10, 0), expand = expand, breaks = c(10, 5, 0)) +
    coord_flip() +
    theme_minimal() +
    theme(panel.grid = element_blank(),
          axis.title = element_blank(),
          axis.line.x = element_line(linewidth = 0.5),
          axis.ticks.x = element_line(linewidth =0.5),
          axis.ticks.length.x = unit(-0.1, "cm"),
          axis.text.x = element_text(margin = margin(t = -12), size = 10),
          axis.text.y = element_blank(),
          legend.position = "none",
          plot.margin = margin(0, 0, 0, 0))
  
  # set up color palettes
  if (color) {
    palette <- c("#cc4949","#435bab","#999999","#3a8fbd","#8a542c","#edba4c","#fa9b6b",
                 "#eba7b8","#9bd8fa","#d1ab92","#dde38f","#e07e2d","#ccb9ed","#dddddd")
  } else {
    set.seed(123)
    palette <- randomcoloR::distinctColorPalette(length(unique(gap_df$group)))
  }
  
  current_color <- rev(palette[gap_df$group[which(gap_df$mya == mya)]])
  p_trait <- p_trait + scale_color_manual(values = current_color)
  p_x <- p_x + scale_colour_manual(values = palette)
  p_y <- p_y + scale_colour_manual(values = palette)
  
  p <- plot_grid(p_y, p_trait, NULL, p_x, 
            rel_widths = c(1, 2.5, 2.5), rel_heights = c(1, 1, 1), ncol = 2,
            align = "hv", axis = "tblr")
  return(p)
}

# ---------------------------------------
# Wrapper that calls summary + plotting
# ---------------------------------------
plot_gap_history_print_summary <- function(dataset, mya = 0, MAX_MYA, THRES_PERSIST, THRES_DIST, THRES_SIZE, color = FALSE) {
  path_base <- file.path("data", dataset)
  
  # get hole information
  df_gaps <- get_gap_summary(path_base, MAX_MYA, THRES_PERSIST, THRES_DIST, THRES_SIZE)$df_gaps
  
  # import data
  tda <- readRDS(file.path(path_base, paste0("tda_", mya, "mya.rds")))
  data_modern <- read.csv(file.path(path_base, "trait_0mya.csv"))
  data_now <- read.csv(file.path(path_base, paste0("trait_", mya, "mya.csv")))
  
  ## prepare data for the scatter hole plot (top right plot)
  idx_mya <- df_gaps$idx[df_gaps$mya == mya] # index of gap appearing in x Mya
  
  # create data frame of all gap coordinates with columns indication gap id and size
  segment_data <- bind_rows(lapply(seq_along(idx_mya), function(i) {
    coords <- get_gap_coordinates(tda)[as.numeric(sub(".*_", "", idx_mya[i]))]
    as.data.frame(coords[[1]]) %>%
      mutate(idx = idx_mya[i], size = NA)
  })) %>%
    left_join(df_gaps %>% select(idx, group), by = "idx") %>%
    mutate(group = factor(group)) %>%
    setNames(c("x1","x2","y1","y2","z1","z2","w1","w2","idx","size","group"))
  
  # set axes limits
  lim_pc1 <- range(data_modern$pc1)
  lim_pc2 <- range(data_modern$pc2)
  
  p <- plot_gap_history(data_now, mya, segment_data, df_gaps, lim_pc1, lim_pc2, THRES_PERSIST, color)
  
  return(list(df_gaps = df_gaps, plot = p))
}
