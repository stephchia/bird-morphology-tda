# Recover original trait values from the focal topological gap; 
# identify and compare non-passeroid gap occupants and passeroids near the gap.

library(phytools)
library(ggtree)
library(dplyr)
library(tidyr)

source("scripts/functions_gap_computation.R")

# ---------------------------------------------
# Extract focal gap geometry
# ---------------------------------------------
tda <- readRDS("data/empirical/tda_0mya.rds")
fg_coords <- get_gap_coordinates(tda)[[468]] # the focal gap coordinates
centroid_size <- compute_gap_centroid_size(fg_coords)
centroid <- centroid_size$centroid
size <- centroid_size$size

# ---------------------------------------------
# Recover original trait values from centroid
# ---------------------------------------------
pca <- readRDS("data/processed/pca.rds")
passeroid_scaled <- readRDS("data/processed/trait_passeroid_raw.rds")
non_passeroid <- readRDS("data/processed/trait_non_passeroid_raw.rds")

# Inverse PCA and scaling to recover original (log-transformed) trait values, including range informed by gap size
gap_trait <- data.frame(centroid = t(centroid %*% t(pca$rotation[, 1:4])),
                        size = size * rowSums((pca$rotation[, 1:4])^2))
gap_trait$lower <- gap_trait$centroid - gap_trait$size/2
gap_trait$upper <- gap_trait$centroid + gap_trait$size/2
gap_trait_org_log <- data.frame(apply(gap_trait, 2, function(x) x * attr(passeroid_scaled, "scaled:scale") + attr(passeroid_scaled, "scaled:center")))
gap_trait_org_log$name <- rownames(gap_trait_org_log)

# Reverse log-transformation (except HWI)
trait_org <- gap_trait_org_log
trait_org[c(1:7,9,10), 1:4] <- exp(trait_org[c(1:7,9,10), 1:4])
trait_org 

# Reverse scaling
passeroid <- t(apply(passeroid_scaled, 1, function(x) x * attr(passeroid_scaled, "scaled:scale") + attr(passeroid_scaled, "scaled:center")))
passeroid_plot <- reshape2::melt(passeroid)

# ---------------------------------------------
# Plot gap traits in relation to overall trait distribution (Figure 4A)
# ---------------------------------------------
plot_trait_violin <- function(gap_trait_org_log, trait, log = TRUE) {
  gap_log <- gap_trait_org_log[trait, , drop = FALSE]
  gap_trait <- data.frame(x = factor(rownames(gap_log), levels = rownames(gap_log)), y = gap_log$centroid)
  
  df <- filter(passeroid_plot, Var2 %in% gap_trait$x)
  df$Var2 <- factor(df$Var2, levels = levels(gap_trait$x))

  breaks <- c(2.5, 5, 10, 25, 50, 100, 250, 500)

  p <- ggplot(df, aes(x = Var2, y = value)) +
  geom_violin(fill = "#3eb08f", width = 0.8, size = 0) +
  geom_crossbar(gap_log, mapping = aes(x = name, y = centroid, ymin = lower, ymax = upper), 
                width = 0.8, size = 0, fill = "white", alpha = 0.7, inherit.aes = FALSE) +
  geom_violin(df, mapping = aes(x = Var2, y = value),
              color = "#3eb08f", fill = NA, width = 0.8, size = 0.3) +
  geom_crossbar(gap_log, mapping = aes(x = name, y = centroid, ymin = centroid, ymax = centroid), 
                width = 0.8, size = 0.3, color = "black", inherit.aes = FALSE) +
  coord_flip() +
  theme_classic() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.title = element_blank(),
        axis.ticks.x = element_line(),
        axis.text.y = element_blank())
  if (log) p <- p + scale_y_log10(breaks = log(breaks), labels = breaks)
  return(p)
}

p1 <- plot_trait_violin(gap_trait_org_log, c("Mass"))
p2 <- plot_trait_violin(gap_trait_org_log, c("Tarsus.Length", "Beak.Depth", "Beak.Width", "Beak.Length_Nares", "Beak.Length_Culmen"))
p3 <- plot_trait_violin(gap_trait_org_log, c("Tail.Length", "Secondary1", "Wing.Length"))
p4 <- plot_trait_violin(gap_trait_org_log, c("Hand.Wing.Index"), log = FALSE)
p <- cowplot::plot_grid(nrow = 4, rel_heights = c(2, 6, 4, 2), p1, p2, p3, p4)
ggsave(filename = "output/gap_trait.pdf", plot = p, width = 3, height = 4, device = pdf)

# ---------------------------------------------
# Project non-passeroid species into PC space
# ---------------------------------------------
passeroid_pca <- read.csv("data/processed/trait_passeroid_pca.csv", row.name = 1)

# Apply the same scaling and projection to non-passeroids
non_passeroid_scaled <- scale(non_passeroid,
                              center = attr(passeroid_scaled, "scaled:center"),
                              scale = attr(passeroid_scaled, "scaled:scale"))
# non_passeroid_scaled <- t(apply(non_passeroid, 1, function(x) (x - attr(passeroid_scaled, "scaled:center")) / attr(passeroid_scaled, "scaled:scale")))
non_passeroid_pca <- as.data.frame(as.matrix(non_passeroid_scaled) %*% pca$rotation[, 1:4])

# --------------------------------------------------------
# Find species within or near the focal gap
# --------------------------------------------------------
# Non-passeroids within gap (distance < size/2) (species in Table S3)
dist_centroid_nonpas <- sapply(1:nrow(non_passeroid_pca), function(x) sqrt(sum((non_passeroid_pca[x, ] - centroid)^2)))
nonpas_within_gap <- non_passeroid_pca[which(dist_centroid_nonpas < size/2), ]

# Passeroids near gap (distance < size) (species in Table S4)
dist_centroid_pas <- sapply(1:nrow(passeroid_pca), function(x) sqrt(sum((passeroid_pca[x, ] - centroid)^2)))
pas_near_gap <- passeroid_pca[which(dist_centroid_pas < size), ]

#-------------------------------------------
# Process data for trait comparison
#-------------------------------------------
## All species
dt <- allsp <- read.csv("data/raw/AVONET_birdtree.csv")
rownames(dt) <- gsub(" ", "_", dt$Species3)

## Territoriality
# import data
terri.raw <- read.csv("data/raw/Tobias_etal_2016_territoriality.csv")

# match scientific names
terri.raw$Species[terri.raw$Species == "Emberiza lathami"] <- "Melophus lathami"
terri.raw$Species[terri.raw$Species == "Calendulauda barlowi"] <- "Certhilauda barlowi"
terri.raw$Species[terri.raw$Species == "Calendulauda erythrochlamys"] <- "Certhilauda erythrochlamys"
terri.raw$Species[terri.raw$Species == "Automolus roraimae"] <- "Syndactyla roraimae"
terri.raw$Species[terri.raw$Species == "Symposiarchus verticalis"] <- "Monarcha verticalis"
terri.raw$Species[terri.raw$Species == "Acritillas indica"] <- "Iole indica"
terri.raw$Species[terri.raw$Species == "Iole palawanensis"] <- "Ixos palawanensis"

passeroid_families <- c("Passeridae","Prunellidae","Urocynchramidae","Estrildidae","Ploceidae","Viduidae",
                        # Nine-primaried oscines
                        "Motacillidae","Peucedramidae","Fringillidae",
                        # Emberizoidea superfamily (New World nine-primaried oscines)
                        "Icteridae","Parulidae","Icteriidae","Phaenicophilidae","Zeledoniidae","Teretistridae","Thraupidae","Mitrospingidae",
                        "Rhodinocichlidae","Calyptophilidae","Nesospingidae","Spindalidae","Cardinalidae","Emberizidae","Passerellidae","Calcariidae")

terri <- data.frame(species = gsub("_", " ", c(rownames(pas_near_gap), rownames(nonpas_within_gap)))) %>%
  left_join(terri.raw, by = join_by(species == Species)) %>%
  select(species, Order, Family, Territory)

## Trait comprison datasets
trait_pas_near <- dt[which(rownames(dt) %in% rownames(pas_near_gap)), ] %>%
  left_join(terri, by = join_by(Species3 == species))
trait_nonpas_in <- dt[which(rownames(dt) %in% rownames(nonpas_within_gap)), ] %>%
  left_join(terri, by = join_by(Species3 == species))

## Trait comprison datasets (South America)
trait_pas_near_sa <- trait_pas_near %>% 
  filter(Centroid.Latitude < 15 & Centroid.Longitude > -82 & Centroid.Longitude < -30 &
    !(Species3 %in% c("Chlorospingus inornatus")))
trait_nonpas_in_sa <- trait_nonpas_in %>% 
  filter(Centroid.Latitude < 15 & Centroid.Longitude > -82 & Centroid.Longitude < -30)

#-------------------------------------------
# Species trait data table (Table S3 & S4)
#-------------------------------------------
# Table S3
select(trait_nonpas_in, Species3, Family3, Order3, Trophic.Niche, Primary.Lifestyle, Habitat, Territory)
# Table S4
select(trait_pas_near, Species3, Family3, Order3, Trophic.Niche, Primary.Lifestyle, Habitat, Territory)

#-------------------------------------------
# Geographic distribution (Figure 6A)
#-------------------------------------------
library(ggforce)

# passeroids & non-passeroids overlaps
pdf("output/geo_distri.pdf", width = 10, height = 4.5)
ggplot() +
  geom_polygon(map_data("world"), mapping = aes(x = long, y = lat, group = group),
               fill = "gray90", color = "gray90", linewidth = 0) +
  geom_circle(trait_nonpas_in, mapping = aes(x0 = Centroid.Longitude, y0 = Centroid.Latitude, r = sqrt(Range.Size/pi)/100),
              fill = alpha("cornflowerblue", 0.1), color = alpha("cornflowerblue", .8), linewidth = 0.3) +
  geom_circle(trait_pas_near, mapping = aes(x0 = Centroid.Longitude, y0 = Centroid.Latitude, r = sqrt(Range.Size/pi)/100), 
              fill = alpha("#cc4949", 0.1), color = alpha("#cc4949", .8), linewidth = 0.3) +
  # scale_size(range = c(1,20)) +
  coord_fixed(xlim = c(-180, 180), ylim = c(-62, 88), expand = 0) +
  xlab("Longitude") + ylab("Latitude") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        text = element_text(size = 12),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())
dev.off()

#-----------------------------------------
# Compare ecological traits (Figure 6B)
#-----------------------------------------
process_eco_trait <- function(data) {
  mutate(data, 
         Trophic.Niche = ifelse(!(Trophic.Niche %in% c("Frugivore","Granivore","Invertivore","Omnivore")), "Others", Trophic.Niche),
         Trophic.Niche = factor(Trophic.Niche, levels = c("Others","Granivore","Frugivore","Omnivore","Invertivore")),
         Territory = factor(Territory),
         Primary.Lifestyle = factor(Primary.Lifestyle, levels = c("Generalist","Aerial","Terrestrial","Insessorial")))
}

trait_compare <- rbind(cbind(trait_pas_near, group = "pas_near"),
                       cbind(trait_nonpas_in, group = "nonpas")) %>% 
  process_eco_trait %>%
  mutate(group = factor(group, levels = rev(c("nonpas","pas_near"))))

trait_compare_sa <- rbind(cbind(trait_pas_near_sa, group = "pas_near_overlap"),
                          cbind(trait_nonpas_in_sa, group = "nonpas_overlap")) %>%
  process_eco_trait %>%
  mutate(group = factor(group, levels = rev(c("nonpas_overlap","pas_near_overlap"))))

## tropical niche
# Figure 6B1
dir.create("output/eco_trait", showWarnings = FALSE)
pdf("output/eco_trait/trophic.pdf", width = 7, height = 1.2)
trait_compare %>% filter(!is.na(Trophic.Niche)) %>% 
  count(Trophic.Niche, group) %>%
  ggplot(aes(x = group, y = n, group = Trophic.Niche, fill = Trophic.Niche)) +
  geom_bar(stat = "identity", position = "fill") +
  geom_text(aes(label = n), position = position_fill(vjust = 0.5), size = 4, color = "white") + 
  scale_x_discrete(expand = c(0, 0)) +
  scale_fill_manual(values = c("gray45","#f2c76b","#d9762b","gray65","#2f6e4e")) +
  coord_flip() +
  ylab("") + xlab("") +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        text = element_text(size = 16))
dev.off()

# Figure 6B2
pdf("output/eco_trait/trophic_southamerica.pdf", width = 7, height = 1.2)
trait_compare_sa %>% filter(!is.na(Trophic.Niche)) %>% 
  count(Trophic.Niche, group) %>%
  ggplot(aes(x = group, y = n, group = Trophic.Niche, fill = Trophic.Niche)) +
  geom_bar(stat = "identity", position = "fill") +
  geom_text(aes(label = n), position = position_fill(vjust = 0.5), size = 4, color = "white") + 
  scale_x_discrete(expand = c(0, 0)) +
  scale_fill_manual(values = c("gray45","#f2c76b","#d9762b","gray65","#2f6e4e")) +
  coord_flip() +
  ylab("") + xlab("") +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        text = element_text(size = 16))
dev.off()

## primary lifestyle (Figure 6B3)
pdf("output/eco_trait/lifestyle.pdf", width = 7, height = 1.2)
trait_compare %>% filter(!is.na(Primary.Lifestyle)) %>% 
  count(Primary.Lifestyle, group) %>%
  ggplot(aes(x = group, y = n, group = Primary.Lifestyle, fill = Primary.Lifestyle)) +
  geom_bar(stat = "identity", position = "fill") +
  geom_text(aes(label = n), position = position_fill(vjust = 0.5), size = 4, color = "white") + 
  scale_x_discrete(expand = c(0, 0)) +
  scale_fill_manual(values = c("gray65", "#78583d", "#9cc977")) +
  coord_flip() +
  ylab("") + xlab("") +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        text = element_text(size = 16))
dev.off()

## territoriality (Figure 6B4)
pdf("output/eco_trait/territory.pdf", width = 7, height = 1.2)
trait_compare %>% filter(!is.na(Territory)) %>%
  count(Territory, group) %>%
  complete(Territory, group, fill = list(n = 0)) %>%
  ggplot(aes(x = group, y = n, group = Territory, fill = Territory)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_x_discrete(expand = c(0, 0)) +
  scale_fill_manual(values = c("#83c9c4", "#4a918d", "#1f4d49")) +
  geom_text(aes(label = n), position = position_fill(vjust = 0.5), size = 4, color = "white") + 
  coord_flip() +
  xlab("") + ylab("") +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        text = element_text(size = 16))
dev.off()

#-------------------------------------------
# Phylogenetic relationship (Figure S4)
#-------------------------------------------
tree <- readRDS("data/raw/consensus_tree.rds") %>%
  drop.tip(.$tip.label[!.$tip.label %in% rownames(dt)])
is_tip <- tree$edge[, 2] <= length(tree$tip.label)
ordered_tips <- tree$edge[is_tip, 2]

heat <- dt[match(tree$tip.label[ordered_tips], rownames(dt)), ]
heat$group <- "NA"
heat$group[rownames(heat) %in% rownames(pas_near_gap)] <- "passeroid_near_gap"
heat$group[rownames(heat) %in% rownames(nonpas_within_gap)] <- "nonpasseroid_within_gap"

# plot phylogenetic tree with highlighted passeroid_near_gap (red) and nonpasseroid_within_gap (blue)
p <- ggtree(tree, ladderize = F, layout = "circular", size = 0.07)
gheatmap(p, heat %>% select(group), width = 0.06, offset = -2, color = NULL, font.size = 0) +
  scale_fill_manual(values = c("gray95", "cornflowerblue", "#cc4949")) +
  ggnewscale::new_scale_fill()
