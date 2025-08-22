#-------------------------------------------
# Phylogenetic relationship
#-------------------------------------------
dt <- allsp <- read.csv("data_input/AVONET_birdtree.csv") %>%
  filter(is.na(Reference.species)) %>% 
  filter(Order3 == "Passeriformes")
rownames(dt) <- gsub(" ", "_", dt$Species3)

tree <- readRDS("data_input/consensus_tree.rds") %>%
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

#-----------------------------------------
# Compare ecological traits
#-----------------------------------------
# Trait datasets
trait_pas_near <- dt[which(rownames(dt) %in% rownames(pas_near_gap)), ]
trait_nonpas_in <- dt[which(rownames(dt) %in% rownames(nonpas_within_gap)), ]

# South America datasets
trait_pas_near_sa <- trait_pas_near %>% filter(Centroid.Latitude < 15 & Centroid.Longitude > -82 & Centroid.Longitude < -30 &
                                                 !(Species3 %in% c("Chlorospingus inornatus")))
trait_nonpas_in_sa <- trait_nonpas_in %>% filter(Centroid.Latitude < 15 & Centroid.Longitude > -82 & Centroid.Longitude < -30)


sp.set <- rbind(
  cbind(trait_pas_near, group = "pas_near"),
  cbind(trait_nonpas_in, group = "nonpas"),
  cbind(trait_pas_near_sa, group = "pas_near_overlap"),
  cbind(trait_nonpas_in_sa, group = "nonpas_overlap")
) %>%
  mutate(Trophic.Niche = ifelse(!(Trophic.Niche %in% c("Frugivore","Granivore","Invertivore","Omnivore")), "Others", Trophic.Niche),
         Trophic.Niche = factor(Trophic.Niche, levels = c("Others","Granivore","Frugivore","Omnivore","Invertivore")),
         Territory = factor(Territory),
         Primary.Lifestyle = factor(Primary.Lifestyle, levels = c("Generalist","Aerial","Terrestrial","Insessorial")),
         # group = factor(group, levels = rev(c("nonpas","nonpas_overlap","pas_near_overlap","pas_near")))
         # group = factor(group, levels = rev(c("Passeriformes","nonpas","nonpas_overlap","pas_near_overlap","pas_near","Paseroidea")))
         # group = factor(group, levels = rev(c("nonpas","pas_near")))
         group = factor(group, levels = rev(c("nonpas_overlap","pas_near_overlap")))
  )

## tropical niche
# pdf("img/trait_trophic_6.pdf", width = 12, height = 3.2)
pdf("img/trait_trophic.pdf", width = 7, height = 1.2)
pdf("img/trait_trophic_southamerica.pdf", width = 7, height = 1.2)
sp.set %>% filter(!is.na(Trophic.Niche)) %>% 
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

## primary lifestyle
# pdf("img/trait_lifestyle_6.pdf", width = 12, height = 3.2)
pdf("img/trait_lifestyle.pdf", width = 7, height = 1.2)
sp.set %>% filter(!is.na(Primary.Lifestyle)) %>% 
  count(Primary.Lifestyle, group) %>%
  ggplot(aes(x = group, y = n, group = Primary.Lifestyle, fill = Primary.Lifestyle)) +
  geom_bar(stat = "identity", position = "fill") +
  geom_text(aes(label = n), position = position_fill(vjust = 0.5), size = 4, color = "white") + 
  scale_x_discrete(expand = c(0, 0)) +
  # scale_fill_manual(values = c("gray55", "#805132", "#90d184")) +
  # scale_fill_manual(values = c("gray55", "#98bbed", "#805132", "#90d184")) +
  scale_fill_manual(values = c("gray65", "#78583d", "#9cc977")) +
  # scale_fill_manual(values = c("gray50", "#4169a3", "#e39c62", "#215c3f")) +
  coord_flip() +
  ylab("") + xlab("") +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        text = element_text(size = 16))
dev.off()

## territoriality
# pdf("img/trait_territory_6.pdf", width = 12, height = 3.2)
pdf("img/trait_territory.pdf", width = 7, height = 1.2)
sp.set %>% filter(!is.na(Territory)) %>%
  count(Territory, group) %>%
  complete(Territory, group, fill = list(n = 0)) %>%
  ggplot(aes(x = group, y = n, group = Territory, fill = Territory)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_x_discrete(expand = c(0, 0)) +
  # scale_fill_manual(values = c("#6bb7ed", "#1b6294", "#03365c")) +
  # scale_fill_manual(values = c("gray70", "gray40", "gray25")) +
  # scale_fill_manual(values = c("#b8967f", "#59493e", "#29221d")) +
  # scale_fill_manual(values = c("#cca88b", "#78583d", "#573b24")) +
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

#### geographic distribution
library(ggforce)

# passeroids & non-passeroids overlaps
pdf("img/geo_distri.pdf", width = 10, height = 4.5)
ggplot() +
  geom_polygon(map_data("world"), mapping = aes(x = long, y = lat, group = group),
               fill = "gray90", color = "gray90", linewidth = 0) +
  geom_circle(sp.nonpas.hole, mapping = aes(x0 = Centroid.Longitude, y0 = Centroid.Latitude, r = sqrt(Range.Size/pi)/100),
              fill = alpha("cornflowerblue", 0.1), color = alpha("cornflowerblue", .8), linewidth = 0.3) +
  geom_circle(sp.hole.near, mapping = aes(x0 = Centroid.Longitude, y0 = Centroid.Latitude, r = sqrt(Range.Size/pi)/100), 
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

# All Passeroidea
ggplot(passeroidea.raw, aes(x = Centroid.Longitude, y = Centroid.Latitude, size = Range.Size)) +
  geom_point(shape = 21, fill = alpha("#cc4949", 0.02), color = "#cc4949", stroke = .2) +
  scale_size(range = c(1,20)) +
  theme_bw() +
  theme(panel.grid = element_blank())

a <- rbind(cbind(sp.hole.near, group = "pas_near"), 
           cbind(sp.nonpas.hole, group = "nonpas")) %>% 
  filter(Centroid.Longitude > -100 & 
           Centroid.Longitude < -40 & Centroid.Latitude < 10 & 
           Range.Size > 1000000)


## species trait data table
species.hole.df <- rbind(sp.nonpas.hole, sp.hole.near) %>% 
  select(Species3, Family3, Order3, Trophic.Niche, Primary.Lifestyle, Habitat) %>% 
  left_join(df.terri, by = join_by(Species3 == species))



