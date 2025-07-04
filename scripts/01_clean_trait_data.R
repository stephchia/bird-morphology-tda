# 01_clean_trait_data.R
# ---------------------
# This script loads and cleans AVONET morphological trait data.

# Load required packages
library(dplyr)
library(ggplot2)

# Filter species
allsp <- read.csv("data/raw/AVONET_birdtree.csv") %>%
  filter(is.na(Reference.species)) %>% # remove species whose traits were referenced from other species
  filter(!grepl("^Apteryx", Species3)) # remove kiwis (outliers)

# Filter to Passeroidea
passeroid_families <- c("Passeridae","Prunellidae","Urocynchramidae","Estrildidae","Ploceidae","Viduidae",
                        # Nine-primaried oscines
                        "Motacillidae","Peucedramidae","Fringillidae",
                        # Emberizoidea superfamily (New World nine-primaried oscines)
                        "Icteridae","Parulidae","Icteriidae","Phaenicophilidae","Zeledoniidae","Teretistridae","Thraupidae","Mitrospingidae",
                        "Rhodinocichlidae","Calyptophilidae","Nesospingidae","Spindalidae","Cardinalidae","Emberizidae","Passerellidae","Calcariidae")
passeroid <- allsp %>% filter(Family3 %in% passeroid_families)
non_passeroid <- allsp %>% filter(!(Family3 %in% passeroid_families))

# Set rownames and process data
process_trait_data <- function(data) {
  rownames(data) <- sub(" ", "_", data$Species3)
  data <- data %>% 
    select(Beak.Length_Culmen, Beak.Length_Nares, Beak.Width, Beak.Depth, Tarsus.Length, 
           Wing.Length, Secondary1, Hand.Wing.Index, Tail.Length, Mass) %>%
    mutate(across(-Hand.Wing.Index, log))
  return(data)
}

passeroid <- process_trait_data(passeroid) %>% scale
non_passeroid <- process_trait_data(non_passeroid)

# save processed passeroid trait data
saveRDS(passeroid, "data/processed/trait_passeroid.rds")
saveRDS(non_passeroid, "data/processed/trait_non_passeroid_raw.rds")

