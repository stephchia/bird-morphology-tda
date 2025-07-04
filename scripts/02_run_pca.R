# 02_run_pca.R
# ------------------
# This script performs PCA on Passeroidea, and saves PC scores for downstream analysis.

# Load required packages
library(dplyr)
library(ggplot2)
library(factoextra)

# Load trait data
passeroid <- readRDS("data/processed/trait_passeroid.rds")

# Run PCA
pca <- prcomp(passeroid, scale. = TRUE)

# Flip sign of selected loadings
flip_pc_sign <- function(pca, pc) {
  pca$x[, pc] <- -pca$x[, pc]
  pca$rotation[, pc] <- -pca$rotation[, pc]
  pca
}

if (pca$rotation["Mass", "PC1"] < 0) pca <- flip_pc_sign(pca, "PC1")
if (pca$rotation["Beak.Width", "PC2"] < 0) pca <- flip_pc_sign(pca, "PC2")
if (pca$rotation["Hand.Wing.Index", "PC3"] < 0) pca <- flip_pc_sign(pca, "PC3")
if (pca$rotation["Beak.Length_Culmen", "PC4"] < 0) pca <- flip_pc_sign(pca, "PC4")

# Add proportion of variance explained as the 11th row
pca_summary <- rbind(pca$rotation, 
                     "Variance (%)" = (pca$sdev)^2/sum((pca$sdev)^2)*100)
pca_summary
sum(pca_summary["Variance (%)", 1:4]) # total variance explained by first 4 PCs
saveRDS(pca, "data/processed/pca.rds")

# Save projected PCA scores
passeroid_pca <- as.data.frame(as.matrix(passeroid) %*% pca$rotation[, 1:4])
colnames(passeroid_pca) <- c("pc1", "pc2", "pc3", "pc4")
write.csv(passeroid_pca, "data/processed/trait_passeroid_pca.csv")

# PCA biplot (PC1 vs PC2)
pca_biplot <- fviz_pca_biplot(pca, label=c("var","quali"), repel = TRUE)
scaler <- 1.7
ggplot() +
  geom_segment(data = ggplot_build(pca_biplot)$data[[5]],
               mapping = aes(x = scaler * x, y = scaler * y, xend = scaler * xend, yend = scaler * yend),
               arrow = arrow(length = unit(0.2, "cm")), 
               color = "#75d9ad", linewidth = 0.6) +
  geom_point(data = passeroid_pca, 
             mapping = aes(x = pc1, y = pc2),
             color = "gray15", stroke = 0, alpha = 0.6, size = 1.4) +
  scale_x_continuous(limits = c(-6.5, 14.3)) +
  scale_y_continuous(limits = c(-4.5, 6), breaks = seq(-4, 6, 2)) +
  theme_classic() +
  theme(panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.line = element_line(linewidth = 0.2))
ggsave("output/trait_pca.pdf", width = 8, height = 4)
