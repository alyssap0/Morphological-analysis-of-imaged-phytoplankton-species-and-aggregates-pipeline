library(vegan)
library(dplyr)
library(readr)
library(parallelDist)
library(RColorBrewer)
library(factoextra)
library(ggrepel)
library(ggplot2)
library(viridis)
library(patchwork)  # For combining plots

# ---- Load and filter data (from your existing code) ----
agg_data <- read.csv("agg_data.csv")

agg_filtered <- agg_data %>%
  mutate(
    Porosity_percent = Porosity * 100,
    Wobbliness            = Perimeter_mm / Equivalent_Spherical_Diameter_mm,
    Fractal_Dimension = 2 * log10(Perimeter_mm) / log10(Area_mm2),
    Species = case_when(
      Species == "S1N100" ~ "Thalassionema sp.",
      Species == "M32C67" ~ "M32C68",
      TRUE ~ Species
    )
  ) %>%
  filter(!Microcosm.No. %in% c("M2_D97C3.1", "M4_D96C4", "M9_D97C3.2"))

agg_filtered <- agg_filtered %>%
  group_by(Species) %>%
  mutate(
    Replicate_No = as.numeric(factor(Microcosm.No.)),
    Species_Replicate = case_when(
      is.na(Species) | Species == "" ~ Microcosm.No.,
      n_distinct(Microcosm.No.) > 1 ~ paste0(Species, "_Rep", Replicate_No),
      TRUE ~ Species
    )
  ) %>%
  ungroup()

  

# Calculate means per replicate
agg_mean <- agg_filtered %>%
  group_by(Species_Replicate) %>%
  summarise(across(c(Equivalent_Spherical_Diameter_mm, Area_mm2, 
                     Circularity, Aspect_Ratio, Wobbliness, Perimeter_mm, Fractal_Dimension, Elongation, Porosity_percent), 
                   ~mean(.x, na.rm = TRUE)), .groups = "drop")

# Add Species column (extract base species name)
agg_mean <- agg_mean %>%
  mutate(
    Species = case_when(
      grepl("Chaetoceros", Species_Replicate) ~ "Chaetoceros sp.",
      grepl("Ditylum", Species_Replicate) ~ "Ditylum sp.",
      grepl("Melosira", Species_Replicate) ~ "Melosira sp.",
      grepl("Skeletonema", Species_Replicate) ~ "Skeletonema sp.",
      grepl("Thalassionema", Species_Replicate) ~ "Thalassionema sp.",
      TRUE ~ Species_Replicate  # Mixed assemblages
    ),
    Species_Abbrev = case_when(
      grepl("Chaetoceros", Species_Replicate) ~ gsub("Chaetoceros sp.", "C", Species_Replicate),
      grepl("Ditylum", Species_Replicate) ~ gsub("Ditylum sp.", "D", Species_Replicate),
      grepl("Melosira", Species_Replicate) ~ gsub("Melosira sp.", "M", Species_Replicate),
      grepl("Skeletonema", Species_Replicate) ~ gsub("Skeletonema sp.", "S", Species_Replicate),
      grepl("Thalassionema", Species_Replicate) ~ gsub("Thalassionema sp.", "T", Species_Replicate),
      TRUE ~ Species_Replicate
    ),
    Type = ifelse(grepl("_Rep", Species_Replicate), "Monoculture", "Mixed")
  )

# ---- Variables and pretty names ----
wanted_vars <- c("Equivalent_Spherical_Diameter_mm","Area_mm2", "Perimeter_mm",
                 "Circularity","Elongation","Aspect_Ratio","Porosity_percent","Wobbliness","Fractal_Dimension")

name_lookup <- c(
  Equivalent_Spherical_Diameter_mm = "ESD (mm)",
  Perimeter_mm = "Perimeter (mm)",
  Area_mm2 = "Area (mm²)",
  Aspect_Ratio = "Aspect Ratio",
  Circularity = "Circularity",
  Elongation = "Elongation",
  Porosity_percent = "2D-porosity (%)",
  Wobbliness = "Wobbliness",
  Fractal_Dimension = "Fractal Dimension"
)

# ---- Run ONE PCA on all data ----
X <- as.data.frame(agg_mean %>% select(all_of(wanted_vars)))
colnames(X) <- name_lookup[colnames(X)]

pca_all <- prcomp(X, center = TRUE, scale. = TRUE)
pc_var <- (pca_all$sdev^2) / sum(pca_all$sdev^2) * 100

# Add PC scores to data
agg_mean$PC1 <- pca_all$x[,1]
agg_mean$PC2 <- pca_all$x[,2]

# ---- Prepare data for LEFT panel (monoculture replicates) ----
left_data <- agg_mean %>% filter(Type == "Monoculture")

# ---- Prepare data for RIGHT panel (mixtures + monoculture means) ----
# Calculate monoculture species means
mono_means <- agg_mean %>%
  filter(Type == "Monoculture") %>%
  group_by(Species) %>%
  summarise(
    PC1 = mean(PC1),
    PC2 = mean(PC2),
    .groups = "drop"
  ) %>%
  mutate(
    Species_Abbrev = case_when(
      Species == "Chaetoceros sp." ~ "C",
      Species == "Ditylum sp." ~ "D",
      Species == "Melosira sp." ~ "M",
      Species == "Skeletonema sp." ~ "S",
      Species == "Thalassionema sp." ~ "T"
    ),
    Type = "Monoculture_Mean"
  )

# Mixtures
mixed_data <- agg_mean %>% filter(Type == "Mixed")

# Combine for right panel
right_data <- bind_rows(mono_means, mixed_data)

# ---- Create color and shape mappings ----
# Get all unique species/treatments for consistent mapping
all_treatments_left <- unique(left_data$Species_Abbrev)
all_treatments_right <- unique(right_data$Species_Abbrev)

# Create consistent color mapping using viridis
n_colors <- max(length(all_treatments_left), length(all_treatments_right))
viridis_colors <- viridis(n_colors, option = "D")

# Shapes
open_shapes <- c(15,16,0,1,2,3,4,5,6,17,7,8,9,10,11,12,13,18,19)

# LEFT panel mappings
color_map_left <- setNames(viridis_colors[1:length(all_treatments_left)], all_treatments_left)
shape_map_left <- setNames(rep(open_shapes, length.out = length(all_treatments_left)), all_treatments_left)

# RIGHT panel mappings
color_map_right <- setNames(viridis_colors[1:length(all_treatments_right)], all_treatments_right)
shape_map_right <- setNames(rep(open_shapes, length.out = length(all_treatments_right)), all_treatments_right)

# ---- Extract loadings for arrows ----
loadings <- as.data.frame(pca_all$rotation[, 1:2])
loadings$Variable <- rownames(loadings)
arrow_scale <- 4  # Adjust as needed


# ---- LEFT PLOT: Monoculture replicates ----
plot_left <- ggplot(left_data, aes(x = PC1, y = PC2, color = Species_Abbrev, shape = Species_Abbrev)) +
  geom_point(size = 4.5) +
  geom_segment(data = loadings, 
               aes(x = 0, y = 0, xend = PC1 * arrow_scale, yend = PC2 * arrow_scale),
               arrow = arrow(length = unit(0.3, "cm")), 
               color = "black", inherit.aes = FALSE) +
  geom_text_repel(data = loadings,
                  aes(x = PC1 * arrow_scale, y = PC2 * arrow_scale, label = Variable),
                  color = "black", size = 5, inherit.aes = FALSE) +
  scale_color_manual(values = color_map_left, name = "Treatment") +
  scale_shape_manual(values = shape_map_left, name = "Treatment") +
  labs(x = paste0("PC1 (", round(pc_var[1], 1), "%)"),
       y = paste0("PC2 (", round(pc_var[2], 1), "%)"),
       tag = "(a)") +
  theme_minimal(base_size = 16) +
  theme(
    axis.title.x = element_text(size = 20, face = "bold", margin = margin(t = 6)),
    axis.title.y = element_text(size = 20, face = "bold", margin = margin(r = 6)),
    axis.text = element_text(size = 16),
    legend.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 14),
    legend.key.size = unit(0.8, "cm"),
    panel.grid = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA),
    plot.tag = element_text(size = 24, face = "bold")
  )

# ---- RIGHT PLOT: Mixtures + Monoculture means ----
plot_right <- ggplot(right_data, aes(x = PC1, y = PC2, color = Species_Abbrev, shape = Species_Abbrev)) +
  geom_point(size = 4.5) +
  geom_segment(data = loadings, 
               aes(x = 0, y = 0, xend = PC1 * arrow_scale, yend = PC2 * arrow_scale),
               arrow = arrow(length = unit(0.3, "cm")), 
               color = "black", inherit.aes = FALSE) +
  geom_text_repel(data = loadings,
                  aes(x = PC1 * arrow_scale, y = PC2 * arrow_scale, label = Variable),
                  color = "black", size = 5, inherit.aes = FALSE) +
  scale_color_manual(values = color_map_right, name = "Treatment") +
  scale_shape_manual(values = shape_map_right, name = "Treatment") +
  labs(x = paste0("PC1 (", round(pc_var[1], 1), "%)"),
       y = paste0("PC2 (", round(pc_var[2], 1), "%)"),
       tag = "(b)") +
  theme_minimal(base_size = 16) +
  theme(
    axis.title.x = element_text(size = 20, face = "bold", margin = margin(t = 6)),
    axis.title.y = element_text(size = 20, face = "bold", margin = margin(r = 6)),
    axis.text = element_text(size = 16),
    legend.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 14),
    legend.key.size = unit(0.8, "cm"),
    panel.grid = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA),
    plot.tag = element_text(size = 24, face = "bold")
  )

# ---- Combine plots side by side ----
combined_plot <- plot_left + plot_right + plot_layout(ncol = 2)

# Save
ggsave("Fig8_PCA_SideBySide.jpg", combined_plot, width = 24, height = 12, dpi = 300)


# ---- PERMANOVA on morphological data ----

# Prepare the data matrix
morph_matrix <- agg_mean %>% 
  select(all_of(wanted_vars)) %>%
  scale()

# Calculate distance matrix
morph_dist <- dist(morph_matrix, method = "euclidean")

# Option 4: Nested design - Species nested within Type
# This tests if species differ within monocultures, and if mixtures differ
adonis_nested <- adonis2(morph_dist ~ Type/Species, 
                         data = agg_mean, 
                         permutations = 999)
print("PERMANOVA - Type and Species nested:")
print(adonis_nested)

# Check homogeneity of variances (important assumption!)
# Betadisper tests if groups have similar multivariate spread
betadisper_type <- betadisper(morph_dist, agg_mean$Type)
print("Betadisper - Type:")
print(anova(betadisper_type))

betadisper_species <- betadisper(morph_dist, agg_mean$Species)
print("Betadisper - Species:")
print(anova(betadisper_species))

# ---- PERMANOVA mirroring PCA panels ----
# Prepare the data matrix (same as what went into PCA)
morph_matrix <- agg_mean %>% 
  select(all_of(wanted_vars)) %>%
  scale()

# Calculate distance matrix
morph_dist <- dist(morph_matrix, method = "euclidean")
# Option 4: Nested design - Species nested within Type
# This tests if species differ within monocultures, and if mixtures differ
adonis_nested <- adonis2(morph_dist ~ Type/Species, 
                         data = agg_mean, 
                         permutations = 999)
print("PERMANOVA - Type and Species nested:")
print(adonis_nested)


##---- Species means (for your chosen workflow) ----
  agg_mean <- agg_filtered %>%
  group_by(Species) %>%
  summarise(across(all_of(vars), ~mean(.x, na.rm = TRUE)), .groups = "drop")

# Scale the means
morph_scaled <- scale(agg_mean %>% select(-Species))

# Distance + PERMANOVA (distance ~ Species with matching rows in agg_mean)
morph_dist <- dist(morph_scaled, method = "euclidean")
adonis_result <- adonis2(morph_dist ~ Species, data = agg_mean)
print(adonis_result)

cor(agg_mean$Perimeter_mm, agg_mean$Porosity_percent, use = "complete.obs")
