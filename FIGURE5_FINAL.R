# Load required packages
library(ggplot2)
library(dplyr)
library(tidyr)
library(reshape2)
library(scales)
library(patchwork)
library(viridis)

setwd("C:/Users/Cashies/OneDrive - University of Tasmania/AP_THESIS/Data")

# Read CSV file
all_data <- read.csv("agg_data.csv", header = TRUE)

all_data <- all_data %>%
  mutate(Porosity_correct = Total_Hole_Area_mm2 / Area_mm2)

# Selecting variables of interest and excluding dots
dat <- all_data %>%
  select(Species, Equivalent_Spherical_Diameter_mm, Perimeter_mm, Area_mm2, Circularity, Aspect_Ratio, Porosity_correct) %>%
  filter(Area_mm2 > 0, Perimeter_mm > 0) %>% 
  mutate(
    Species = as.factor(Species),
    Porosity_percent = Porosity_correct * 100,
    Wobble = Perimeter_mm / Equivalent_Spherical_Diameter_mm,
    Fractal_Dimension = 2 * log10(Perimeter_mm) / log10(Area_mm2)
  ) %>%
  filter(Species %in% c("Chaetoceros sp.", "Ditylum sp.", "Melosira sp.", "Skeletonema sp.", "Thalassionema sp."),
         (Porosity_percent > 0))

# Regression-derived FD by species
fd_by_species <- dat %>%
  group_by(Species) %>%
  do({
    model <- lm(log10(Perimeter_mm) ~ log10(Area_mm2), data = .)
    slope <- coef(model)["log10(Area_mm2)"]
    data.frame(Fractal_Dimension = 2 * slope)
  })

# Reshape the data to long format
dat_long <- melt(dat, id.vars = "Species")

#####ADD SIGNIFICANCE LETTERS 

# Define species levels in the order they appear in your data
species_levels <- c("Chaetoceros sp.", "Ditylum sp.", "Melosira sp.", "Skeletonema sp.", "Thalassionema sp.")

metrics_to_test <- c(
  "Equivalent_Spherical_Diameter_mm",
  "Perimeter_mm",
  "Area_mm2",
  "Circularity",
  "Aspect_Ratio",
  "Porosity_percent",
  "Wobble",
  "Fractal_Dimension"
)

# Create dataframe with all combinations
letters_grad <- expand.grid(
  Species = factor(species_levels, levels = species_levels),
  variable = metrics_to_test,
  stringsAsFactors = FALSE
)

# MANUAL ENTER OF ESD LETTERS 

letters_grad$Letter[ letters_grad$variable == "Equivalent_Spherical_Diameter_mm" & letters_grad$Species == "Chaetoceros sp." ] <- "a" 
letters_grad$Letter[ letters_grad$variable == "Equivalent_Spherical_Diameter_mm" & letters_grad$Species == "Ditylum sp." ] <- "b" 
letters_grad$Letter[ letters_grad$variable == "Equivalent_Spherical_Diameter_mm" & letters_grad$Species == "Melosira sp." ] <- "c" 
letters_grad$Letter[ letters_grad$variable == "Equivalent_Spherical_Diameter_mm" & letters_grad$Species == "Skeletonema sp." ] <- "a" 
letters_grad$Letter[ letters_grad$variable == "Equivalent_Spherical_Diameter_mm" & letters_grad$Species == "Thalassionema sp." ] <- "d" 

# MANUAL ENTER OF PERIMETER LETTERS 

letters_grad$Letter[ letters_grad$variable == "Perimeter_mm" & letters_grad$Species == "Chaetoceros sp." ] <- "a" 
letters_grad$Letter[ letters_grad$variable == "Perimeter_mm" & letters_grad$Species == "Ditylum sp." ] <- "b" 
letters_grad$Letter[ letters_grad$variable == "Perimeter_mm" & letters_grad$Species == "Melosira sp." ] <- "c" 
letters_grad$Letter[ letters_grad$variable == "Perimeter_mm" & letters_grad$Species == "Skeletonema sp." ] <- "d" 
letters_grad$Letter[ letters_grad$variable == "Perimeter_mm" & letters_grad$Species == "Thalassionema sp." ] <- "e" 

# MANUAL ENTER OF AREA LETTERS 

letters_grad$Letter[ letters_grad$variable == "Area_mm2" & letters_grad$Species == "Chaetoceros sp." ] <- "a" 
letters_grad$Letter[ letters_grad$variable == "Area_mm2" & letters_grad$Species == "Ditylum sp." ] <- "b" 
letters_grad$Letter[ letters_grad$variable == "Area_mm2" & letters_grad$Species == "Melosira sp." ] <- "c" 
letters_grad$Letter[ letters_grad$variable == "Area_mm2" & letters_grad$Species == "Skeletonema sp." ] <- "a" 
letters_grad$Letter[ letters_grad$variable == "Area_mm2" & letters_grad$Species == "Thalassionema sp." ] <- "d" 

# MANUAL ENTER OF Circularity LETTERS 

letters_grad$Letter[ letters_grad$variable == "Circularity" & letters_grad$Species == "Chaetoceros sp." ] <- "a" 
letters_grad$Letter[ letters_grad$variable == "Circularity" & letters_grad$Species == "Ditylum sp." ] <- "b" 
letters_grad$Letter[ letters_grad$variable == "Circularity" & letters_grad$Species == "Melosira sp." ] <- "c" 
letters_grad$Letter[ letters_grad$variable == "Circularity" & letters_grad$Species == "Skeletonema sp." ] <- "d" 
letters_grad$Letter[ letters_grad$variable == "Circularity" & letters_grad$Species == "Thalassionema sp." ] <- "e" 

# MANUAL ENTER OF Aspect Ratio LETTERS 

letters_grad$Letter[ letters_grad$variable == "Aspect_Ratio" & letters_grad$Species == "Chaetoceros sp." ] <- "a" 
letters_grad$Letter[ letters_grad$variable == "Aspect_Ratio" & letters_grad$Species == "Ditylum sp." ] <- "b" 
letters_grad$Letter[ letters_grad$variable == "Aspect_Ratio" & letters_grad$Species == "Melosira sp." ] <- "c" 
letters_grad$Letter[ letters_grad$variable == "Aspect_Ratio" & letters_grad$Species == "Skeletonema sp." ] <- "d" 
letters_grad$Letter[ letters_grad$variable == "Aspect_Ratio" & letters_grad$Species == "Thalassionema sp." ] <- "e" 

# MANUAL ENTER OF Porosity LETTERS 

letters_grad$Letter[ letters_grad$variable == "Porosity_percent" & letters_grad$Species == "Chaetoceros sp." ] <- "a" 
letters_grad$Letter[ letters_grad$variable == "Porosity_percent" & letters_grad$Species == "Ditylum sp." ] <- "b" 
letters_grad$Letter[ letters_grad$variable == "Porosity_percent" & letters_grad$Species == "Melosira sp." ] <- "b" 
letters_grad$Letter[ letters_grad$variable == "Porosity_percent" & letters_grad$Species == "Skeletonema sp." ] <- "c" 
letters_grad$Letter[ letters_grad$variable == "Porosity_percent" & letters_grad$Species == "Thalassionema sp." ] <- "d" 

# MANUAL ENTER OF Wobble LETTERS 

letters_grad$Letter[ letters_grad$variable == "Wobble" & letters_grad$Species == "Chaetoceros sp." ] <- "a" 
letters_grad$Letter[ letters_grad$variable == "Wobble" & letters_grad$Species == "Ditylum sp." ] <- "b" 
letters_grad$Letter[ letters_grad$variable == "Wobble" & letters_grad$Species == "Melosira sp." ] <- "c" 
letters_grad$Letter[ letters_grad$variable == "Wobble" & letters_grad$Species == "Skeletonema sp." ] <- "d" 
letters_grad$Letter[ letters_grad$variable == "Wobble" & letters_grad$Species == "Thalassionema sp." ] <- "e" 

# MANUAL ENTER OF Fractal Dimension LETTERS 

letters_grad$Letter[ letters_grad$variable == "Fractal_Dimension" & letters_grad$Species == "Chaetoceros sp." ] <- "a" 
letters_grad$Letter[ letters_grad$variable == "Fractal_Dimension" & letters_grad$Species == "Ditylum sp." ] <- "b" 
letters_grad$Letter[ letters_grad$variable == "Fractal_Dimension" & letters_grad$Species == "Melosira sp." ] <- "b" 
letters_grad$Letter[ letters_grad$variable == "Fractal_Dimension" & letters_grad$Species == "Skeletonema sp." ] <- "c" 
letters_grad$Letter[ letters_grad$variable == "Fractal_Dimension" & letters_grad$Species == "Thalassionema sp." ] <- "d" 

# Create letters_long for use in plotting functions
letters_long <- letters_grad

# ─────────────────────────────────────────────
# UPDATED Violin-plot function with proper clipping
# ─────────────────────────────────────────────
create_violinplot <- function(dat, metric, y_label, y_limits = NULL, tag = "",
                              show_x_labels = FALSE, letters_df = NULL, letter_size = 6,
                              add_points = FALSE, point_alpha = 0.35, point_size = 1.7) {
  dat_sub <- dat %>% filter(variable == metric)
  
  p <- ggplot(dat_sub, aes(x = Species, y = value, fill = Species)) +
    geom_violin(trim = FALSE, scale = "width", alpha = 0.8, color = "black") +
    # optional raw points
    { if (add_points) geom_jitter(width = 0.15, alpha = point_alpha, size = point_size) } +
    # median tick mark
    stat_summary(fun = median, geom = "point", shape = 95, size = 8, color = "black") +
    scale_fill_viridis_d(option = "D") +
    labs(y = y_label, x = NULL, tag = tag) +
    theme_bw() +
    theme(
      legend.position = "none",
      axis.text.x = if (show_x_labels) element_text(angle = 45, hjust = 1, vjust = 1, color = "black", size = 18) else element_blank(),
      axis.ticks.x = if (show_x_labels) element_line() else element_blank(),
      axis.text.y = element_text(color = "black"),
      axis.line   = element_line(color = "black"),
      text        = element_text(size = 20, color = "black"),
      panel.grid  = element_blank(),
      plot.tag    = element_text(face = "bold", size = 16),
      plot.margin = margin(t = 30, r = 10, b = 10, l = 10)
    )
 
  # Apply y-axis limits
  if (!is.null(y_limits)) {
    
    ymax <- y_limits[2]
    
    p <- p +
      scale_y_continuous(
        limits = y_limits,
        expand = c(0, 0)
      )
    
  } else {
    
    ymax <- max(dat_sub$value, na.rm = TRUE)
    
  }
  
  # Add letters if provided and significant for this metric
  if (!is.null(letters_df)) {
    letters_sub <- letters_df %>%
      filter(variable == metric) %>%
      mutate(Species = factor(Species, levels = levels(dat_sub$Species)))
    
    if (nrow(letters_sub) > 0) {
      # Position letters at the top of the y-limit range
      p <- p +
        geom_text(
          data = letters_sub,
          aes(x = Species, y = ymax, label = Letter),
          inherit.aes = FALSE,
          vjust = -0.6,
          size = letter_size,
          fontface = "bold",
          color = "black"
        ) +
        coord_cartesian(clip = "off")
      
      return(p)
    }
  }
  
  # If no letters, just apply coord_cartesian for consistency
  p <- p + coord_cartesian(clip = "off")
  p
}

# UPDATED Scatterplot function with matching letter positioning
create_scatterplot <- function(fd_data, y_label, y_limits = NULL, tag = "",
                               show_x_labels = FALSE, letters_df = NULL,
                               letter_size = 6) {
  
  p <- ggplot(fd_data, aes(x = Species, y = Fractal_Dimension, color = Species)) +
    geom_point(size = 5) +
    scale_color_viridis_d(option = "D") +
    labs(y = y_label, x = NULL, tag = tag) +
    theme_bw() +
    theme(
      legend.position = "none",
      axis.text.x = if (show_x_labels) element_text(angle = 45, hjust = 1, vjust = 1,
                                                    color = "black", size = 18) else element_blank(),
      axis.ticks.x = if (show_x_labels) element_line() else element_blank(),
      axis.text.y = element_text(color = "black"),
      axis.line   = element_line(color = "black"),
      text        = element_text(size = 20, color = "black"),
      panel.grid  = element_blank(),
      plot.tag    = element_text(face = "bold", size = 16),
      plot.margin = margin(t = 30, r = 10, b = 10, l = 10)
    )
  
  # Apply y-limits if specified
  if (!is.null(y_limits)) {
    ymin <- y_limits[1]
    ymax <- y_limits[2]
    p <- p + scale_y_continuous(limits = c(ymin, ymax), expand = c(0, 0))
  } else {
    ymin <- min(fd_data$Fractal_Dimension, na.rm = TRUE)
    ymax <- max(fd_data$Fractal_Dimension, na.rm = TRUE)
  }
  
  # Add letters if significant
  if (!is.null(letters_df)) {
    letters_sub <- letters_df %>%
      filter(variable == "Fractal_Dimension") %>%
      mutate(Species = factor(Species, levels = levels(fd_data$Species)))
    
    if (nrow(letters_sub) > 0) {
      # Position letters at the top of the y-limit range
      p <- p +
        geom_text(
          data = letters_sub,
          aes(x = Species, y = ymax, label = Letter),
          inherit.aes = FALSE,
          vjust = -0.6,
          size = letter_size,
          fontface = "bold",
          color = "black"
        ) +
        coord_cartesian(clip = "off")
      
      return(p)
    }
  }
  
  # If no letters, just apply coord_cartesian for consistency
  p <- p + coord_cartesian(clip = "off")
  p
}

# Build panels (now with proper clipping)
p1 <- create_violinplot(dat_long, "Equivalent_Spherical_Diameter_mm", "ESD (mm)", c(0, 4),  "(a)", letters_df = letters_long)
p2 <- create_violinplot(dat_long, "Perimeter_mm",                    "Perimeter (mm)", c(0, 20), "(b)", letters_df = letters_long)
p3 <- create_violinplot(dat_long, "Area_mm2",                        "Area (mm²)", c(0, 4),     "(c)", letters_df = letters_long)
p4 <- create_violinplot(dat_long, "Circularity",                     "Circularity", c(0, 1),    "(d)", letters_df = letters_long)
p5 <- create_violinplot(dat_long, "Aspect_Ratio",                    "Aspect Ratio", c(0, 6),   "(e)", letters_df = letters_long)
p6 <- create_violinplot(dat_long, "Porosity_percent",                "2D-Porosity (%)", c(0, 20),  "(f)", letters_df = letters_long)
p7 <- create_violinplot(dat_long, "Wobble",                          "Wobbliness", c(0, 15),        "(g)", show_x_labels = TRUE, letters_df = letters_long)

# FD scatter with updated function
p8 <- create_scatterplot(fd_by_species, "Fractal Dimension", c(0, 2), "(h)", show_x_labels = TRUE, letters_df = letters_long)

final_fig <- (p1 + p2 + p3 + p4 + p5 + p6 + p7 + p8) +
  plot_layout(ncol = 2)


ggsave("Fig5_Aggs_Violins.jpg", final_fig, width = 14, height = 12, dpi = 300)


