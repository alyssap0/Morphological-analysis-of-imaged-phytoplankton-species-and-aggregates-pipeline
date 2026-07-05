# ============================================================
# 0. Load packages 
# ============================================================

library(ggplot2)
library(tidyr)
library(scales)
library(ggpubr)
library(patchwork)
library(viridis)
library(dplyr)

## ============================================================
## 1. SETTINGS — ONLY Mel-Chaet gradient in this script
## ============================================================
speciespct_levels <- c("0%", "15%", "51%", "54%", "59%", "63%", "68%", "94%", "100%")
species0_label   <- "Melosira sp."
species100_label <- "Chaetoceros sp."

all_data <- read.csv("agg_data.csv", header = TRUE)

## ============================================================
## 2. BASE DATA
## ============================================================
dat <- all_data %>%
  select(
    Species,
    Equivalent_Spherical_Diameter_mm,
    Perimeter_mm,
    Area_mm2,
    Circularity,
    Aspect_Ratio,
    Porosity
  ) %>%
  dplyr::filter(
    Area_mm2 > 0,
    Perimeter_mm > 0
  ) %>%
  dplyr::mutate(
    Porosity_Percent  = Porosity * 100,
    Wobbliness            = Perimeter_mm / Equivalent_Spherical_Diameter_mm,
    Fractal_Dimension = 2 * log10(Perimeter_mm) / log10(Area_mm2)) %>%
  dplyr::filter(Porosity_Percent > 0)

selected_species <- c(
  "Melosira sp.", "M85C15", "M49C51", "M46C54",
  "M41C59", "M37C63", "M32C67", "M6C94", "Chaetoceros sp."
)

dat2 <- dat %>%
  dplyr::filter(Species %in% selected_species) %>%
  dplyr::mutate(Species = factor(Species, levels = selected_species))


## ============================================================
## 3. CONVERT TO LONG FORMAT
## ============================================================
convert_gradient <- function(dat2) {
  recode_map <- c(
    "Melosira sp."   = "0%",
    "M85C15"             = "15%",
    "M49C51"            = "51%",
    "M46C54"            = "54%",
    "M41C59"            = "59%",
    "M37C63"            = "63%",
    "M32C67"            = "68%",
    "M6C94"            = "94%",
    "Chaetoceros sp." = "100%"
  )
  
  dat2 %>%
    mutate(
      SpeciesPct      = recode(Species, !!!recode_map),
      SpeciesPct      = factor(SpeciesPct, levels = speciespct_levels),
      Percent_numeric = as.numeric(gsub("%","", SpeciesPct))
    ) %>%
    pivot_longer(
      cols = -c(Species, SpeciesPct, Percent_numeric),
      names_to = "variable",
      values_to = "value"
    )
}

dat2_long <- convert_gradient(dat2)
levs_speciespct <- levels(dat2_long$SpeciesPct)

## ============================================================
## 4. PLACEHOLDER LETTERS (REPLACE LATER WITH REAL STATS)
## ============================================================
## ============================================================
## 4. PLACEHOLDER LETTERS (EDITABLE, WORKING VERSION)
## ============================================================

metrics_to_test <- c(
  "Equivalent_Spherical_Diameter_mm",
  "Perimeter_mm",
  "Area_mm2",
  "Circularity",
  "Aspect_Ratio",
  "Porosity_Percent",
  "Wobbliness",
  "Fractal_Dimension"
)

letters_grad <- expand.grid( SpeciesPct = factor(speciespct_levels, levels = speciespct_levels), variable = metrics_to_test, stringsAsFactors = FALSE ) 

# MANUAL ENTER OF ESD LETTERS 

letters_grad$Letter[ letters_grad$variable == "Equivalent_Spherical_Diameter_mm" & letters_grad$SpeciesPct == "0%" ] <- "a" 
letters_grad$Letter[ letters_grad$variable == "Equivalent_Spherical_Diameter_mm" & letters_grad$SpeciesPct == "15%" ] <- "a" 
letters_grad$Letter[ letters_grad$variable == "Equivalent_Spherical_Diameter_mm" & letters_grad$SpeciesPct == "51%" ] <- "c" 
letters_grad$Letter[ letters_grad$variable == "Equivalent_Spherical_Diameter_mm" & letters_grad$SpeciesPct == "54%" ] <- "d" 
letters_grad$Letter[ letters_grad$variable == "Equivalent_Spherical_Diameter_mm" & letters_grad$SpeciesPct == "59%" ] <- "e" 
letters_grad$Letter[ letters_grad$variable == "Equivalent_Spherical_Diameter_mm" & letters_grad$SpeciesPct == "63%" ] <- "f" 
letters_grad$Letter[ letters_grad$variable == "Equivalent_Spherical_Diameter_mm" & letters_grad$SpeciesPct == "68%" ] <- "g" 
letters_grad$Letter[ letters_grad$variable == "Equivalent_Spherical_Diameter_mm" & letters_grad$SpeciesPct == "94%" ] <- "h" 
letters_grad$Letter[ letters_grad$variable == "Equivalent_Spherical_Diameter_mm" & letters_grad$SpeciesPct == "100%" ] <- "b" 

#MANUAL ENTER OF PERIMETER LETTERS 
letters_grad$Letter[ letters_grad$variable == "Perimeter_mm" & letters_grad$SpeciesPct == "0%" ] <- "a" 
letters_grad$Letter[ letters_grad$variable == "Perimeter_mm" & letters_grad$SpeciesPct == "15%" ] <- "c" 
letters_grad$Letter[ letters_grad$variable == "Perimeter_mm" & letters_grad$SpeciesPct == "51%" ] <- "d" 
letters_grad$Letter[ letters_grad$variable == "Perimeter_mm" & letters_grad$SpeciesPct == "54%" ] <- "d" 
letters_grad$Letter[ letters_grad$variable == "Perimeter_mm" & letters_grad$SpeciesPct == "59%" ] <- "d" 
letters_grad$Letter[ letters_grad$variable == "Perimeter_mm" & letters_grad$SpeciesPct == "63%" ] <- "e" 
letters_grad$Letter[ letters_grad$variable == "Perimeter_mm" & letters_grad$SpeciesPct == "68%" ] <- "f" 
letters_grad$Letter[ letters_grad$variable == "Perimeter_mm" & letters_grad$SpeciesPct == "94%" ] <- "g" 
letters_grad$Letter[ letters_grad$variable == "Perimeter_mm" & letters_grad$SpeciesPct == "100%" ] <- "b" 
# MANUAL ENTER OF LETTERS AREA 
letters_grad$Letter[ letters_grad$variable == "Area_mm2" & letters_grad$SpeciesPct == "0%" ] <- "a" 
letters_grad$Letter[ letters_grad$variable == "Area_mm2" & letters_grad$SpeciesPct == "15%" ] <- "a" 
letters_grad$Letter[ letters_grad$variable == "Area_mm2" & letters_grad$SpeciesPct == "51%" ] <- "c" 
letters_grad$Letter[ letters_grad$variable == "Area_mm2" & letters_grad$SpeciesPct == "54%" ] <- "d" 
letters_grad$Letter[ letters_grad$variable == "Area_mm2" & letters_grad$SpeciesPct == "59%" ] <- "e" 
letters_grad$Letter[ letters_grad$variable == "Area_mm2" & letters_grad$SpeciesPct == "63%" ] <- "f" 
letters_grad$Letter[ letters_grad$variable == "Area_mm2" & letters_grad$SpeciesPct == "68%" ] <- "g" 
letters_grad$Letter[ letters_grad$variable == "Area_mm2" & letters_grad$SpeciesPct == "94%" ] <- "h" 
letters_grad$Letter[ letters_grad$variable == "Area_mm2" & letters_grad$SpeciesPct == "100%" ] <- "b" 

### CIRCULARITY LETTERS 
letters_grad$Letter[ letters_grad$variable == "Circularity" & letters_grad$SpeciesPct == "0%" ] <- "a" 
letters_grad$Letter[ letters_grad$variable == "Circularity" & letters_grad$SpeciesPct == "15%" ] <- "c" 
letters_grad$Letter[ letters_grad$variable == "Circularity" & letters_grad$SpeciesPct == "51%" ] <- "d" 
letters_grad$Letter[ letters_grad$variable == "Circularity" & letters_grad$SpeciesPct == "54%" ] <- "e" 
letters_grad$Letter[ letters_grad$variable == "Circularity" & letters_grad$SpeciesPct == "59%" ] <- "f" 
letters_grad$Letter[ letters_grad$variable == "Circularity" & letters_grad$SpeciesPct == "63%" ] <- "g" 
letters_grad$Letter[ letters_grad$variable == "Circularity" & letters_grad$SpeciesPct == "68%" ] <- "e" 
letters_grad$Letter[ letters_grad$variable == "Circularity" & letters_grad$SpeciesPct == "94%" ] <- "h" 
letters_grad$Letter[ letters_grad$variable == "Circularity" & letters_grad$SpeciesPct == "100%" ] <- "b" 

############## ASPECT RATIO LETTERS 
letters_grad$Letter[ letters_grad$variable == "Aspect_Ratio" & letters_grad$SpeciesPct == "0%" ] <- "a" 
letters_grad$Letter[ letters_grad$variable == "Aspect_Ratio" & letters_grad$SpeciesPct == "15%" ] <- "c" 
letters_grad$Letter[ letters_grad$variable == "Aspect_Ratio" & letters_grad$SpeciesPct == "51%" ] <- "de" 
letters_grad$Letter[ letters_grad$variable == "Aspect_Ratio" & letters_grad$SpeciesPct == "54%" ] <- "b" 
letters_grad$Letter[ letters_grad$variable == "Aspect_Ratio" & letters_grad$SpeciesPct == "59%" ] <- "f" 
letters_grad$Letter[ letters_grad$variable == "Aspect_Ratio" & letters_grad$SpeciesPct == "63%" ] <- "d" 
letters_grad$Letter[ letters_grad$variable == "Aspect_Ratio" & letters_grad$SpeciesPct == "68%" ] <- "e" 
letters_grad$Letter[ letters_grad$variable == "Aspect_Ratio" & letters_grad$SpeciesPct == "94%" ] <- "f" 
letters_grad$Letter[ letters_grad$variable == "Aspect_Ratio" & letters_grad$SpeciesPct == "100%" ] <- "b" 

###### POROSITY LETTERS 
letters_grad$Letter[ letters_grad$variable == "Porosity_Percent" & letters_grad$SpeciesPct == "0%" ] <- "a" 
letters_grad$Letter[ letters_grad$variable == "Porosity_Percent" & letters_grad$SpeciesPct == "15%" ] <- "c" 
letters_grad$Letter[ letters_grad$variable == "Porosity_Percent" & letters_grad$SpeciesPct == "51%" ] <- "d" 
letters_grad$Letter[ letters_grad$variable == "Porosity_Percent" & letters_grad$SpeciesPct == "54%" ] <- "ef" 
letters_grad$Letter[ letters_grad$variable == "Porosity_Percent" & letters_grad$SpeciesPct == "59%" ] <- "e" 
letters_grad$Letter[ letters_grad$variable == "Porosity_Percent" & letters_grad$SpeciesPct == "63%" ] <- "f" 
letters_grad$Letter[ letters_grad$variable == "Porosity_Percent" & letters_grad$SpeciesPct == "68%" ] <- "g" 
letters_grad$Letter[ letters_grad$variable == "Porosity_Percent" & letters_grad$SpeciesPct == "94%" ] <- "h" 
letters_grad$Letter[ letters_grad$variable == "Porosity_Percent" & letters_grad$SpeciesPct == "100%" ] <- "b" 

#####WOBBLE LETTERS 
letters_grad$Letter[ letters_grad$variable == "Wobbliness" & letters_grad$SpeciesPct == "0%" ] <- "a" 
letters_grad$Letter[ letters_grad$variable == "Wobbliness" & letters_grad$SpeciesPct == "15%" ] <- "c" 
letters_grad$Letter[ letters_grad$variable == "Wobbliness" & letters_grad$SpeciesPct == "51%" ] <- "d" 
letters_grad$Letter[ letters_grad$variable == "Wobbliness" & letters_grad$SpeciesPct == "54%" ] <- "e" 
letters_grad$Letter[ letters_grad$variable == "Wobbliness" & letters_grad$SpeciesPct == "59%" ] <- "f" 
letters_grad$Letter[ letters_grad$variable == "Wobbliness" & letters_grad$SpeciesPct == "63%" ] <- "g" 
letters_grad$Letter[ letters_grad$variable == "Wobbliness" & letters_grad$SpeciesPct == "68%" ] <- "e" 
letters_grad$Letter[ letters_grad$variable == "Wobbliness" & letters_grad$SpeciesPct == "94%" ] <- "h" 
letters_grad$Letter[ letters_grad$variable == "Wobbliness" & letters_grad$SpeciesPct == "100%" ] <- "b" 

######FRACTAL DIMENSION LETTERS 
letters_grad$Letter[ letters_grad$variable == "Fractal_Dimension" & letters_grad$SpeciesPct == "0%" ] <- "a" 
letters_grad$Letter[ letters_grad$variable == "Fractal_Dimension" & letters_grad$SpeciesPct == "15%" ] <- "a" 
letters_grad$Letter[ letters_grad$variable == "Fractal_Dimension" & letters_grad$SpeciesPct == "51%" ] <- "c" 
letters_grad$Letter[ letters_grad$variable == "Fractal_Dimension" & letters_grad$SpeciesPct == "54%" ] <- "c" 
letters_grad$Letter[ letters_grad$variable == "Fractal_Dimension" & letters_grad$SpeciesPct == "59%" ] <- "d" 
letters_grad$Letter[ letters_grad$variable == "Fractal_Dimension" & letters_grad$SpeciesPct == "63%" ] <- "d" 
letters_grad$Letter[ letters_grad$variable == "Fractal_Dimension" & letters_grad$SpeciesPct == "68%" ] <- "e" 
letters_grad$Letter[ letters_grad$variable == "Fractal_Dimension" & letters_grad$SpeciesPct == "94%" ] <- "c" 
letters_grad$Letter[ letters_grad$variable == "Fractal_Dimension" & letters_grad$SpeciesPct == "100%" ] <- "b" 


## ============================================================
## 5. UPDATED FUNCTION TO ADD LETTERS
## ============================================================
add_letters_above <- function(p, dsub, metric, x_offset = 0) {
  
  df <- letters_grad[letters_grad$variable == metric, ]
  if (nrow(df) == 0) return(p)
  
  df$x <- as.numeric(df$SpeciesPct) + x_offset
  
  p +
    annotate(
      "text",
      x = df$x,
      y = Inf,           # always at top of panel
      vjust = -0.6,      # raise above panel
      label = df$Letter,
      size = 5.5,
      fontface = "bold",
      color = "black"
    ) +
    coord_cartesian(clip = "off")
}




## ============================================================
## ============================================================
## 6. UPDATED VIOLIN PLOT FUNCTION (FIXED Y-AXIS LIMITING)
## ============================================================

vir_fun <- col_numeric(viridis(256, option="D"), domain=c(0,100))

create_violin2 <- function(dat_long, metric, y_label,
                           plot_tag="(a)", show_x_axis=TRUE,
                           ymax_limit = NULL, 
                           x_offset = 0) {
  
  dsub <- dat_long %>% dplyr::filter(variable == metric)
  
  # Calculate the actual data range for positioning letters
  data_max <- max(dsub$value, na.rm = TRUE)
  
  ref <- dsub %>%
    dplyr::filter(SpeciesPct %in% c("0%", "100%")) %>%
    dplyr::group_by(SpeciesPct) %>%
    dplyr::summarise(median_value = median(value), .groups="drop") %>%
    dplyr::mutate(RefColor = SpeciesPct)
  
  p <- ggplot(dsub, aes(x = SpeciesPct, y = value, fill = Percent_numeric)) +
    geom_violin(scale="width", trim=TRUE, color="black", alpha=0.85) +
    stat_summary(fun=median, geom="crossbar", width=0.4)
  
  if (nrow(ref) > 0) {
    p <- p + geom_hline(
      data = ref,
      aes(yintercept = median_value, color = RefColor),
      linetype = "dashed",
      linewidth = 0.9
    )
  }
  
  p <- p +
    scale_fill_viridis_c(option="D", name="Treatment",
                         breaks=c(0,100),
                         labels=c(species0_label, species100_label)) +
    scale_color_manual(name="Reference",
                       values=c("0%"=vir_fun(0),"100%"=vir_fun(100)))
  
  # --- SET Y-AXIS LIMITS (this clips violins properly) ---
  if (!is.null(ymax_limit)) {
    p <- p + scale_y_continuous(limits = c(0, ymax_limit), expand = c(0, 0))
    # Add letters based on the VISIBLE range (ymax_limit)
    letter_df <- letters_grad[letters_grad$variable == metric, ]
    if (nrow(letter_df) > 0) {
      letter_df$x <- as.numeric(letter_df$SpeciesPct) + x_offset
      letter_df$y <- ymax_limit  # Position at top of visible range
      
      p <- p + geom_text(
        data = letter_df,
        aes(x = x, y = y, label = Letter),
        inherit.aes = FALSE,
        vjust = -0.6,
        size = 5.5,
        fontface = "bold",
        color = "black"
      )
    }
  }
  
  p <- p +
    labs(x=NULL, y=y_label, tag=plot_tag) +
    coord_cartesian(clip = "off") +
    theme_bw() +
    theme(
      legend.position="right",
      axis.text.x = if(show_x_axis) element_text(angle=45,hjust=1,size=16,color="black")
      else element_blank(),
      axis.text.y = element_text(size=16,color="black"),
      text        = element_text(size=18, family="Arial"),
      panel.grid  = element_blank(),
      plot.margin = margin(t=35, r=10, b=10, l=10)
    )
  
  return(p)
}

## ============================================================
## 7. BUILD PANELS
## ============================================================
p19 <- create_violin2(dat2_long, "Equivalent_Spherical_Diameter_mm", "ESD (mm)", "(a)", FALSE, ymax_limit = 4 )
p20 <- create_violin2(dat2_long, "Perimeter_mm", "Perimeter (mm)", "(b)", FALSE, ymax_limit = 20)
p21 <- create_violin2(dat2_long, "Area_mm2", "Area (mm²)", "(c)", FALSE, ymax_limit = 3)
p22 <- create_violin2(dat2_long, "Circularity", "Circularity", "(d)", FALSE, ymax_limit = 1)
p23 <- create_violin2(dat2_long, "Aspect_Ratio", "Aspect Ratio", "(e)", FALSE, ymax_limit = 3)
p24 <- create_violin2(dat2_long, "Porosity_Percent", "2D-Porosity (%)", "(f)", FALSE, ymax_limit = 10)
p25 <- create_violin2(dat2_long, "Wobbliness", "Wobbliness", "(g)", TRUE, ymax_limit = 15)

## ============================================================
## 8. FRACTAL DIMENSION PANEL (h)
## ============================================================
create_fractal_panel <- function(dat,
                                 selected_species,
                                 speciespct_levels,
                                 species0_label,
                                 species100_label,
                                 ymax_limit = NULL,
                                 plot_tag = "(h)",
                                 y_label = "Fractal Dimension",
                                 x_offset = 0) {
  
  # ----- 1. COMPUTE FRACTAL DIMENSIONS -----
  fd_by_species <- dat %>%
    dplyr::filter(Species %in% selected_species) %>%
    dplyr::group_by(Species) %>%
    dplyr::summarise(
      Fractal_Dimension = {
        model <- lm(log10(Perimeter_mm) ~ log10(Area_mm2), na.action = na.omit)
        slope <- coef(model)[2]
        2 * slope
      },
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      SpeciesPct = case_when(
        Species == "Melosira sp."    ~ "0%",
        Species == "M85C15"          ~ "15%",
        Species == "M49C51"          ~ "51%",
        Species == "M46C54"          ~ "54%",
        Species == "M41C59"          ~ "59%",
        Species == "M37C63"          ~ "63%",
        Species == "M32C67"          ~ "68%",
        Species == "M6C94"           ~ "94%",
        Species == "Chaetoceros sp." ~ "100%",
        TRUE ~ NA_character_
      ),
      SpeciesPct      = factor(SpeciesPct, levels = speciespct_levels),
      Percent_numeric = as.numeric(gsub("%","", SpeciesPct)),
      variable        = "Fractal_Dimension",
      value           = Fractal_Dimension
    )
  
  # ----- 2. REFERENCE LINES -----
  ref_fd <- fd_by_species %>%
    dplyr::filter(SpeciesPct %in% c("0%", "100%")) %>%
    dplyr::mutate(Reference = SpeciesPct)
  
  # ----- 3. BASE PLOT -----
  p26 <- ggplot(fd_by_species,
                aes(x = SpeciesPct, y = value, fill = Percent_numeric)) +
    geom_point(shape = 21, size = 6, stroke = 0.4, color = "black") +
    geom_hline(data = ref_fd,
               aes(yintercept = value, color = SpeciesPct),
               linetype = "dashed", linewidth = 0.9,
               show.legend = FALSE) +
    scale_color_manual(values = c("0%" = vir_fun(0), "100%" = vir_fun(100))) +
    scale_fill_viridis_c(
      option = "D",
      name   = "Treatment",
      breaks = c(0, 100),
      labels = c(species0_label, species100_label)
    )
  
  # ----- 4. ADD LETTERS USING y = Inf (never clipped) -----
  letter_df <- letters_grad[letters_grad$variable == "Fractal_Dimension", ]
  
  if (nrow(letter_df) > 0) {
    letter_df$x <- as.numeric(letter_df$SpeciesPct) + x_offset
    
    p26 <- p26 +
      annotate(
        "text",
        x = letter_df$x,
        y = Inf,
        vjust = -0.6,
        label = letter_df$Letter,
        size = 5.5,
        fontface = "bold",
        color = "black"
      )
  }
  
  # ----- 5. Y-AXIS -----
  p26 <- p26 +
    scale_y_continuous(
      limits = c(0, ymax_limit),
      breaks = scales::pretty_breaks(),
      expand = c(0, 0)
    ) +
    coord_cartesian(clip = "off")
  
  # ----- 6. FINAL THEMING -----
  p26 +
    labs(x = NULL, y = y_label, tag = plot_tag) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1,
                                 size = 16, color = "black"),
      axis.text.y = element_text(size = 16, color = "black"),
      text        = element_text(size = 18, family = "Arial"),
      panel.grid  = element_blank(),
      plot.margin = margin(t = 35, r = 10, b = 10, l = 10)
    )
}


p26 <- create_fractal_panel(
  dat               = dat,
  selected_species  = selected_species,
  speciespct_levels = speciespct_levels,
  species0_label    = species0_label,
  species100_label  = species100_label,
  ymax_limit        = 2.2,
  plot_tag          = "(h)"
)

p26
## ============================================================
## 9. COMBINE PANELS
## ============================================================
final_fig <-
  (p19 + p20 + p21 + p22 + p23 + p24 + p25 + p26) +
  plot_layout(ncol = 2, guides = "collect") &
  scale_fill_viridis_c(
    option = "D",
    name   = "Treatment",
    limits = c(0,100),
    breaks = c(0,100),
    labels = c(species0_label, species100_label)
  )

#print(final_fig)

## ============================================================
## 10. SAVE FIGURE
## ============================================================
ggsave("Fig6_Gradient1_Violins.jpg", final_fig,
       width=14, height=12, dpi=600)

