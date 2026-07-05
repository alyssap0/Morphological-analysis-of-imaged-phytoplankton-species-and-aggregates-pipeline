# ============================================================
# 0. Load packages 
# ============================================================

library(ggplot2)
library(tidyr)
library(scales)
library(ggpubr)
library(patchwork)
library(FSA)
library(rstatix)
library(multcompView)
library(viridis)
library(dplyr)


## ============================================================
## 1. SETTINGS — ONLY Skel → Thal gradient in this script
## ============================================================
speciespct_levels <- c("0%", "1%", "12%", "18%", "36%", "61%", "76%", "100%")
species0_label   <- "Skeletonema sp."
species100_label <- "Thalassionema sp."

setwd("C:/Users/Cashies/OneDrive - University of Tasmania/AP_THESIS/Data")
all_data <- read.csv("agg_data.csv", header = TRUE)

## ============================================================
## 2. BASE DATA
## ============================================================
dat <- all_data %>%
  dplyr::select(
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
    Porosity_percent  = Porosity * 100,
    Wobble            = Perimeter_mm / Equivalent_Spherical_Diameter_mm,
    Fractal_Dimension = 2 * log10(Perimeter_mm) / log10(Area_mm2)
  ) %>%
  dplyr::filter(Porosity_percent > 0)

selected_species <- c(
  "Skeletonema sp.", "S99N1", "S88N12", "S82N18",
  "S64N36", "S49N61", "S24N76", "S1N100", "Thalassionema sp."
)

dat2 <- dat %>%
  dplyr::filter(Species %in% selected_species) %>%
  dplyr::mutate(Species = factor(Species, levels = selected_species))


## ============================================================
## 3. CONVERT TO LONG FORMAT
## ============================================================
convert_gradient <- function(dat2) {
  recode_map <- c(
    "Skeletonema sp."   = "0%",
    "S99N1"             = "1%",
    "S88N12"            = "12%",
    "S82N18"            = "18%",
    "S64N36"            = "36%",
    "S49N61"            = "61%",
    "S24N76"            = "76%",
    "S1N100"            = "100%",
    "Thalassionema sp." = "100%"
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
## 4. SIGNIFICANCE LETTERS
## ============================================================
metrics_violin <- c(
  "Equivalent_Spherical_Diameter_mm",
  "Perimeter_mm",
  "Area_mm2",
  "Circularity",
  "Aspect_Ratio",
  "Porosity_percent",
  "Wobble"
)

vars_for_letters <- c(metrics_violin, "Fractal_Dimension")

letters_list <- list()

for (var in vars_for_letters) {
  dsub <- dat2_long %>% dplyr::filter(variable == var)
  if (nrow(dsub) == 0) next
  
  kw <- kruskal_test(value ~ SpeciesPct, data = dsub)
  
  if (kw$p < 0.05) {
    dunn <- dunn_test(dsub, value ~ SpeciesPct, p.adjust.method = "BH")
    
    groups <- sort(unique(c(dunn$group1, dunn$group2)))
    mat <- matrix(1, length(groups), length(groups),
                  dimnames=list(groups, groups))
    
    for (i in seq_len(nrow(dunn))) {
      g1 <- dunn$group1[i]; g2 <- dunn$group2[i]
      mat[g1,g2] <- dunn$p.adj[i]
      mat[g2,g1] <- dunn$p.adj[i]
    }
    
    letters <- multcompLetters(mat, threshold = 0.05)$Letters
    
    letters_list[[var]] <- data.frame(
      SpeciesPct = factor(names(letters), levels = levs_speciespct),
      Letter     = unname(letters),
      variable   = var
    )
  }
}

letters_grad <- bind_rows(letters_list)

## ============================================================
## 5. UPDATED FUNCTION TO ADD LETTERS
## ============================================================
add_letters_above <- function(p, dsub, metric, ymax_limit, x_offset = 0) {
  
  df <- letters_grad[letters_grad$variable == metric, ]
  if (nrow(df) == 0) return(p)
  
  df$x <- as.numeric(df$SpeciesPct) + x_offset
  df$y <- ymax_limit * 1.16  # Position letters above the plot
  
  p + 
    geom_text(
      data = df,
      aes(x = x, y = y, label = Letter),
      inherit.aes = FALSE,
      size = 5.5,
      fontface = "bold",
      color = "black",
      vjust = 0
    ) +
    scale_y_continuous(
      limits = c(0, ymax_limit * 1.16),  # Expand plot area for letters
      breaks = scales::pretty_breaks(),
      expand = c(0, 0)
    ) +
    coord_cartesian(clip = "off")  # Allow letters to show in margin
}

## ============================================================
## 6. UPDATED VIOLIN PLOT FUNCTION
## ============================================================

vir_fun <- col_numeric(viridis(256, option="D"), domain=c(0,100))

create_violin2 <- function(dat_long, metric, y_label,
                           plot_tag="(a)", show_x_axis=TRUE,
                           ymax_limit = NULL, 
                           x_offset = 0) {
  
  dsub <- dat_long %>% dplyr::filter(variable == metric)
  
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
                       values=c("0%"=vir_fun(0),"100%"=vir_fun(100))) +
    labs(x=NULL, y=y_label, tag=plot_tag) +
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
  
  # --- ADD LETTERS ABOVE ---
  if (!is.null(ymax_limit)) {
    p <- add_letters_above(p, dsub, metric, ymax_limit, x_offset)
  }
  
  return(p)
}

## ============================================================
## 7. BUILD PANELS
## ============================================================
p11 <- create_violin2(dat2_long, "Equivalent_Spherical_Diameter_mm", "ESD (mm)", "(a)", FALSE, ymax_limit = 4)
p12 <- create_violin2(dat2_long, "Perimeter_mm", "Perimeter (mm)", "(b)", FALSE, ymax_limit = 20)
p13 <- create_violin2(dat2_long, "Area_mm2", "Area (mm²)", "(c)", FALSE, ymax_limit = 5)
p14 <- create_violin2(dat2_long, "Circularity", "Circularity", "(d)", FALSE, ymax_limit = 1)
p15 <- create_violin2(dat2_long, "Aspect_Ratio", "Aspect Ratio", "(e)", FALSE, ymax_limit = 3)
p16 <- create_violin2(dat2_long, "Porosity_percent", "2D-Porosity (%)", "(f)", FALSE, ymax_limit = 15)
p17 <- create_violin2(dat2_long, "Wobble", "Wobbliness", "(g)", TRUE, ymax_limit = 15)

## ============================================================
## 8. FRACTAL DIMENSION PANEL (h)
## ============================================================
create_fractal_panel <- function(dat,
                                 selected_species,
                                 speciespct_levels,
                                 species0_label,
                                 species100_label,
                                 ymax_limit = NULL,
                                 y_offset = 1.16,    # controls height of letters
                                 plot_tag = "(h)",
                                 y_label = "Fractal Dimension") {
  
  ## ------------------------------------------------------------
  ## 1. Compute FD for each species (slope of log10(P) ~ log10(A))
  ## ------------------------------------------------------------
  fd_by_species <- dat %>%
    dplyr::filter(Species %in% selected_species) %>%
    dplyr::mutate(
      Species = case_when(
        Species %in% c("S1N100", "Thalassionema sp.") ~ "100pct_mix",
        TRUE ~ as.character(Species)
      )
    ) %>%
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
        Species == "Skeletonema sp." ~ "0%",
        Species == "S99N1"           ~ "1%",
        Species == "S88N12"          ~ "12%",
        Species == "S82N18"          ~ "18%",
        Species == "S64N36"          ~ "36%",
        Species == "S49N61"          ~ "61%",
        Species == "S24N76"          ~ "76%",
        Species == "100pct_mix"      ~ "100%",
        TRUE ~ NA_character_
      ),
      SpeciesPct      = factor(SpeciesPct, levels = speciespct_levels),
      Percent_numeric = as.numeric(gsub("%","", SpeciesPct)),
      variable        = "Fractal_Dimension",
      value           = Fractal_Dimension
    )
  
  ## ------------------------------------------------------------
  ## 2. Pull letters for FD and reference lines
  ## ------------------------------------------------------------
  ref_fd <- fd_by_species %>%
    dplyr::filter(SpeciesPct %in% c("0%", "100%")) %>%
    dplyr::mutate(Reference = SpeciesPct)
  
  
  letter_df <- letters_grad[letters_grad$variable == "Fractal_Dimension", ]
  
  ## ------------------------------------------------------------
  ## 3. Compute height for letters
  ## ------------------------------------------------------------
  y_top <- ymax_limit * y_offset
  
  ## ------------------------------------------------------------
  ## 4. Build base plot
  ## ------------------------------------------------------------
  p18 <- ggplot(fd_by_species,
                aes(x = SpeciesPct, y = value, fill = Percent_numeric)) +
    geom_point(shape = 21, size = 6, stroke = 0.4, color = "black") +
    geom_hline(data = ref_fd, aes(yintercept = value, color = SpeciesPct), linetype = "dashed",
               linewidth = 0.9, show.legend = FALSE)+
    scale_color_manual(name="Reference",
                       values=c("0%"=vir_fun(0),"100%"=vir_fun(100))) +
    scale_fill_viridis_c(
      option = "D",
      name = "Treatment",
      breaks = c(0,100),
      labels = c(species0_label, species100_label)
    )
  
  ## ------------------------------------------------------------
  ## 5. Add letters directly aligned to each SpeciesPct
  ## ------------------------------------------------------------
  if (nrow(letter_df) > 0) {
    p18 <- p18 +
      geom_text(
        data = letter_df,
        aes(x = SpeciesPct, y = y_top, label = Letter),
        inherit.aes = FALSE,
        size = 5.5,
        fontface = "bold",
        color = "black",
        vjust = 0
      )
  }
  
  ## ------------------------------------------------------------
  ## 6. Apply y-limits with extra space for letters
  ## ------------------------------------------------------------
  p18 <- p18 + 
    scale_y_continuous(
      limits = c(0, ymax_limit * 1.16),
      breaks = scales::pretty_breaks(),
      expand = c(0, 0)
    ) +
    coord_cartesian(clip = "off")
  
  ## ------------------------------------------------------------
  ## 7. Theme & labels
  ## ------------------------------------------------------------
  p18 +
    labs(x = NULL, y = y_label, tag = plot_tag) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 16, color = "black"),
      axis.text.y = element_text(size = 16, color = "black"),
      text        = element_text(size = 18, family = "Arial"),
      panel.grid  = element_blank(),
      plot.margin = margin(t = 35, r = 10, b = 10, l = 10)
    )
}

## ------------------------------------------------------------
## CALL IT
## ------------------------------------------------------------
p18 <- create_fractal_panel(
  dat               = dat,
  selected_species  = selected_species,
  speciespct_levels = speciespct_levels,
  species0_label    = species0_label,
  species100_label  = species100_label,
  ymax_limit        = 2.2,     # your chosen range
  y_offset          = 1.16,    # raise/lower letters
  plot_tag          = "(h)"
)

## ============================================================
## 9. COMBINE PANELS + APPLY *ONE* GLOBAL FILL SCALE
## ============================================================
final_fig <-
  (p11 + p12 + p13 + p14 + p15 + p16 + p17 + p18) +
  plot_layout(ncol = 2, guides = "collect") &
  scale_fill_viridis_c(
    option = "D",
    name   = "Treatment",
    limits = c(0,100),
    breaks = c(0,100),
    labels = c(species0_label, species100_label)
  )


print(final_fig)


## ============================================================
## 10. SAVE FIGURE
## ============================================================

ggsave("Fig7_Gradient_Violins.jpg", final_fig,
       width=14, height=12, dpi=600)
