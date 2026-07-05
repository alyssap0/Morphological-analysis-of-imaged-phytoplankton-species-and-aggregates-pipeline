# ===================================================================
# DATA LOADING + REPLICATE TRAJECTORY PLOT + SPECIES PEAK-DAY BARPLOT
# (extracted from the full aggregate morphology LMM script)
# ===================================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(myplotfunction)

# ===================================================================
# 0. LOAD DATA, FILTER EXPERIMENTS, CREATE REPLICATES
# ===================================================================
setwd("C:/Users/aarfer/OneDrive - NOC/Documents/Alyssa_manuscript")
all_data <- read.csv("agg_data.csv", header = TRUE)

all_species <- c("Chaetoceros sp.", "Melosira sp.", "Skeletonema sp.",
                 "Ditylum sp.", "Thalassionema sp.")

# Canonical species-experiment filtering
all_species_data <- all_data %>%
  filter(
    (Species == "Chaetoceros sp."   & Experiment.No. == "Exp1") |
      (Species == "Thalassionema sp." & Experiment.No. == "Exp2") |
      (Species == "Melosira sp."      & Experiment.No. == "Exp3") |
      (Species == "Ditylum sp."       & Experiment.No. == "Exp1") |
      (Species == "Skeletonema sp."   & Experiment.No. == "Exp3")
  ) %>%
  filter(Species %in% all_species) %>%
  mutate(
    Species = factor(Species)
  )

# ===============================================================
# Create replicate IDs WITHIN each species
# ===============================================================
all_species_data <- all_species_data %>%
  mutate(
    Microcosm_num = as.numeric(
      gsub("M|_.*", "", Microcosm.No.)
    )
  ) %>%
  group_by(Species) %>%
  mutate(
    Replicate = paste0(
      "R",
      match(Microcosm.No., unique(Microcosm.No.))
    )
  ) %>%
  ungroup() %>%
  mutate(
    Replicate = factor(Replicate)
  ) %>%
  droplevels()


# ===================================================================
# REPLICATE TRAJECTORY PLOT
# Mean trait value per replicate per day, faceted by Species x Trait
# Uses all_species_data (all 5 species; only species with repeated
# days will show more than a single point per line)
# ===================================================================

# -------------------------------------------------------------------
# 1. Calculate per-Microcosm per-day means (+ SD/SE) on the RAW scale
#    (raw values are more biologically interpretable than log/logit
#    for a figure about replicate consistency)
# -------------------------------------------------------------------
traits_plot <- c("Area_mm2", "Circularity",
                 "Equivalent_Spherical_Diameter_mm", "Aspect_Ratio")

# Nice labels for facet strips
trait_labels <- c(
  Area_mm2                         = "Area (mm²)",
  Circularity                      = "Circularity",
  Equivalent_Spherical_Diameter_mm = "ESD (mm)",
  Aspect_Ratio                     = "Aspect Ratio"
)

microcosm_means <- all_species_data %>%
  group_by(Species, Replicate, Day.No.) %>%
  summarise(
    across(
      all_of(traits_plot),
      list(
        mean = ~ mean(.x, na.rm = TRUE),
        sd   = ~ sd(.x,   na.rm = TRUE),
        n    = ~ sum(!is.na(.x))
      ),
      .names = "{.col}__{.fn}"
    ),
    .groups = "drop"
  ) %>%
  pivot_longer(
    cols      = -c(Species, Replicate, Day.No.),
    names_to  = c("Trait", ".value"),
    names_sep = "__"
  ) %>%
  mutate(
    se = sd / sqrt(n)
  ) %>%
  rename(Mean_value = mean) %>%
  mutate(
    Microcosm = factor(
      case_when(
        Replicate == "R1" ~ "M1",
        Replicate == "R2" ~ "M2",
        Replicate == "R3" ~ "M3"
      ),
      levels = c("M1", "M2", "M3")
    ),
    Trait = factor(Trait,
                   levels = c("Area_mm2",
                              "Equivalent_Spherical_Diameter_mm",
                              "Aspect_Ratio",
                              "Circularity",
                              "Elongation"),
                   labels = c("Area (mm²)", "ESD (mm)",
                              "Aspect Ratio", "Circularity", "Elongation"))
  )

microcosm_means <- microcosm_means %>%
  mutate(
    Species = factor(
      Species,
      levels = c(
        "Chaetoceros sp.",
        "Melosira sp.",
        "Skeletonema sp.",
        "Ditylum sp.",
        "Thalassionema sp."
      )
    )
  )

# -------------------------------------------------------------------
# 2. Plot
# -------------------------------------------------------------------
# Colour palette: 3 microcosms, colourblind-friendly
microcosm_colours <- c(M1 = "#E69F00", M2 = "#56B4E9", M3 = "#009E73")

p_trajectories <- ggplot(
  microcosm_means,
  aes(x     = Day.No.,
      y     = Mean_value,
      colour = Microcosm,
      fill   = Microcosm,
      group  = Microcosm)
) +
  # --- Error ribbon: choose SD (active) or SE (commented out) ---
  #  geom_ribbon(
  #    aes(ymin = Mean_value - sd, ymax = Mean_value + sd),
  #    alpha = 0.15, colour = NA
  #  ) +
  geom_ribbon(
    aes(ymin = Mean_value - se, ymax = Mean_value + se),
    alpha = 0.15, colour = NA
  ) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.5) +
  scale_x_continuous(breaks = 1:4) +
  scale_colour_manual(values = microcosm_colours) +
  scale_fill_manual(values = microcosm_colours, guide = "none") +
  facet_grid(Trait ~ Species, scales = "free", switch = "y") +
  labs(
    x      = "Day",
    y      = NULL,
    colour = "Microcosm",
  ) +
  theme_custom(base_size = 11) +
  theme(
    strip.placement    = "outside",
    strip.background   = element_blank(),
    strip.text.x       = element_text(face = "italic", size = 11),
    strip.text.y       = element_text(face = "bold", size = 10, angle = 90),
    panel.spacing      = unit(0.8, "lines"),
    legend.position    = "bottom",
    legend.title       = element_text(face = "bold"),
    plot.title         = element_text(face = "bold", size = 12, hjust = 0.5),
    axis.title.x       = element_text(margin = margin(t = 8))
  )

print(p_trajectories)

ggsave(
  file.path("LMM_diagnostic_plots", "replicate_trajectories_ALL.png"),
  p_trajectories,
  width  = 9,
  height = 9,
  dpi    = 600
)

cat("\nReplicate trajectory plot saved to LMM_diagnostic_plots/replicate_trajectories.png\n")


# ===================================================================
# SPECIES MORPHOLOGY AT PEAK CUMULATIVE-AREA DAY (bar chart)
# ===================================================================
# For each microcosm, find the day with the LARGEST CUMULATIVE
# Area_mm2 (sum of all aggregate areas that day). Keep only that
# tank's own peak-day rows, then summarise trait means/SD per species.
# ===================================================================

# -------------------------------------------------------------------
# 1. Identify each replicate's peak cumulative-Area day
# -------------------------------------------------------------------
area_by_day <- all_species_data %>%
  group_by(Species, Replicate, Day.No.) %>%
  summarise(
    cumulative_Area_mm2 = sum(Area_mm2, na.rm = TRUE),
    n_aggregates = n(),
    .groups = "drop"
  )

peak_area_day <- area_by_day %>%
  group_by(Species, Replicate) %>%
  slice_max(cumulative_Area_mm2, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  rename(peak_Day = Day.No., peak_cumulative_Area_mm2 = cumulative_Area_mm2) %>%
  arrange(Species, Replicate)

# -------------------------------------------------------------------
# 2. Subset original data to just each tank's peak day
# -------------------------------------------------------------------
peak_area_data <- all_species_data %>%
  inner_join(
    peak_area_day %>% select(Species, Replicate, peak_Day),
    by = c("Species", "Replicate")
  ) %>%
  filter(Day.No. == peak_Day) %>%
  droplevels()

# -------------------------------------------------------------------
# 3. Summarise trait means/SD per species at their peak day
# -------------------------------------------------------------------
traits_plot <- c("Area_mm2", "Circularity",
                 "Equivalent_Spherical_Diameter_mm", "Length_mm")

trait_labels <- c(
  Area_mm2                         = "Area (mm²)",
  Circularity                      = "Circularity",
  Equivalent_Spherical_Diameter_mm = "ESD (mm)",
  Length_mm                        = "Length (mm)"
)

species_peak_summary <- peak_area_data %>%
  group_by(Species) %>%
  summarise(
    across(
      all_of(traits_plot),
      list(
        mean = ~ mean(.x, na.rm = TRUE),
        sd   = ~ sd(.x, na.rm = TRUE)
      ),
      .names = "{.col}_{.fn}"
    ),
    .groups = "drop"
  ) %>%
  pivot_longer(
    cols = -Species,
    names_to = c("Trait", "Statistic"),
    names_pattern = "(.*)_(mean|sd)",
    values_to = "Value"
  ) %>%
  pivot_wider(
    names_from = Statistic,
    values_from = Value
  ) %>%
  mutate(
    Trait = factor(
      Trait,
      levels = traits_plot,
      labels = trait_labels
    )
  )

# -------------------------------------------------------------------
# 4. Plot
# -------------------------------------------------------------------
p_species_peak_bars <- ggplot(
  species_peak_summary,
  aes(
    x = Species,
    y = mean
  )
) +
  geom_col(
    width = 0.7
  ) +
  geom_errorbar(
    aes(
      ymin = mean - sd,
      ymax = mean + sd
    ),
    width = 0.2,
    linewidth = 0.6
  ) +
  facet_wrap(
    ~Trait,
    scales = "free_y",
    ncol = 2
  ) +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme_custom(base_size = 11) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    panel.spacing = unit(1, "lines")
  )

print(p_species_peak_bars)

ggsave(
  file.path("LMM_diagnostic_plots", "species_peak_day_barplot.png"),
  p_species_peak_bars,
  width = 8,
  height = 6,
  dpi = 600
)