# ===================================================================
# AGGREGATE MORPHOLOGY LMMs: SPECIES x DAY
# ===================================================================
# Structure:
#   0. Load data, filter Species/Experiment.No. pairing ONCE, used
#      throughout the whole script
#   1. Shared helper functions (assumption checks, plots, variance
#      breakdown, significant-term extraction)
#   2. MODEL 1 - Species * Day (CONTINUOUS) + random slopes
#   3. MODEL 2 - Species * factor(Day) + random intercept
#   4. MODEL 3 - Species only, each microcosm's own peak-cumulative-
#      -Area day (all 5 species)
#   5. Notes on interpretation
#
# WHY THREE MODELS (see prior discussion in this analysis):
#   - Continuous Day assumes a straight-line trend over time. The
#     cumulative-Area-by-day plots showed this is false (non-linear,
#     non-monotonic, inconsistent direction across replicates -
#     clearest in Melosira). Kept here as a baseline / for comparison,
#     but NOT the primary model for testing whether Species matters.
#   - factor(Day.No.) makes no linearity assumption, comparing species
#     at each discrete calendar day - more honest about the shape of
#     change, but still compares species/tanks at the same CALENDAR
#     day, which is not necessarily the same DEVELOPMENTAL stage
#     (Day 1 = first day aggregates were detectable, which differs by
#     species/tank; some tanks also drop out before Day 4 once
#     aggregates sink, so later-day estimates rest on a non-random
#     subset of slower-developing replicates).
#   - The peak-cumulative-Area model sidesteps the calendar-day
#     mismatch by letting each tank define its own reference point
#     (the day it had the most aggregated biomass), at the cost of
#     comparing across days/tanks with very different absolute
#     magnitudes of aggregation (see prior discussion - Melosira R1's
#     peak is ~5-10x larger than R2/R3/R4's peaks).
#
# All three are kept and reported together: they test subtly
# different things, and the agreement/disagreement between them is
# itself part of the result (see Section 5).
# ===================================================================

library(lme4)
library(myplotfunction)
library(lmerTest)
library(dplyr)
library(tidyr)
library(MuMIn)        # r.squaredGLMM - fast R2, works with random slopes
library(ggplot2)
library(emmeans)      # estimated means/trends for result plots

# Keep emmeans' default behaviour of SKIPPING Kenward-Roger/Satterthwaite
# df adjustments above ~3000 rows (the warning you may see is
# informational, not an error). Recalculating these on large datasets
# would be extremely slow and isn't needed here, since emmeans is only
# used below for PLOTTING point estimates and asymptotic CIs, not for
# formal hypothesis tests (those come from anova() on the model itself
# via lmerTest). Do NOT raise lmerTest.limit/pbkrtest.limit to Inf.

# ===================================================================
# 0. LOAD DATA, FILTER EXPERIMENTS, CREATE REPLICATES
# ===================================================================
setwd("C:/Users/aarfer/OneDrive - NOC/Documents/Alyssa_manuscript")
all_data <- read.csv("agg_data.csv", header = TRUE)

all_species <- c("Chaetoceros sp.", "Melosira sp.", "Skeletonema sp.",
                 "Ditylum sp.", "Thalassionema sp.")

multi_day_species <- c("Chaetoceros sp.", "Melosira sp.", "Skeletonema sp.")


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
all_species_data <- all_species_data %>%
  mutate(
    Wobbliness = Perimeter_mm / Equivalent_Spherical_Diameter_mm
  )
range(all_species_data$Equivalent_Spherical_Diameter_mm)

all_species_data %>%
  group_by(Species) %>%
  summarise(
    min_ESD = min(Equivalent_Spherical_Diameter_mm, na.rm = TRUE),
    max_ESD = max(Equivalent_Spherical_Diameter_mm, na.rm = TRUE),
    range_ESD = max_ESD - min_ESD,
    .groups = "drop"
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


# Check mapping
cat("\nReplicate mapping:\n")
print(
  all_species_data %>%
    distinct(Species, Microcosm.No., Replicate) %>%
    arrange(Species, Replicate)
)


# Model data = only species with repeated days
model_data <- all_species_data %>%
  filter(Species %in% multi_day_species) %>%
  droplevels()


model_data$Day_c <- model_data$Day.No. - 1

traits <- c("Area_mm2", "Circularity",
            "Equivalent_Spherical_Diameter_mm", "Length_mm", "Aspect_Ratio", "Wobbliness")

# -------------------------------------------------------------
# 0a. Per-trait transform lookup
# -------------------------------------------------------------
# log() is appropriate for unbounded-positive, right-skewed traits,
# but Porosity is a proportion bounded in [0, 1]: log() compresses
# values near 1, explodes as the value -> 0, and is undefined
# (-Inf) at exactly 0. The correct transform for a bounded
# proportion is the logit, not the log - hence "logit" below.
#
# Circularity and Porosity are both bounded proportions in [0,1].
# For bounded responses, a logit transformation is generally more
# appropriate than a log transform because:
#
# - log() is undefined at 0
# - log() does not respect the upper boundary at 1
# - variance is typically non-constant near the boundaries
#
# Both traits therefore use a boundary-safe logit transform.
# -------------------------------------------------------------

trait_transform <- c(
  Area_mm2                          = "log",
  Circularity                       = "logit",
  Equivalent_Spherical_Diameter_mm  = "log",
  Aspect_Ratio                      = "log",
  Length_mm                         = "log",
  Wobbliness                        = "log"
)

# Boundary-safe logit: qlogis(0) = -Inf and qlogis(1) = +Inf, which
# would silently drop or break rows if Porosity ever hits an exact
# 0 or 1. The Smithson & Verkuilen (2006) correction nudges only
# values AT the boundary just inside (0,1); values already strictly
# inside are left untouched. eps shrinks as n grows, so the
# correction is as small as possible for your sample size.
safe_logit <- function(p, eps = NULL) {
  n <- sum(!is.na(p))
  if (is.null(eps)) eps <- 1 / (2 * n)
  p_adj <- pmin(pmax(p, eps), 1 - eps)
  qlogis(p_adj)
}


# ===================================================================
# 1. SHARED HELPER FUNCTIONS
# ===================================================================

# -------------------------------------------------------------
# 1a. Assumption checks
# -------------------------------------------------------------
# Linear mixed models assume: (1) linearity in the modelled predictor,
# (2) normal residuals, (3) homogeneity of variance, (4) normal
# random effects, (5) no problematic influential points. Traits like
# Porosity/Circularity are bounded 0-1 and several earlier model
# summaries showed scaled residuals with extreme max values (>100),
# so violations are likely for at least some traits - these checks
# quantify how badly.
# -------------------------------------------------------------

run_assumption_checks <- function(model, response, label) {
  
  cat("\n--- Assumption checks:", label, "-", response, "---\n")
  
  resid_fitted_df <- data.frame(fitted = fitted(model), resid = resid(model))
  
  # geom_hex + GAM smoother instead of geom_point + LOESS: with large
  # N, a raw scatter is unreadable and LOESS hangs (~O(n^2))
  p_resid_fitted <- ggplot(resid_fitted_df, aes(x = fitted, y = resid)) +
    geom_hex(bins = 60) +
    scale_fill_viridis_c(trans = "log10", name = "count") +
    geom_hline(yintercept = 0, colour = "red", linetype = "dashed") +
    geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"),
                colour = "blue", se = FALSE) +
    labs(title = paste0(label, ": Residuals vs Fitted (", response, ")"),
         x = "Fitted values", y = "Residuals") +
    theme_classic()
  
  # Subsample to 20,000 points for the Q-Q plot - visually identical
  # at this density, much faster to render
  qq_sample <- resid_fitted_df
  if (nrow(qq_sample) > 20000) {
    qq_sample <- qq_sample[sample(nrow(qq_sample), 20000), ]
  }
  
  p_qq <- ggplot(qq_sample, aes(sample = resid)) +
    stat_qq(alpha = 0.15, size = 0.5) +
    stat_qq_line(colour = "red") +
    labs(title = paste0(label, ": Q-Q plot of residuals (", response, ")")) +
    theme_classic()
  
  # Random-effects normality: limited power with few tank-level
  # groups - treat as a rough sanity check only
  re <- ranef(model)[[1]]
  re_df <- data.frame(value = re[, 1], term = "Intercept")
  if (ncol(re) > 1) {
    re_df <- bind_rows(re_df, data.frame(value = re[, 2], term = colnames(re)[2]))
  }
  
  p_re_qq <- ggplot(re_df, aes(sample = value)) +
    stat_qq() +
    stat_qq_line(colour = "red") +
    facet_wrap(~ term, scales = "free") +
    labs(title = paste0(label, ": Q-Q plot of random effects (", response, ")")) +
    theme_classic()
  
  resid_summary <- summary(resid(model))
  cat("Residual summary:\n")
  print(resid_summary)
  
  list(
    resid_fitted_plot = p_resid_fitted,
    qq_plot            = p_qq,
    re_qq_plot          = p_re_qq,
    resid_summary       = resid_summary
  )
}

# -------------------------------------------------------------
# 1b. Result plots
# -------------------------------------------------------------

standardise_ci_cols <- function(emm_df) {
  # emmeans names CI columns differently depending on whether it uses
  # t-based CIs (lower.CL/upper.CL) or asymptotic z-based CIs
  # (asymp.LCL/asymp.UCL, used automatically above 3000 rows)
  if (!"lower.CL" %in% names(emm_df)) {
    emm_df$lower.CL <- emm_df$asymp.LCL
    emm_df$upper.CL <- emm_df$asymp.UCL
  }
  emm_df
}

plot_continuous_trend <- function(model, response, y_label) {
  emm_df <- emmeans(model, ~ Species | Day_c,
                    at = list(Day_c = seq(0, 3, by = 1))) %>%
    as.data.frame() %>%
    standardise_ci_cols()
  
  ggplot(emm_df, aes(x = Day_c, y = emmean, colour = Species, fill = Species,
                     group = Species)) +
    geom_ribbon(aes(ymin = lower.CL, ymax = upper.CL), alpha = 0.15, colour = NA) +
    geom_line() +
    geom_point(size = 2) +
    labs(x = "Day (centred, Day 1 = 0)", y = y_label,
         title = paste0("Model-estimated trend (continuous Day): ", response)) +
    theme_classic()
}

plot_factor_trend <- function(model, response, y_label) {
  emm_df <- emmeans(model, ~ Species | Day.No.) %>%
    as.data.frame() %>%
    standardise_ci_cols()
  
  ggplot(emm_df, aes(x = Day.No., y = emmean, colour = Species, fill = Species,
                     group = Species)) +
    geom_ribbon(aes(ymin = lower.CL, ymax = upper.CL), alpha = 0.15, colour = NA) +
    geom_line() +
    geom_point(size = 2) +
    labs(x = "Day", y = y_label,
         title = paste0("Model-estimated trend (factor Day): ", response)) +
    theme_classic()
}

plot_species_means <- function(model, response, y_label, title_suffix) {
  emm_df <- emmeans(model, ~ Species) %>%
    as.data.frame() %>%
    standardise_ci_cols()
  
  ggplot(emm_df, aes(x = Species, y = emmean, colour = Species)) +
    geom_pointrange(aes(ymin = lower.CL, ymax = upper.CL), size = 0.8) +
    labs(x = NULL, y = y_label,
         title = paste0(title_suffix, ": ", response)) +
    theme_classic() +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 20, hjust = 1))
}

# -------------------------------------------------------------
# 1c. Variance breakdown ("how much does each variable explain")
# -------------------------------------------------------------
# R2m (fixed effects only) / R2c (fixed + random) from r.squaredGLMM,
# plus VarCorr variance components converted to proportions of total
# variance. Species stays a FIXED effect throughout (only 3-5 levels -
# too few for a stable random-effect variance estimate, rule of thumb
# 5-10+ levels), with replicate/tank structure captured via random
# effects nested in Species:Replicate, so tank-level
# non-independence is properly accounted for rather than treating
# every aggregate as an independent observation.
# -------------------------------------------------------------

get_variance_breakdown <- function(model) {
  r2 <- r.squaredGLMM(model)
  vc <- as.data.frame(VarCorr(model))
  vc$prop <- vc$vcov / sum(vc$vcov, na.rm = TRUE)
  list(r2 = r2, variance_components = vc)
}

# -------------------------------------------------------------
# 1d. Extract which Species:Day terms are significant
# -------------------------------------------------------------
# Pulls significant (p < 0.05) fixed-effect rows directly out of
# summary(model)$coefficients, so you don't have to manually scan
# console output to see which specific Species:Day combinations
# differ - prints a clean filtered table per model/trait.
# -------------------------------------------------------------

get_significant_terms <- function(model, alpha = 0.05) {
  coefs <- as.data.frame(summary(model)$coefficients)
  coefs$Term <- rownames(coefs)
  sig <- coefs[coefs[["Pr(>|t|)"]] < alpha, ]
  sig[, c("Term", "Estimate", "Std. Error", "df", "t value", "Pr(>|t|)")]
}


# ===================================================================
# 2. MODEL 1 - SPECIES x DAY (CONTINUOUS) + RANDOM SLOPES
# ===================================================================
# Species * Day_c + (Day_c | Species:Replicate). Random slopes
# let each tank have its own intercept and rate of change, properly
# accounting for tank-level non-independence (the alternative,
# random-intercept-only, was shown earlier in this analysis to
# produce severely anti-conservative p-values for Day-related terms -
# do not substitute that structure here).
#
# Interpretation: tests whether there is a SINGLE LINEAR slope that
# differs by species, on log(Area_mm2)-type traits where relevant.
# Given the cumulative-Area-by-day plots showed non-linear,
# non-monotonic, replicate-inconsistent trajectories (clearest for
# Melosira), a null or weak result here should NOT be read as "Day
# doesn't matter" - it may mean "Day's effect isn't a straight line",
# which this model structurally cannot detect even if present.
# ===================================================================

run_model_continuous <- function(response, data) {
  transform <- trait_transform[[response]]
  if (transform == "log") {
    
    formula_lhs <- paste0("log(", response, ")")
    axis_label  <- paste0("log(", response, ")")
    
  } else if (transform == "logit") {
    
    new_var <- paste0(response, "_logit")
    
    data[[new_var]] <- safe_logit(data[[response]])
    
    formula_lhs <- new_var
    axis_label  <- paste0("logit(", response, ")")
    
  } else {
    
    formula_lhs <- response
    axis_label  <- response
    
  }
  
  model <- lmer(
    as.formula(paste0(formula_lhs, " ~ Species * Day_c + (Day_c | Species:Replicate)")),
    data = data, REML = TRUE
  )
  
  checks      <- run_assumption_checks(model, response, label = "Model 1: Continuous Day")
  trend_plot  <- plot_continuous_trend(model, response, y_label = axis_label)
  var_break   <- get_variance_breakdown(model)
  sig_terms   <- get_significant_terms(model)
  
  list(
    model         = model,
    summary       = summary(model),
    anova         = as.data.frame(anova(model, type = 3)) %>%
      mutate(Response = response, Model = "ContinuousDay_RandomSlopes"),
    r2            = var_break$r2,
    var_comp      = var_break$variance_components %>% mutate(Response = response),
    sig_terms     = sig_terms,
    singular      = isSingular(model),
    checks        = checks,
    trend_plot    = trend_plot
  )
}

# RESULTS STORED IN: continuous_results
continuous_results <- lapply(traits, function(tr) {
  
  cat("\n==========================================\n")
  cat("MODEL 1 (Continuous Day) - Trait:", tr, "\n")
  cat("==========================================\n")
  
  res <- run_model_continuous(tr, model_data)
  
  print(res$summary)
  cat("\n--- Type III ANOVA ---\n")
  print(res$anova)
  cat("\n--- R2m (fixed only) / R2c (fixed + random) ---\n")
  print(res$r2)
  cat("\n--- Variance components (proportion of total variance) ---\n")
  print(res$var_comp %>% select(grp, var1, vcov, prop))
  cat("\n--- Significant fixed-effect terms (p < 0.05) ---\n")
  print(res$sig_terms)
  cat("\nSingular fit:", res$singular, "\n")
  
  print(res$checks$resid_fitted_plot)
  print(res$checks$qq_plot)
  print(res$checks$re_qq_plot)
  print(res$trend_plot)
  
  res
})
names(continuous_results) <- traits

# --- Export Model 1 results ---
continuous_anova <- bind_rows(lapply(traits, function(tr) continuous_results[[tr]]$anova))
write.csv(continuous_anova, "Model1_ContinuousDay_ANOVA.csv", row.names = FALSE)

continuous_r2 <- bind_rows(lapply(traits, function(tr) {
  data.frame(Response = tr, t(continuous_results[[tr]]$r2))
}))
write.csv(continuous_r2, "Model1_ContinuousDay_R2.csv", row.names = FALSE)

continuous_var_comp <- bind_rows(lapply(traits, function(tr) continuous_results[[tr]]$var_comp))
write.csv(continuous_var_comp, "Model1_ContinuousDay_VarianceComponents.csv", row.names = FALSE)

continuous_sig_terms <- bind_rows(lapply(traits, function(tr) {
  st <- continuous_results[[tr]]$sig_terms
  if (nrow(st) > 0) st$Response <- tr
  st
}))
write.csv(continuous_sig_terms, "Model1_ContinuousDay_SignificantTerms.csv", row.names = FALSE)

for (tr in traits) {
  r <- continuous_results[[tr]]
  ggsave(file.path("LMM_diagnostic_plots", paste0(tr, "_m1_resid_fitted.png")),
         r$checks$resid_fitted_plot, width = 7, height = 5, dpi = 200)
  ggsave(file.path("LMM_diagnostic_plots", paste0(tr, "_m1_qq.png")),
         r$checks$qq_plot, width = 6, height = 6, dpi = 200)
  ggsave(file.path("LMM_diagnostic_plots", paste0(tr, "_m1_re_qq.png")),
         r$checks$re_qq_plot, width = 7, height = 5, dpi = 200)
  ggsave(file.path("LMM_diagnostic_plots", paste0(tr, "_m1_trend.png")),
         r$trend_plot, width = 7, height = 5, dpi = 200)
}

cat("\nModel 1 (Continuous Day) complete - CSVs and plots saved.\n")


# ===================================================================
# 3. MODEL 2 - SPECIES x DAY (FACTOR) + RANDOM INTERCEPT
# ===================================================================
# Species * factor(Day.No.) + (1 | Species:Replicate)
#
# No linearity assumption - estimates a separate mean response
# for each day.
#
# Random effect:
# (1 | Species:Replicate)
#
# Each replicate microcosm receives its own intercept, accounting
# for repeated observations collected from the same tank through
# time.
#
# Unlike Model 1, there is no continuous Day slope being estimated,
# therefore no replicate-specific slope term is required.
#
# Interpretation:
# tests whether species differ at specific calendar days and whether
# those species differences change through time.
#
# CAVEAT: comparing at the same calendar day still does not guarantee
# comparing the same developmental stage across species/tanks (see
# Model 3 for the alternative that addresses this).
# ===================================================================

run_model_factor <- function(response, data) {
  
  transform <- trait_transform[[response]]
  
  if (transform == "log") {
    
    formula_lhs <- paste0("log(", response, ")")
    axis_label  <- paste0("log(", response, ")")
    
  } else if (transform == "logit") {
    
    new_var <- paste0(response, "_logit")
    
    data[[new_var]] <- safe_logit(data[[response]])
    
    formula_lhs <- new_var
    axis_label  <- paste0("logit(", response, ")")
    
  } else {
    
    formula_lhs <- response
    axis_label  <- response
    
  }
  
  model <- lmer(
    as.formula(paste0(formula_lhs,
                      " ~ Species * factor(Day.No.) + (1 | Species:Replicate)")),
    data = data, REML = TRUE
  )
  
  checks      <- run_assumption_checks(model, response, label = "Model 2: Factor Day")
  trend_plot  <- plot_factor_trend(model, response, y_label = axis_label)
  var_break   <- get_variance_breakdown(model)
  sig_terms   <- get_significant_terms(model)
  
  list(
    model         = model,
    summary       = summary(model),
    anova         = as.data.frame(anova(model, type = 3)) %>%
      mutate(Response = response, Model = "FactorDay_RandomIntercept"),
    r2            = var_break$r2,
    var_comp      = var_break$variance_components %>% mutate(Response = response),
    sig_terms     = sig_terms,
    singular      = isSingular(model),
    checks        = checks,
    trend_plot    = trend_plot
  )
}

# RESULTS STORED IN: factor_results
factor_results <- lapply(traits, function(tr) {
  
  cat("\n==========================================\n")
  cat("MODEL 2 (Factor Day) - Trait:", tr, "\n")
  cat("==========================================\n")
  
  res <- run_model_factor(tr, model_data)
  
  print(res$summary)
  cat("\n--- Type III ANOVA ---\n")
  print(res$anova)
  cat("\n--- R2m (fixed only) / R2c (fixed + random) ---\n")
  print(res$r2)
  cat("\n--- Variance components (proportion of total variance) ---\n")
  print(res$var_comp %>% select(grp, var1, vcov, prop))
  cat("\n--- Significant fixed-effect terms (p < 0.05) ---\n")
  print(res$sig_terms)
  cat("\nSingular fit:", res$singular, "\n")
  
  print(res$checks$resid_fitted_plot)
  print(res$checks$qq_plot)
  print(res$checks$re_qq_plot)
  print(res$trend_plot)
  
  res
})
names(factor_results) <- traits

# --- Export Model 2 results ---
factor_anova <- bind_rows(lapply(traits, function(tr) factor_results[[tr]]$anova))
write.csv(factor_anova, "Model2_FactorDay_ANOVA.csv", row.names = FALSE)

factor_r2 <- bind_rows(lapply(traits, function(tr) {
  data.frame(Response = tr, t(factor_results[[tr]]$r2))
}))
write.csv(factor_r2, "Model2_FactorDay_R2.csv", row.names = FALSE)

factor_var_comp <- bind_rows(lapply(traits, function(tr) factor_results[[tr]]$var_comp))
write.csv(factor_var_comp, "Model2_FactorDay_VarianceComponents.csv", row.names = FALSE)

factor_sig_terms <- bind_rows(lapply(traits, function(tr) {
  st <- factor_results[[tr]]$sig_terms
  if (nrow(st) > 0) st$Response <- tr
  st
}))
write.csv(factor_sig_terms, "Model2_FactorDay_SignificantTerms.csv", row.names = FALSE)

for (tr in traits) {
  r <- factor_results[[tr]]
  ggsave(file.path("LMM_diagnostic_plots", paste0(tr, "_m2_resid_fitted.png")),
         r$checks$resid_fitted_plot, width = 7, height = 5, dpi = 200)
  ggsave(file.path("LMM_diagnostic_plots", paste0(tr, "_m2_qq.png")),
         r$checks$qq_plot, width = 6, height = 6, dpi = 200)
  ggsave(file.path("LMM_diagnostic_plots", paste0(tr, "_m2_re_qq.png")),
         r$checks$re_qq_plot, width = 7, height = 5, dpi = 200)
  ggsave(file.path("LMM_diagnostic_plots", paste0(tr, "_m2_trend.png")),
         r$trend_plot, width = 7, height = 5, dpi = 200)
}

cat("\nModel 2 (Factor Day) complete - CSVs and plots saved.\n")


# ===================================================================
# 4. MODEL 3 - SPECIES ONLY, PEAK-CUMULATIVE-AREA DAY (ALL 5 SPECIES)
# ===================================================================
# For each microcosm, finds the day with the LARGEST CUMULATIVE
# Area_mm2 (sum of all aggregate areas that day - reflects total
# imaged biomass-like signal, driven by both individual aggregate
# size AND particle count). Keeps only that tank's own peak-day rows,
# then tests Species alone. Uses explore_data (not model_data) so ALL
# 5 species are included - Ditylum/Thalassionema only have Day 1, so
# their "peak day" is trivially Day 1, which is still a valid
# single-day snapshot for this model (it doesn't need multi-day data).
#
# CAVEAT: peak cumulative Area is NOT on a comparable absolute scale
# across species/tanks (e.g. Melosira's largest replicate peak is
# roughly an order of magnitude smaller than its other replicates,
# and far smaller than Chaetoceros's typical values - see the
# cumulative-Area-by-day plots). This model aligns tanks on
# DEVELOPMENTAL STAGE (their own point of maximum aggregation) rather
# than calendar day, fixing the Model 1/2 alignment problem, but at
# the cost of comparing across very different absolute magnitudes of
# aggregation. Treat agreement/disagreement with Models 1-2 as part
# of the overall result, not as one model being simply "correct".
#
# Uses all_species_data (NOT raw explore_data) - this applies the same
# canonical Species:Experiment.No. filter as model_data, just without
# restricting to multi_day_species, so all 5 species are included
# without risking contamination from the wrong Experiment.No. for any
# species that appears under more than one experiment in the raw CSV.
# ===================================================================

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

print(peak_area_day)
#write.csv(peak_area_day, "Peak_day_cumulative_Area.csv", row.names = FALSE)

peak_area_plot <- ggplot(
  area_by_day,
  aes(x = Day.No., y = cumulative_Area_mm2, colour = Replicate,
      group = interaction(Species, Replicate))
) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2) +
  geom_point(
    data = peak_area_day,
    aes(x = peak_Day, y = peak_cumulative_Area_mm2),
    shape = 8, size = 3, colour = "black", inherit.aes = FALSE
  ) +
  facet_wrap(~ Species, scales = "free_y") +
  labs(x = "Day", y = expression("Cumulative Area (mm"^2*")"),
       title = "Cumulative aggregate area through time, peak day marked (*)") +
  theme_classic()

print(peak_area_plot)
#ggsave("LMM_diagnostic_plots/peak_day_cumulative_Area_trend.png", peak_area_plot,
#       width = 8, height = 6, dpi = 200)

peak_area_data <- all_species_data %>%
  inner_join(
    peak_area_day %>% select(Species, Replicate, peak_Day),
    by = c("Species", "Replicate")
  ) %>%
  filter(Day.No. == peak_Day) %>%
  droplevels()

cat("\nSpecies present in peak_area_data (should be all 5):\n")
print(table(peak_area_data$Species))

cat("\nRows contributing to peak-cumulative-Area-day model, per species/microcosm:\n")
print(peak_area_data %>% count(Species, Replicate, Day.No.))


run_model_peak_area <- function(response, data) {
  
  transform <- trait_transform[[response]]
  if (transform == "log") {
    
    formula_lhs <- paste0("log(", response, ")")
    axis_label  <- paste0("log(", response, ")")
    
  } else if (transform == "logit") {
    
    new_var <- paste0(response, "_logit")
    
    data[[new_var]] <- safe_logit(data[[response]])
    
    formula_lhs <- new_var
    axis_label  <- paste0("logit(", response, ")")
    
  } else {
    
    formula_lhs <- response
    axis_label  <- response
    
  }
  
  model <- lmer(
    as.formula(
      paste0(formula_lhs, " ~ Species + (1 | Species:Replicate)")
    ),
    data = data,
    REML = TRUE
  )
  n_tanks <- data %>% distinct(Species, Replicate) %>% nrow()
  
  checks      <- run_assumption_checks(model, response, label = "Model 3: Peak cumulative-Area day")
  means_plot  <- plot_species_means(
    model,
    response,
    y_label = axis_label,
    title_suffix = "Peak-cumulative-Area-day estimated means"
  )
  var_break   <- get_variance_breakdown(model)
  sig_terms   <- get_significant_terms(model)
  
  list(
    model      = model,
    summary    = summary(model),
    anova      = as.data.frame(anova(model, type = 3)) %>%
      mutate(Response = response, Model = "PeakCumulativeAreaDay_SpeciesOnly"),
    r2         = var_break$r2,
    var_comp   = var_break$variance_components %>% mutate(Response = response),
    sig_terms  = sig_terms,
    singular   = isSingular(model),
    n_tanks    = n_tanks,
    checks     = checks,
    means_plot = means_plot
  )
}

# RESULTS STORED IN: peak_area_results
peak_area_results <- lapply(traits, function(tr) {
  
  cat("\n==========================================\n")
  cat("MODEL 3 (Peak cumulative-Area day) - Trait:", tr, "\n")
  cat("==========================================\n")
  
  res <- run_model_peak_area(tr, peak_area_data)
  
  print(res$summary)
  cat("\n--- Type III ANOVA ---\n")
  print(res$anova)
  cat("\n--- R2m (fixed only) / R2c (fixed + random) ---\n")
  print(res$r2)
  cat("\n--- Variance components (proportion of total variance) ---\n")
  print(res$var_comp %>% select(grp, var1, vcov, prop))
  cat("\n--- Significant fixed-effect terms (p < 0.05) ---\n")
  print(res$sig_terms)
  cat("\nSingular fit:", res$singular, "| Tanks contributing:", res$n_tanks, "\n")
  
  print(res$checks$resid_fitted_plot)
  print(res$checks$qq_plot)
  print(res$means_plot)
  
  res
})
names(peak_area_results) <- traits

# --- Export Model 3 results ---
peak_area_anova <- bind_rows(lapply(traits, function(tr) peak_area_results[[tr]]$anova))
write.csv(peak_area_anova, "Model3_PeakCumulativeAreaDay_ANOVA.csv", row.names = FALSE)

peak_area_r2 <- bind_rows(lapply(traits, function(tr) {
  data.frame(Response = tr, t(peak_area_results[[tr]]$r2))
}))
write.csv(peak_area_r2, "Model3_PeakCumulativeAreaDay_R2.csv", row.names = FALSE)

peak_area_var_comp <- bind_rows(lapply(traits, function(tr) peak_area_results[[tr]]$var_comp))
write.csv(peak_area_var_comp, "Model3_PeakCumulativeAreaDay_VarianceComponents.csv", row.names = FALSE)

peak_area_sig_terms <- bind_rows(lapply(traits, function(tr) {
  st <- peak_area_results[[tr]]$sig_terms
  if (nrow(st) > 0) st$Response <- tr
  st
}))
write.csv(peak_area_sig_terms, "Model3_PeakCumulativeAreaDay_SignificantTerms.csv", row.names = FALSE)

diagnostics_table <- bind_rows(lapply(traits, function(tr) {
  data.frame(
    Response                = tr,
    Model1_singular          = continuous_results[[tr]]$singular,
    Model2_singular          = factor_results[[tr]]$singular,
    Model3_singular          = peak_area_results[[tr]]$singular,
    Model3_n_tanks           = peak_area_results[[tr]]$n_tanks
  )
}))
print(diagnostics_table)
write.csv(diagnostics_table, "All_models_diagnostics.csv", row.names = FALSE)

for (tr in traits) {
  r <- peak_area_results[[tr]]
  ggsave(file.path("LMM_diagnostic_plots", paste0(tr, "_m3_resid_fitted.png")),
         r$checks$resid_fitted_plot, width = 7, height = 5, dpi = 200)
  ggsave(file.path("LMM_diagnostic_plots", paste0(tr, "_m3_qq.png")),
         r$checks$qq_plot, width = 6, height = 6, dpi = 200)
  ggsave(file.path("LMM_diagnostic_plots", paste0(tr, "_m3_means.png")),
         r$means_plot, width = 5, height = 5, dpi = 200)
}

cat("\nModel 3 (Peak cumulative-Area day) complete - CSVs and plots saved.\n")
cat("\nAll three models complete. See LMM_diagnostic_plots/ for figures.\n")


# ===================================================================
# 5. NOTES ON INTERPRETING AND COMPARING THE THREE MODELS
# ===================================================================
# - Model 1 (continuous Day): tests for a single LINEAR trend per
#   species. A null Day or Species:Day result here means "no straight-
#   line trend detected" - given the non-linear, non-monotonic,
#   replicate-inconsistent trajectories visible in the cumulative-Area
#   plots, this should NOT be read as "Day doesn't matter".
#
# - Model 2 (factor Day): no linearity assumption - can detect
#   non-monotonic Species:Day patterns Model 1 cannot. Check
#   Model2_FactorDay_SignificantTerms.csv for which specific
#   Species:Day(level) combinations differ. Still compares
#   species/tanks at the same CALENDAR day, not necessarily the same
#   developmental stage.
#
# - Model 3 (peak cumulative-Area day): tests Species alone, aligned
#   on each tank's own developmental peak rather than calendar day -
#   addresses the calendar-day mismatch, but compares across very
#   different absolute magnitudes of peak aggregation between species/
#   replicates (see peak_day_cumulative_Area_trend.png).
#
# Circularity and Porosity are both bounded [0,1] traits and are
# therefore analysed using a boundary-safe logit transformation.
#
# - If Species is significant in ALL THREE models, that is reasonably
#   strong, robust evidence that species predicts the trait
#   REGARDLESS of how Day is handled.
#
# - If the Species:Day interaction differs in significance between
#   Model 1 and Model 2, that itself indicates the species difference
#   in TRAJECTORY (not mean level) is sensitive to the linearity
#   assumption - report this explicitly rather than picking whichever
#   model gives the "cleaner" result.
#
# - Residuals vs Fitted: look for a roughly flat, even band around
#   zero. A funnel shape indicates heteroscedasticity - common for
#   Area/Porosity. A curve indicates a non-linear relationship Model 1
#   specifically cannot capture.
#
# - Q-Q plot of residuals: expect heavy right-skew/long tails for
#   Area and Porosity given prior diagnostics. Area_mm2,
#   Equivalent_Spherical_Diameter_mm, and Length_mm are log-
#   transformed (trait_transform). Porosity is logit-transformed
#   (bounded [0,1] - log() is not appropriate for a proportion).
#   Circularity is also bounded [0,1] and is currently still on logit
#   
# - Q-Q plot of random effects: limited power with few tank-level
#   groups - treat as a rough sanity check only, not a strict test.
# ===================================================================


library(lmerTest)

m_full <- lmer(log(Area_mm2) ~ Species * factor(Day.No.) + (1 | Species:Replicate),
               data = model_data, REML = FALSE)

ranova(m_full)


m_full <- lmer(logit(Circularity) ~ Species * factor(Day.No.) + (1 | Species:Replicate),
               data = model_data, REML = FALSE)

ranova(m_full)


m_full <- lmer(log(Length_mm) ~ Species * factor(Day.No.) + (1 | Species:Replicate),
               data = model_data, REML = FALSE)

ranova(m_full)


m_full <- lmer(log(Equivalent_Spherical_Diameter_mm) ~ Species * factor(Day.No.) + (1 | Species:Replicate),
               data = model_data, REML = FALSE)

ranova(m_full)











# ===================================================================
# REPLICATE TRAJECTORY PLOT
# Mean trait value per replicate per day, faceted by Species x Trait
# Excludes Porosity. Uses model_data (3 multi-day species only).
# ===================================================================
library(dplyr)
library(tidyr)
library(ggplot2)
library(myplotfunction)
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





# -------------------------------------------------------------------
# Species morphology at peak cumulative-Area day (Model 3)
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
# Plot
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







