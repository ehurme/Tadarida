# =============================================================================
# Script 2: Displacement Threshold and Behavioral Classification
# Tadarida teniotis Fall Migration
#
# Adapted from CommonNoctuleSpringMigration/2_burst_flying_threshold.R
# and 5_behavior_thresholds.R
#
# KEY DIFFERENCE FROM NOCTULE PIPELINE:
# Noctule tags carry accelerometers (VeDBA) + GPS → flying vs. roosting
# distinguished at fix level. Sigfox/Nanofox tags on T. teniotis have no
# accelerometer. Classification must be based entirely on:
#   (1) Sigfox location accuracy  (sigfox_computed_location_radius)
#   (2) Daily net displacement     (straight-line distance day-to-day)
#   (3) Daily mean ground speed    (from move2 step metrics if available)
#   (4) Net squared displacement   (NSD — trajectory-level migratory signal)
#
# OUTPUTS:
#   displacement_threshold.robj   — filtered bats_daily with behavior labels
#   Plots/displacement_threshold_gmm.png
#   Plots/nsd_by_individual.png
#   Plots/behavior_summary.png
# =============================================================================

library(tidyverse)
library(move2)
library(sf)
library(mclust)      # GMM classification (same as noctule pipeline script 5)
library(scales)
library(patchwork)
library(geosphere)   # distGeo for great-circle distances

# Load data from Script 1
load("../../../Dropbox/MPI/Tadarida/Data/rdata/move_icarus_tt_daytime_daily.robj")


# =============================================================================
# SECTION 1: Location accuracy filter
# =============================================================================
# The sigfox_computed_location_radius is the TDOA multilateration error
# estimate in metres. Positions with very large radii are unreliable and
# should be excluded before calculating any movement metrics.
# Inspect the distribution first; do NOT hard-code a cutoff blindly.

# --- 1.1 Inspect radius distribution ---
p_radius <- bats_full %>%
  as_tibble() %>%
  filter(!is.na(sigfox_computed_location_radius)) %>%
  ggplot(aes(sigfox_computed_location_radius / 1000)) +
  geom_histogram(bins = 80, fill = "steelblue", col = "white", alpha = 0.8) +
  geom_vline(xintercept = 30, linetype = "dashed", col = "firebrick", linewidth = 0.8) +
  scale_x_log10(labels = scales::comma) +
  labs(
    title = "Sigfox location radius distribution",
    subtitle = "Dashed line = 30 km candidate threshold",
    x = "Location radius (km, log scale)",
    y = "Count"
  ) +
  theme_minimal()

# Cumulative: what fraction of fixes are retained at candidate thresholds?
radius_summary <- bats_full %>%
  as_tibble() %>%
  filter(!is.na(sigfox_computed_location_radius)) %>%
  summarise(
    n_total      = n(),
    pct_lt_10km  = mean(sigfox_computed_location_radius < 10000) * 100,
    pct_lt_20km  = mean(sigfox_computed_location_radius < 20000) * 100,
    pct_lt_30km  = mean(sigfox_computed_location_radius < 30000) * 100,
    pct_lt_50km  = mean(sigfox_computed_location_radius < 50000) * 100
  )
print(radius_summary)

# --- 1.2 Apply threshold ---
# Default 30 km matches noctule Sigfox studies; adjust after inspecting above.
RADIUS_THRESHOLD_M <- 30000

bats_loc_filt <- bats_loc %>%
  filter(sigfox_computed_location_radius <= RADIUS_THRESHOLD_M |
           is.na(sigfox_computed_location_radius))

cat(sprintf(
  "Location filter: retained %d / %d fixes (%.1f%%)\n",
  nrow(bats_loc_filt), nrow(bats_loc),
  100 * nrow(bats_loc_filt) / nrow(bats_loc)
))


# =============================================================================
# SECTION 2: Compute daily displacement metrics
# =============================================================================
# bats_daily already has one row per individual per day from import_nanofox_movebank.
# We need to add:
#   - net daily displacement  (km between consecutive daily positions)
#   - mean ground speed       (km/h, if step metrics exist in bats_daily)
#   - cumulative NSD          (squared displacement from deployment location)

# --- 2.1 Extract coordinates and sort ---
# Coordinates may be in an sf geometry column or explicit lon/lat columns
# depending on SigfoxTagPrep version. Handle both cases.

if (!"lon" %in% names(bats_daily)) {
  coords <- sf::st_coordinates(bats_daily)
  bats_daily <- bats_daily %>%
    mutate(
      lon = coords[, 1],
      lat = coords[, 2]
    )
}

bats_daily <- bats_daily %>%
  arrange(individual_local_identifier, timestamp)

# --- 2.2 Daily net displacement (great-circle, km) ---
# Calculated as step distance from the previous day's position.
# First fix per individual gets NA (no "previous" position).

bats_daily <- bats_daily %>%
  group_by(individual_local_identifier) %>%
  mutate(
    lon_lag = lag(lon),
    lat_lag = lag(lat),
    daily_displacement_km = case_when(
      is.na(lon_lag) ~ NA_real_,
      TRUE ~ geosphere::distGeo(
        cbind(lon_lag, lat_lag),
        cbind(lon,     lat)
      ) / 1000
    )
  ) %>%
  ungroup()

# --- 2.3 Net squared displacement from deployment location (NSD) ---
# NSD is a classical migration metric (Bunnefeld et al. 2011).
# We compute it as (great-circle distance from deployment location)^2 in km^2.

# Get deployment coordinates per individual from track_data
deploy_coords <- mt_track_data(bats_full) %>%
  as_tibble() %>%
  select(individual_local_identifier, deploy_on_location) %>%
  mutate(
    deploy_lon = sf::st_coordinates(deploy_on_location)[, 1],
    deploy_lat = sf::st_coordinates(deploy_on_location)[, 2]
  ) %>%
  select(individual_local_identifier, deploy_lon, deploy_lat)

bats_daily <- bats_daily %>%
  left_join(deploy_coords, by = "individual_local_identifier") %>%
  mutate(
    nsd_km2 = (geosphere::distGeo(
      cbind(deploy_lon, deploy_lat),
      cbind(lon, lat)
    ) / 1000)^2
  )

# --- 2.4 Inspect displacement distribution before classification ---
p_disp_hist <- bats_daily %>%
  filter(!is.na(daily_displacement_km)) %>%
  ggplot(aes(daily_displacement_km)) +
  geom_histogram(bins = 80, fill = "steelblue", col = "white", alpha = 0.8) +
  scale_x_log10(labels = scales::comma, breaks = c(1, 5, 10, 50, 100, 300, 500)) +
  labs(
    title = "Daily displacement distribution (all individuals)",
    x = "Daily displacement (km, log scale)",
    y = "Count (individual-days)"
  ) +
  theme_minimal()

p_disp_ind <- bats_daily %>%
  filter(!is.na(daily_displacement_km), daily_displacement_km > 0) %>%
  ggplot(aes(
    x     = reorder(individual_local_identifier, daily_displacement_km, median),
    y     = daily_displacement_km
  )) +
  geom_boxplot(outlier.size = 0.5, fill = "steelblue", alpha = 0.6) +
  coord_flip() +
  scale_y_log10() +
  labs(
    title = "Daily displacement by individual",
    x     = NULL,
    y     = "Daily displacement (km, log scale)"
  ) +
  theme_minimal(base_size = 9)

print(p_disp_hist / p_disp_ind)


# =============================================================================
# SECTION 3: GMM-based behavioral classification
# =============================================================================
# The noctule script 5 fits a Gaussian Mixture Model (via mclust) to VeDBA
# to find the stationary/flying threshold. Here we substitute log10 daily
# displacement as the primary classification variable, which is ecologically
# analogous: stationary days cluster at low displacement, migration nights
# at high displacement.
#
# We use a 2-component GMM (stationary vs. migrating) as default.
# A 3-component model (stationary / local movement / migrating) is also
# fitted and compared by BIC — choose based on your data and ecology.

# --- 3.1 Prepare data for GMM ---
# Use only days with valid displacement > 0 (0 km = same roost both days,
# unambiguously stationary; can be assigned directly without GMM)

gmm_data <- bats_daily %>%
  filter(!is.na(daily_displacement_km), daily_displacement_km > 0) %>%
  mutate(log10_disp = log10(daily_displacement_km))

x_gmm <- gmm_data$log10_disp

# --- 3.2 Fit 1-, 2-, 3-component GMMs and compare BIC ---
set.seed(42)
gmm_fits <- lapply(1:3, function(g) mclust::Mclust(x_gmm, G = g, modelNames = "V"))
bic_vals  <- sapply(gmm_fits, function(m) m$BIC)
best_G    <- which.max(bic_vals)

cat(sprintf(
  "GMM BIC: G=1: %.1f | G=2: %.1f | G=3: %.1f → best G=%d\n",
  bic_vals[1], bic_vals[2], bic_vals[3], best_G
))

# --- 3.3 Use 2-component model (stationary vs. migrating) as primary ---
# Even if G=3 wins by BIC, 2-component is more interpretable for
# binary state models downstream. Revisit if you need a "local foraging"
# state between stopover and long-distance migration.

gmm2 <- gmm_fits[[2]]

# Component means (in log10 km) — lower = stationary, higher = migratory
comp_means <- gmm2$parameters$mean
comp_order <- order(comp_means)   # 1 = lower mean (stationary), 2 = higher (migratory)

cat(sprintf(
  "GMM 2-component means: stationary = %.2f km | migratory = %.2f km\n",
  10^comp_means[comp_order[1]],
  10^comp_means[comp_order[2]]
))

# Posterior probabilities
gmm_data <- gmm_data %>%
  mutate(
    gmm_class        = gmm2$classification,
    behavior_gmm     = ifelse(gmm_class == comp_order[1], "stationary", "migrating"),
    p_migrating      = gmm2$z[, comp_order[2]],   # posterior prob of migratory class
    p_stationary     = gmm2$z[, comp_order[1]]
  )

# --- 3.4 Derive the threshold as the crossover point of the two densities ---
# This is the displacement value at which P(migrating|x) = P(stationary|x) = 0.5
# Useful as a reproducible, interpretable single cutoff for reporting.

x_seq <- seq(min(x_gmm), max(x_gmm), length.out = 2000)

dens_stat <- dnorm(x_seq,
                   mean = gmm2$parameters$mean[comp_order[1]],
                   sd   = sqrt(gmm2$parameters$variance$sigmasq[comp_order[1]])) *
  gmm2$parameters$pro[comp_order[1]]

dens_migr <- dnorm(x_seq,
                   mean = gmm2$parameters$mean[comp_order[2]],
                   sd   = sqrt(gmm2$parameters$variance$sigmasq[comp_order[2]])) *
  gmm2$parameters$pro[comp_order[2]]

# Crossover index
crossover_idx  <- which.min(abs(dens_stat - dens_migr))
DISP_THRESHOLD_KM <- 10^x_seq[crossover_idx]

cat(sprintf("Displacement threshold (GMM crossover): %.1f km/day\n", DISP_THRESHOLD_KM))


# =============================================================================
# SECTION 4: Visualize the GMM threshold
# =============================================================================

gmm_plot_df <- tibble(
  x_log10 = x_seq,
  d_stat   = dens_stat,
  d_migr   = dens_migr
)

p_gmm <- ggplot() +
  geom_histogram(
    data    = gmm_data,
    aes(x   = log10_disp, y = after_stat(density)),
    bins    = 60,
    fill    = "grey70",
    col     = "white",
    alpha   = 0.6
  ) +
  geom_line(data = gmm_plot_df, aes(x_log10, d_stat), col = "steelblue",  linewidth = 1) +
  geom_line(data = gmm_plot_df, aes(x_log10, d_migr), col = "firebrick",  linewidth = 1) +
  geom_vline(
    xintercept = log10(DISP_THRESHOLD_KM),
    linetype   = "dashed",
    col        = "black",
    linewidth  = 0.8
  ) +
  annotate(
    "text",
    x     = log10(DISP_THRESHOLD_KM) + 0.05,
    y     = Inf, vjust = 1.5, hjust = 0,
    label = sprintf("threshold = %.1f km", DISP_THRESHOLD_KM),
    size  = 3.5
  ) +
  scale_x_continuous(
    breaks = log10(c(1, 5, 10, 50, 100, 300, 500)),
    labels = c("1", "5", "10", "50", "100", "300", "500")
  ) +
  labs(
    title    = "GMM behavioral classification — T. teniotis fall migration",
    subtitle = "Blue = stationary component | Red = migratory component",
    x        = "Daily displacement (km)",
    y        = "Density"
  ) +
  theme_minimal()

print(p_gmm)

dir.create("../../../Dropbox/MPI/Tadarida/Plots", showWarnings = FALSE, recursive = TRUE)
ggsave(
  "../../../Dropbox/MPI/Tadarida/Plots/displacement_threshold_gmm.png",
  p_gmm, width = 8, height = 5, dpi = 300
)


# =============================================================================
# SECTION 5: Assign behavior labels to all individual-days
# =============================================================================
# Days with displacement = 0 → stationary (same roost)
# Days with displacement > 0 → from GMM posterior
# Days with no fix or filtered out → NA

bats_daily <- bats_daily %>%
  left_join(
    gmm_data %>%
      select(individual_local_identifier, timestamp,
             behavior_gmm, p_migrating, p_stationary),
    by = c("individual_local_identifier", "timestamp")
  ) %>%
  mutate(
    behavior = case_when(
      daily_displacement_km == 0      ~ "stationary",
      !is.na(behavior_gmm)            ~ behavior_gmm,
      is.na(daily_displacement_km)    ~ NA_character_,
      TRUE                            ~ NA_character_
    ),
    # Hard threshold as an alternative/check on the GMM classification
    behavior_thresh = case_when(
      daily_displacement_km == 0                           ~ "stationary",
      daily_displacement_km >= DISP_THRESHOLD_KM           ~ "migrating",
      daily_displacement_km <  DISP_THRESHOLD_KM           ~ "stationary",
      TRUE                                                 ~ NA_character_
    )
  )

# Agreement between GMM and threshold approaches
agree_tbl <- bats_daily %>%
  filter(!is.na(behavior), !is.na(behavior_thresh)) %>%
  count(behavior, behavior_thresh) %>%
  mutate(pct = round(100 * n / sum(n), 1))
print(agree_tbl)


# =============================================================================
# SECTION 6: NSD trajectory plots per individual
# =============================================================================

p_nsd <- bats_daily %>%
  filter(!is.na(nsd_km2)) %>%
  ggplot(aes(
    x   = timestamp,
    y   = sqrt(nsd_km2),   # root-NSD in km for interpretable y-axis
    col = behavior,
    group = individual_local_identifier
  )) +
  geom_line(col = "grey70", linewidth = 0.4) +
  geom_point(size = 1.2, alpha = 0.7) +
  scale_color_manual(
    values = c("stationary" = "steelblue", "migrating" = "firebrick"),
    na.value = "grey50"
  ) +
  facet_wrap(~ individual_local_identifier, scales = "free_y") +
  labs(
    title = "Root-NSD trajectories — T. teniotis fall migration",
    x     = "Date",
    y     = "Distance from deployment (km)",
    col   = "Behavior"
  ) +
  theme_minimal(base_size = 8) +
  theme(
    legend.position  = "top",
    strip.text       = element_text(size = 7),
    axis.text.x      = element_text(angle = 30, hjust = 1, size = 6)
  )

ggsave(
  "../../../Dropbox/MPI/Tadarida/Plots/nsd_by_individual.png",
  p_nsd, width = 14, height = 10, dpi = 300
)


# =============================================================================
# SECTION 7: Behavior summary
# =============================================================================

behavior_summary <- bats_daily %>%
  filter(!is.na(behavior)) %>%
  group_by(individual_local_identifier, behavior) %>%
  summarise(n_days = n(), .groups = "drop") %>%
  pivot_wider(names_from = behavior, values_from = n_days, values_fill = 0) %>%
  mutate(
    total_days      = stationary + migrating,
    pct_migrating   = round(100 * migrating / total_days, 1),
    max_nsd_km      = bats_daily %>%
      group_by(individual_local_identifier) %>%
      summarise(m = max(sqrt(nsd_km2), na.rm = TRUE)) %>%
      pull(m)
  ) %>%
  arrange(desc(pct_migrating))

print(behavior_summary)

p_summary <- behavior_summary %>%
  pivot_longer(c(stationary, migrating), names_to = "behavior", values_to = "days") %>%
  ggplot(aes(
    x    = reorder(individual_local_identifier, pct_migrating),
    y    = days,
    fill = behavior
  )) +
  geom_col(position = "stack") +
  scale_fill_manual(values = c("stationary" = "steelblue", "migrating" = "firebrick")) +
  coord_flip() +
  labs(
    title = "Migration activity by individual",
    x     = NULL,
    y     = "Days",
    fill  = "Behavior"
  ) +
  theme_minimal(base_size = 9) +
  theme(legend.position = "top")

ggsave(
  "../../../Dropbox/MPI/Tadarida/Plots/behavior_summary.png",
  p_summary, width = 8, height = 6, dpi = 300
)


# =============================================================================
# SECTION 8: Save outputs
# =============================================================================

save(
  bats_daily,
  bats_loc_filt,
  gmm2,
  gmm_data,
  DISP_THRESHOLD_KM,
  RADIUS_THRESHOLD_M,
  behavior_summary,
  file = "../../../Dropbox/MPI/Tadarida/Data/rdata/displacement_threshold.robj"
)

cat("Script 2 complete.\n")
cat(sprintf("  Radius filter:        <= %.0f km\n", RADIUS_THRESHOLD_M / 1000))
cat(sprintf("  Displacement threshold: %.1f km/day (GMM crossover)\n", DISP_THRESHOLD_KM))
cat(sprintf(
  "  Individuals classified: %d\n",
  n_distinct(bats_daily$individual_local_identifier[!is.na(bats_daily$behavior)])
))
