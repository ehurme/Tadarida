# =============================================================================
# Script 3 (R wrapper): Export → Python env extraction → Re-import
# =============================================================================
# This script:
#   (1) Exports the classified bats_daily from Script 2 to a CSV that the
#       Python script can read.
#   (2) Adds a movement bearing column if not already present (needed for
#       tailwind/crosswind decomposition in Python).
#   (3) Calls 3_extract_env_data.py via system2().
#   (4) Reads the enriched CSV back into R and merges onto bats_daily.
#   (5) Saves the final enriched object for Script 4+.
# =============================================================================

library(tidyverse)
library(move2)
library(sf)
library(geosphere)

load("../../../Dropbox/MPI/Tadarida/Data/rdata/displacement_threshold.robj")

# ── 1. Compute movement bearing (day-to-day) ─────────────────────────────────
# bearingRhumb() gives the forward azimuth from position(t-1) → position(t),
# which is what we want for tailwind decomposition relative to travel direction.
# This must be done before export; Python does not have the move2 context.

if (!"lon" %in% names(bats_daily)) {
  coords     <- sf::st_coordinates(bats_daily)
  bats_daily <- bats_daily %>% mutate(lon = coords[, 1], lat = coords[, 2])
}

bats_daily <- bats_daily %>%
  arrange(individual_local_identifier, timestamp) %>%
  group_by(individual_local_identifier) %>%
  mutate(
    lon_lag  = lag(lon),
    lat_lag  = lag(lat),
    bearing  = case_when(
      is.na(lon_lag) ~ NA_real_,
      TRUE ~ geosphere::bearingRhumb(
        cbind(lon_lag, lat_lag),
        cbind(lon,     lat)
      )
    )
  ) %>%
  ungroup()

# ── 2. Export to CSV ──────────────────────────────────────────────────────────
export_path <- "../../../Dropbox/MPI/Tadarida/Data/bats_daily_export.csv"

bats_daily %>%
  as_tibble() %>%
  select(
    individual_local_identifier,
    timestamp,
    lon, lat,
    behavior,
    daily_displacement_km,
    bearing,
    nsd_km2
  ) %>%
  # Ensure timestamp is character ISO-8601 UTC (Python reads this cleanly)
  mutate(timestamp = format(timestamp, "%Y-%m-%dT%H:%M:%SZ")) %>%
  write_csv(export_path)

cat(sprintf("Exported %d rows to %s\n", nrow(bats_daily), export_path))

# ── 3. Run Python extraction ──────────────────────────────────────────────────
# Requires: pip install requests pandas tqdm numpy shapely
# Adjust the python path if using conda/venv.

python_bin    <- "python3"   # or "python", or full path to your env
python_script <- here::here("3_extract_env_data.py")

cat("Running Python env extraction...\n")
exit_code <- system2(
  python_bin,
  args   = shQuote(python_script),
  stdout = "",   # print to console
  stderr = ""
)

if (exit_code != 0) {
  stop(sprintf(
    "Python script exited with code %d. Check output above for errors.", exit_code
  ))
}

# ── 4. Read result back ───────────────────────────────────────────────────────
env_path <- "../../../Dropbox/MPI/Tadarida/Data/bats_daily_env.csv"

bats_env <- read_csv(env_path, show_col_types = FALSE) %>%
  mutate(timestamp = as.POSIXct(timestamp, tz = "UTC"))

cat(sprintf(
  "Env data loaded: %d rows × %d columns\n",
  nrow(bats_env), ncol(bats_env)
))

# ── 5. Quick QC ──────────────────────────────────────────────────────────────
# How many individual-days have complete env data?
n_missing_wind <- sum(is.na(bats_env$night_mean_wind_speed_10m))
n_missing_pbl  <- sum(is.na(bats_env$night_mean_boundary_layer_height))

cat(sprintf(
  "Missing night wind speed: %d / %d (%.1f%%)\n",
  n_missing_wind, nrow(bats_env),
  100 * n_missing_wind / nrow(bats_env)
))
cat(sprintf(
  "Missing PBL height: %d / %d (%.1f%%)\n",
  n_missing_pbl, nrow(bats_env),
  100 * n_missing_pbl / nrow(bats_env)
))

# Spot-check: distribution of dusk pressure tendency
summary(bats_env$dusk_pressure_tendency)

# ── 6. Quick visualization ────────────────────────────────────────────────────
p_wind <- bats_env %>%
  filter(!is.na(behavior), !is.na(night_mean_tailwind_10m)) %>%
  ggplot(aes(night_mean_tailwind_10m, fill = behavior)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c("stationary" = "steelblue", "migrating" = "firebrick")) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(
    title = "Nocturnal tailwind by behavioral state",
    subtitle = "Positive = wind assists movement direction",
    x = "Mean nocturnal tailwind (m/s)",
    fill = "Behavior"
  ) +
  theme_minimal()

p_pbl <- bats_env %>%
  filter(!is.na(behavior), !is.na(night_mean_boundary_layer_height)) %>%
  ggplot(aes(night_mean_boundary_layer_height, fill = behavior)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c("stationary" = "steelblue", "migrating" = "firebrick")) +
  labs(
    title = "Planetary boundary layer height by behavioral state",
    subtitle = "Higher PBL → more convective mixing → generally favourable for active flight",
    x = "Mean nocturnal PBL height (m)",
    fill = "Behavior"
  ) +
  theme_minimal()

print(p_wind)
print(p_pbl)

ggsave("../../../Dropbox/MPI/Tadarida/Plots/tailwind_by_behavior.png",
       p_wind, width = 8, height = 5, dpi = 300)
ggsave("../../../Dropbox/MPI/Tadarida/Plots/pbl_by_behavior.png",
       p_pbl, width = 8, height = 5, dpi = 300)

# ── 7. Save ───────────────────────────────────────────────────────────────────
save(
  bats_env,
  file = "../../../Dropbox/MPI/Tadarida/Data/rdata/bats_daily_env.robj"
)

cat("Script 3 complete. Object saved: bats_daily_env\n")
