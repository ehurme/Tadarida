# explore Tadarida teniotis tracks

library(tidyverse)
library(move2)
library(sf)
library(terra)

source("../../../Desktop/movebank_login.R")
source("../SigfoxTagPrep/R/mt_add_start.R")

df <- movebank_download_study(study_id = 4882437204)
data <- mt_track_data(df)

# add start locations
df <- mt_add_start(df)

# add lat long
coords <- st_coordinates(df)
df$lon <- coords[,1]
df$lat <- coords[,2]


# add displacement
library(dplyr)
library(sf)

# Make sure rows are ordered within tracks
df <- df %>%
  arrange(individual_local_identifier, timestamp)

# (Optional) drop empty geometries
df_clean <- df %>%
  filter(!st_is_empty(geometry))

# Compute displacement from first point of *each track*
df_disp <- df_clean %>%
  group_by(individual_local_identifier) %>%
  mutate(
    first_geom = geometry[1],  # first location of this individual
    distance = c(NA, mt_distance(.)),
    displacement = st_distance(geometry, first_geom, by_element = TRUE), # units (m)
    displacement_km = as.numeric(displacement) / 1000
  ) %>%
  ungroup()

df_disp$displacement_km %>% plot()

df_clean <- df_clean %>% mutate(distance = mt_distance(.)/1000)
df_clean$distance
df_clean <- df_clean %>% group_by(individual_local_identifier) %>%
  mutate(displacement = cumsum(distance))

plot(df_clean$displacement, df_disp$displacement_km, col = factor(df_clean$individual_local_identifier))

displacement <- units::set_units(displacement, "km")
NSD <- displacement^2
NSD <- units::set_units(NSD, "km^2")
