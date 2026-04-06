library(tidyverse)
library(move2)
library(sf)
library(terra)
library(mapview)

# source("../SigfoxTagPrep/R/mt_add_start.R")
source("../SigfoxTagPrep/R/import_nanofox_movebank.R")
source("../../../Desktop/movebank_login.R")

y <- import_nanofox_movebank(study_id =
                               c(4882437204), daily_method = "daytime_only"
)
bats_full <- y$full
bats_loc  <- y$location
bats_daily <- y$daily
track_data <- mt_track_data(bats_full)
track_data$deploy_on_location %>% as.data.frame()

bats_full$sigfox_computed_location_radius
ml <- mt_track_lines(bats_full)

mapview(ml, zcol = "individual_local_identifier", legend = FALSE) +
mapview(bats_full, zcol = "individual_local_identifier", cex = "sigfox_computed_location_radius", legend = FALSE)

world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
extent_final <- terra::ext(-10,2,35,45) #terra::ext(bats_loc)
buffer_final <- 1

bats_loc$year <- year(bats_loc$timestamp)
bats_loc$season <- ifelse(month(bats_loc$timestamp) <= 7, "Spring", "Fall")

track_data <- bats_loc %>% mt_track_data()
# plot(track_data$deploy_on_location)
# View(track_data)


ggplot()+
  geom_sf(data = world, col = "black", fill = "lightgrey")+
  geom_path(data = bats_loc, aes(lon, lat, group = individual_local_identifier), col = "black")+
  # geom_sf(data = track_data$deploy_on_location, col = "gold")+
  geom_sf(data = bats_loc, aes(col = individual_local_identifier, alpha = 0.5))+
  geom_sf(data = track_data$deploy_on_location, col = 1, cex = 2)+
  xlab("Longitude")+ylab("Latitude") +
  coord_sf(xlim = c(extent_final[1]-buffer_final, extent_final[2]+buffer_final),
           ylim = c(extent_final[3]-buffer_final, extent_final[4]+buffer_final)) +
  theme_minimal()+
  theme(axis.text = element_text(size = 10, colour = 1),
        legend.text = element_text(size = 10, colour = 1),
        legend.title = element_text(size = 10, colour = 1),
        legend.position = "none") #+
ggsave(filename = "../../../Dropbox/MPI/Tadarida/Plots/Tadarida_map.png")
# facet_wrap(~year+season)

save(y, bats_full, bats_loc, bats_daily, track_data, world,
     file = "../../../Dropbox/MPI/Tadarida/Data/rdata/move_icarus_tt_daytime_daily.robj")
load("../../../Dropbox/MPI/Tadarida/Data/rdata/move_icarus_tt_daytime_daily.robj")

# Check the radius distribution — decide your accuracy cutoff
summary(bats_full$sigfox_computed_location_radius)
hist(bats_full$sigfox_computed_location_radius, breaks = 50)

bats_daily %>%
  group_by(individual_local_identifier) %>%
  summarise(
    n_days = n(),
    first = min(timestamp),
    last  = max(timestamp),
    lat_range = max(lat) - min(lat)
  ) %>%
  arrange(desc(lat_range))
