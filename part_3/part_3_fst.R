# AUTHOR / EMAIL: Aya Surheyao (as4789@princeton.edu)
# CREATION DATE: 03/16/2026
# MODIFIED DATE: 04/06/2026
# PURPOSE: remake PCAs of neutral/unlinked refiltered data


library("ggplot2")
library("sf")
library("dplyr")
library("vegan")
library("geosphere")
library("R.utils")
# library("dartR")
# library("dartRverse")
# library("vcfR")


################# Part 1: Reading in and Cleaning Data ##########################
## working directory
setwd("C:/Users/sophi/OneDrive/Documents/EEB331") #set your wd

## read in shape file, take out colonies we have no observations from
colony_sf <- st_read("part_3/marmot_polygons_wgs84.shp") |>
  select(Site, Area_sqm, geometry)

library("PROJ")

dst <- "+proj=longlat +datum=WGS84"
src <- "+proj=utm +datum=WGS84"


## get centroids
colony_centroids <- st_centroid(colony_sf) |>
  select(Site, Area_sqm)

colony_centroids <- colony_centroids |>
  cbind(colony_centroids, st_coordinates(colony_centroids)) |>
  select(Site, Area_sqm) |>
  proj_trans(target_crs = "EPSG:4326", source_crs = "+proj=utm +zone=13 +datum=WGS84")


## separate longitude and latitude into separate coloumns
colony_centroids <- colony_centroids |>
  cbind(st_coordinates(colony_centroids))

## overlay with google earth (come back to this)

colony_sf |>
  ggplot() +
  geom_sf(aes(fill = Site))

## clean up this plot
centers <- colony_centroids |>
  filter(Site != "Avery") |>
  filter(Site != "Bellview") |>
  filter(Site != "Stonefield") |>
  ggplot() +
  geom_sf(aes(color = Site)) +
  labs(title = "Geography  of Colony Centers") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8)) +
  theme(title = element_text(vjust = -93))

centers

## show on google maps 
# library("ggmap")
# my_map <- ggmap(get_googlemap(center = c( lon = colony_centroids$X[1], lat = colony_centroids$Y[1]), maptype = 'terrain'))

## read in metadata
marm_meta <- read.csv("data/radseq_pedigree_metadata.csv") |>
  select(uid, col_area)

## clean up column names for joining
marm_meta <- marm_meta %>%
  mutate(col_area = gsub("picnic_(lower|upper|middle)$", "Picnic", col_area)) %>%
  mutate(col_area = gsub("mm_(aspen|maintalus)$", "MarmotMeadow", col_area)) %>%
  mutate(col_area = gsub("river_(rivermound|sagemound|southmound|middlemound)$", "River", col_area)) %>%
  mutate(col_area = gsub("northpk", "NorthPicnic", col_area)) %>%
  mutate(col_area = gsub("rvannex", "RiverAnnex", col_area)) %>%
  mutate(col_area = gsub("gothictown", "GothicTown", col_area)) %>%
  mutate(col_area = gsub("cliff_(lower|upper|middle)$", "Cliff", col_area)) %>%
  mutate(col_area = capitalize(col_area))

## combine metadata with colony coordinates -> assigns individuals to their colonies!!
marm_by_colony <- left_join(marm_meta, colony_centroids, by = join_by(col_area == Site))

unique(marm_by_colony$col_area)

## plot  number of marmots by colony

m_distribution <- marm_by_colony |>
  group_by(col_area)|>
  mutate(individuals = n()) |>
  ungroup() |>
  select(col_area,individuals) |>
  ggplot(aes(x = col_area, fill = col_area)) +
  geom_bar() +
  labs(title = "Individuals by Colony", subtitle = "Metadata from 236 Individuals", 
       x = "Collection Area", y = "Number of Marmots in Sample") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))


m_distribution

marm_by_colony |>
  group_by(col_area)|>
  mutate(individuals = n()) |>
  ungroup() |>
  select(col_area,individuals) |>
  distinct() |>
  mutate(total = sum(individuals))

## make matrix of geographic distance
colony_coord <- data.frame(marm_by_colony$X, marm_by_colony$Y)
pair_distances <- data.frame(distm(colony_coord, fun = distHaversine))


# create FST



## do something logistic to look at up vs. down-valley? does it make sense to use F-statistics instead of relatedness?
## perform mantel test


## compare average Fst by colonies and do a chi-squared test


## figure out which metric is the one that we would need, probably:
## lynchrd, quellergt, lynchli, or wang though they're all pretty close

## first plot: relatedness vs. distance
