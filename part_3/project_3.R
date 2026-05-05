# AUTHOR / EMAIL: Aya Surheyao (as4789@princeton.edu)
# CREATION DATE: 03/16/2026
# MODIFIED DATE: 04/30/2026
# PURPOSE: remake PCAs of neutral/unlinked refiltered data


library(ggplot2) ## for plotting
library(sf) ## for presenting shapefile data
library(geosphere) ## for calculate geographic distances
library(R.utils)
library(ecodist) ## for performing mantel test
library(readxl) # for reading Excel files
library(dplyr) # for data manipulation
library(wesanderson) # to get some fun colors
library(paletteer) ## get some more fun colors
library(cowplot) # for plotting multiple plots in one
library(tidyr) # for modifying dataframes
library(PROJ) # converting coordinate system into most common longitude/latitude system

# library("ggmap") 
# library("dartR")
# library("dartRverse")
# library("vcfR")
# library("vegan")

# set up plot aesthetics we can use for all figures
eeb_theme <- function(base_size = 11, base_family = "TT Arial") {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      # Set all text to Arial
      text = element_text(family = "TT Arial"),
      # Remove grid lines and background
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "transparent", color = NA),
      plot.background = element_rect(fill = "transparent", color = NA),
      # Plot title
      plot.title = element_text(size = base_size, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = base_size, hjust = 0.5, margin = margin(b = 10)),
      # Axes
      axis.line = element_line(color = "black"),
      axis.text = element_text(size = 11, color = "black"),
      axis.title = element_text(size = 11, face = "bold"),
      # Legend
      legend.background = element_rect(fill = "transparent", color = NA),
      legend.key = element_blank(),
      legend.text = element_text(size = 11),
      legend.title = element_text(size = 11),
    ) 
}

## working directory
setwd("C:/Users/sophi/OneDrive/Documents/EEB331") #set your wd

####################################################################################
################ Part 0: Using PCA Code from Class  #################################

#### Now read in our newly filtered PCA files ####
pve = read.table("data/pve.txt", header=F)
head(pve) # record % variance explained on PC1 and PC2, control F to add updated values to all PCA plots!
pcs = read.table("data/pcs.txt", header=T)

# how many individuals (rows) do we have now?
dim(pcs) # 221 marmots

############################################################################
#### data cleaning #### 
# we can see here that column 1 is redundant and the first 3 rows don't have valid IDs, so let's drop them
# let's also edit the column of unique animal IDs to match the "uid" format of our metadata for merging
pcs <- pcs %>%
  select(-FID) %>%
  rename(uid = 1) %>%
  slice(-(1:3)) %>%
  mutate(uid = gsub("^UID_|_sorted$", "", uid)) #find the beginning and end of each string -> replace with whitespace in the column 'uid'

View(pcs) ## at this point we have 218 individuals (3 lost with invalid IDs)

# we'll also need to clean up our metadata a bit before we proceed
# we can ignore some of the columns for now, like 'furmark', but no need to drop them
# picnic_lower and picnic_upper are the same colony, so edit this be replacing all white space after picnic_
# follow the same logic to clean up the names for the mm_ and river_ colonies
## southmound, rivermound, sagemound, and middle mound are all included in river colony
## consider cliff as part of Marmot Meadow due to low sample size


## information from originally 236 marmots
marm_meta <- read.csv("data/radseq_pedigree_metadata.csv") |>
  select(uid, col_area)


length(unique(marm_meta$col_area)) ## starting with 17 colony names including redundant ones

## designate up and down-valley by individual
down_valley <- c("Avalanche", "Bench", "RiverAnnex", "River", 
                 "Horsemound", "GothicTown")
up_valley <- c("MarmotMeadow", "Picnic", "NorthPicnic", "Boulder")

## clean up column names for joining
marm_meta_cleaned <- marm_meta %>%
  mutate(col_area = gsub("picnic_(lower|upper|middle)$", "Picnic", col_area)) %>%
  mutate(col_area = gsub("mm_(aspen|maintalus)$", "MarmotMeadow", col_area)) %>%
  mutate(col_area = gsub("river_(rivermound|sagemound|southmound|middlemound)$", "River", col_area)) %>%
  mutate(col_area = gsub("northpk", "NorthPicnic", col_area)) %>%
  mutate(col_area = gsub("rvannex", "RiverAnnex", col_area)) %>%
  mutate(col_area = gsub("gothictown", "GothicTown", col_area)) %>%
  mutate(col_area = gsub("cliff_(lower|upper|middle)$", "MarmotMeadow", col_area)) %>%
  
  ## re-format colony name for consistency and later merging
  mutate(col_area = capitalize(col_area)) |>
  
  ## add up vs. down-valley designations based on colony
  mutate(region = case_when(
    col_area %in% down_valley ~ "down",
    col_area %in% up_valley ~ "up")) |>
  mutate(region = as.factor(region)) 

# head(marm_meta_cleaned)
length(unique(marm_meta_cleaned$col_area)) ## cleaned dataset has 10 distinct sites

# now MERGE new PCA data with new metadata so we can easily make updated plots!
merged_data <- merge(pcs, marm_meta_cleaned, by.x = "uid", by.y = "uid")
View(marm_meta_cleaned)

#### PCA by colony  ####
colony <- ggplot(merged_data, aes(x = PC1, y = PC2)) + 
  geom_point(aes(fill=col_area), color = "black", size = 3, pch = 21) + 
  labs(title = "PCA Analysis by Colony Area",
       x = "PC 1 (6.82%)",
       y = "PC 2 (5.19%)",
       fill = "Colony", caption = "PC1: F = 67.81, p = <2e-16, PC2:  F = 72.43, p = <2e-16") +
  eeb_theme() +
  theme(legend.box.margin = margin(t = 8, r = 8))
colony

# anovas: PCs x colony 
col_aov_pc1 <- aov(merged_data$PC1 ~ as.factor(merged_data$col_area)) 
summary(col_aov_pc1)  # F = 67.81, p = <2e-16 ***
col_aov_pc2 <- aov(merged_data$PC2 ~ as.factor(merged_data$col_area))  
summary(col_aov_pc2)  # F = 72.43, p = <2e-16 ***

## PCA by region (up vs. down) 
region <- ggplot(merged_data, aes(x = PC1, y = PC2)) + 
  geom_point(aes(fill=region), color = "black", size = 3, pch = 21) + 
  labs(title = "PCA Analysis by Region", # rename and add PC %s
       x = "PC 1 (6.82%)",
       y = "PC 2 (5.19%)",
       fill = "Valley Region") +
  eeb_theme()
region

# anovas: PCs x region
reg_aov_pc1 <- aov(merged_data$PC1 ~ as.factor(merged_data$region)) 
summary(col_aov_pc1)  # F = 137.5, p = <2e-16 ***
reg_aov_pc2 <- aov(merged_data$PC2 ~ as.factor(merged_data$region))  
summary(col_aov_pc2)  # F = 65.15, p = 4.8e-14 ***

###############################################################################################
################# Part I: Reading in and Cleaning Up Geospatial Data ##########################


## read in shape file, take out colonies we have no observations from
colony_sf <- st_read("part_3/marmot_polygons_wgs84.shp") |>
  select(Site, Area_sqm, geometry)

## get centroids
colony_centroids <- st_centroid(colony_sf) |>
  ## convert current UTM coordinates into longitude and latitude
  proj_trans(target_crs = "EPSG:4326", source_crs = "+proj=utm +zone=13 +datum=WGS84")

## separate longitude and latitude into separate coloumns
colony_centroids <- colony_centroids |>
  cbind(colony_centroids, st_coordinates(colony_centroids)) |>
  select(Site, Area_sqm, X,Y) 


# write.csv(colony_centroids, "colony_locations.csv")

###################################################################################################
#################### PART II: combining metadata and geographic information #######################


## combine metadata with colony coordinates -> assigns individuals to their colonies!!
marm_by_colony <- merged_data |>
  left_join(colony_centroids, by = join_by(col_area == Site)) |>
  select(uid, col_area, Area_sqm, X, Y, region)

## clean up this plot
centers <- marm_by_colony |>
  select(X, Y, col_area, region) |>
  distinct() |>
  ggplot(aes(x = X, y = Y, color = col_area, shape = region)) +
  geom_point(size = 2) +
  # geom_text(aes(label = col_area,
  # ), size = 2, hjust = -1, vjust = -2) +
  # geom_sf(aes(color = Site)) +
  # geom_sf_label(aes(label = Site), size = 2, nudge_x = 0.01) +
  labs(title = "Geographic Location of Colony Centers", x = "Longitude", y = "Latitude", 
       color = "Colony") +
  theme_minimal() +
  eeb_theme() +
  scale_colour_paletteer_d("ggthemes::Classic_10") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8), 
        legend.key.spacing = unit(0.005, 'cm'))

centers
ggsave("colony_centers.png", plot = centers, w = 5, height = 4, units = "in", dpi = 500, bg = "white")


## get list of IDs for down-valley
up_valleys <- merged_data |>
  filter(region == "up") |>
  select(uid, region)
up_ids <- up_valleys$uid

## get list of IDs for down-valley
down_valleys <- merged_data |>
  filter(region == "down") |>
  select(uid, region)
down_ids <- down_valleys$uid

## EXPLORATORY: plot number of marmots by colony
m_distribution <- marm_by_colony |>
  group_by(col_area)|>
  mutate(individuals = n()) |>
  ungroup() |>
  select(col_area,individuals) |>
  ggplot(aes(x = col_area, fill = col_area)) +
  geom_bar(color = "black") +
  labs(title = "Individuals by Colony", subtitle = "Metadata from 218 Individuals", 
       x = "Collection Area", y = "Number of Marmots in Sample") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8), legend.position = "none")
m_distribution

## EXPLORATORY: show differences in the size of collection areas (maybe useful to visualize this in a different way)
size_of_colony <-  marm_by_colony |>
  select(col_area, Area_sqm) |>
  distinct() |>
  ggplot(aes(x = col_area, y = Area_sqm, fill = col_area)) +
  geom_col(color = "black") +
  labs(title = "Individuals by Colony", subtitle = "Metadata from 218 Individuals", 
       x = "Collection Area", y = "Size of Collection Area (sqm)") +
  eeb_theme() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8), legend.position = "none")

size_of_colony


## create data-frame with correctly formatted latitudes and longitudes to create distance matrix
colony_coord <- data.frame(marm_by_colony$X, marm_by_colony$Y)

## make matrix of pairwise geographic distances!
pair_distances <- data.frame(distm(colony_coord, fun = distHaversine))


##########################################################################################
########## PART III: Genetic Relatedness and Isolation by Distance Analysis ##############

# read in genetic relatedness data provided from Canvas
load("part_3/rel.RData")
re_vals <- get("rel")[1]$relatedness


###### cleaning up indices on related data frame #####
## removing individuals not included in our 218 subset
## changing ID labels to only contain the 9-character center so all ids are of the format: XXXX_XXXX
## removing marmots with missing labels
re_vals <- re_vals |>
  filter(!is.na(ind1.id) | !is.na(ind2.id)) |>
  filter(group == "UIUI") |>
  mutate(ind1.id = gsub("^UID_|_sorted$", "", ind1.id)) |>
  mutate(ind2.id = gsub("^UID_|_sorted$", "", ind2.id)) |>
  mutate(ind1.id = gsub("|_sort$", "", ind1.id)) |>
  mutate(ind2.id = gsub("|_sort$", "", ind2.id)) |>
  mutate(ind1.id = gsub("|_sorte$", "", ind1.id)) |>
  mutate(ind2.id = gsub("|_sorte$", "", ind2.id)) |>
  mutate(ind1.id = gsub("^RMBL_|_(marmot)", "", ind1.id)) |>
  mutate(ind2.id = gsub("^RMBL_|_(marmot)", "", ind2.id)) |>
  mutate(ind1.id = gsub("RMBL_marmot_no_label", "NA", ind1.id)) |>
  mutate(ind2.id = gsub("RMBL_marmot_no_label", "NA", ind2.id))

#### rename the columns of pair_distances ####

## extract all individual marmot id labels
id_names = marm_by_colony$uid
f_id1 = unique(re_vals$ind1.id)
f_id2 = unique(re_vals$ind2.id)
f_id1
f_id2

## assign column names and rownames in the distance pairwise matrix to match individual ID labels
## a given entry in the matrix is the distance betweeen [column name] individuals and [row name] individual
colnames(pair_distances) = id_names
rownames(pair_distances) = id_names

### check that row names were changed and in the correct order
colnames(pair_distances)
rownames(pair_distances)


######## combining data frame with relatedness info and geographic info ######

## numbers of rows in relatedness matrix
rel_length <- length(rownames(re_vals)) ## 23,635

## make new condensed pairwise distances for plotting and mantel testing
ibd_combined <- data.frame()

## iterate through the length of the relatedness data-frame row-wise
for (i in seq(rel_length)) {
  
  ## get id values to match to
  ind1 = re_vals$ind1.id[i]
  ind2 = re_vals$ind2.id[i]
  
  ## find the values corresponding in the distance matrix and record the index
  x_dist <- match(ind1, rownames(pair_distances))
  y_dist <- match(ind2, colnames(pair_distances))
  
  
  ### add new row to data frame containing the two indices, geographic distance, and relatedness
  row <- data.frame(id_1 = ind1, 
                    id_2 = ind2, 
                    distance = pair_distances[x_dist,y_dist],
                    ritland = re_vals$ritland[i], 
                    wang = re_vals$wang[i],  
                    trioml = re_vals$trioml[i], 
                    dyadml = re_vals$dyadml[i]) 
  ibd_combined <- rbind(ibd_combined, row)
}


### re-assign colonies to individuals in new dataset
all_ids <- marm_by_colony$uid
col_assignment <- marm_by_colony$col_area

## re-assign up-down and colony assignments
ibd_combined <- ibd_combined |>
  ## regions
  mutate(ind1_region = case_when(
    id_1%in% up_ids ~ as.factor("up"),
    id_1 %in% down_ids ~ as.factor("down"))) |>
  mutate(ind2_region = case_when(
    id_2 %in% up_ids ~ as.factor("up"), 
    id_2 %in% down_ids ~ as.factor("down"))) |>
  
  ## colonies
  mutate(ind1_colony = col_assignment[match(id_1, all_ids)]) |>
  mutate(ind2_colony = col_assignment[match(id_2, all_ids)]) |>
  
  ## drop any rows with missing values
  drop_na()


#### create different subsets for analysis ####

## subset containing all pairwise comparisons between different-colonies individuals
ibd_diff_colony <- ibd_combined |>
  filter(distance != 0)

## containing all pairwise comparisons between up-valley colonies
ibd_up <- ibd_combined |>
  filter(ind1_region == "up") |>
  filter(ind2_region == "up")

## containing all pairwise comparisons between down-valley colonies
ibd_down <- ibd_combined |>
  filter(ind1_region == "down") |>
  filter(ind2_region == "down")

## containing all comparisons between and up-valley individual and a down-valley individual
ibd_up_down <- ibd_combined |>
  filter(ind1_region != ind2_region)


##################################### Plotting ###################################

## violin plot by relatedness
ibd_combined |>
  ## create three categories: up-up comparisons, down-down comparisons, up-down comparisons
  mutate(comparison = case_when(
    ind1_region != ind2_region ~ as.factor("inter"),
    ind1_region == "up" & ind2_region == "up" ~ as.factor("up"),
    ind1_region == "down" & ind2_region == "down" ~ as.factor("down")
  )) |>
  ggplot(aes(x = comparison, y = ritland, fill = comparison)) +
  geom_violin() +
  eeb_theme()


## plot relatedness vs. geographic distance for all observations
ibd1 <- ibd_combined |>
  filter(!is.na(ind1) && !is.na(ind2)) |>
  distinct() |>
  ggplot(aes(x = distance, y = ritland)) +
  geom_point(color = "slateblue", alpha = 0.5) +
  geom_smooth(formula = "y ~ x",
              method = "lm", color = "black", se = TRUE) +
  labs(title = "Isolation by distance", x = "Distance in meters", 
       y = "Relatedness", 
       subtitle = "All pairwise comparisons between 218 individuals", 
       caption = "Mantel Test: r = -0.4839, p = 0.0001") +
  
  eeb_theme()

ibd1

## plotting for different-colony subset
ibd2 <- ibd_diff_colony |>
  distinct() |>
  ggplot(aes(x = distance, y = Ritland)) +
  geom_point(color = "slateblue", alpha = 0.5) +
  labs(title = "Isolation by distance", x = "Distance in meters", 
       y = "Relatedness", 
       subtitle = "Excluding same-colony comparisons") +
  eeb_theme()
ibd2


### up-valley subset
ibd3 <- ibd_up |>
  filter(!is.na(ind1) && !is.na(ind2)) |>
  distinct() |>
  ggplot(aes(x = distance, y = Ritland)) +
  geom_point(color = "darkgoldenrod", alpha = 0.5) +
  geom_smooth(formula = "y ~ x",
              method = "lm", color = "black", se = TRUE) +
  labs(title = "Isolation by distance", x = "Distance in meters", 
       y = "Relatedness", 
       subtitle = "Comparisons between up-valley individuals", 
       caption = "Mantel Test: r = -0.5029, p = 0.0001") +
  eeb_theme()
ibd3


## down-valley subset
ibd4 <- ibd_down |>
  filter(!is.na(ind1) && !is.na(ind2)) |>
  distinct() |>
  ggplot(aes(x = distance, y = Ritland)) +
  geom_point(color = "seagreen", alpha = 0.5) +
  geom_smooth(formula = "y ~ x",
              method = "lm", color = "black", se = TRUE) +
  labs(title = "Isolation by distance", x = "Distance in meters", 
       y = "Relatedness", 
       subtitle = "Comparisons between down-valley individuals", 
       caption = "Mantel Test: r = -0.3876, p = 0.0001") +
  eeb_theme()
ibd4

## up vs. down comparisons subset
ibd5 <- ibd_up_down |>
  filter(!is.na(ind1) && !is.na(ind2)) |>
  distinct() |>
  ggplot(aes(x = distance, y = Ritland)) +
  geom_point(color = "darkgoldenrod", alpha = 0.5) +
  labs(title = "Isolation by distance", x = "Distance in meters", 
       y = "Relatedness", 
       subtitle = "All up vs. down-colony comparisons") +
  theme_minimal()
ibd5

### create high-quality individual isolation by distance plots 
ggsave("isolation_by_distance_part1.png", plot = ibd1, width = 5, height = 4, units = "in", dpi = 500, bg = "white")
ggsave("isolation_by_distance_part2.png", plot = ibd3, width = 5, height = 4, units = "in", dpi = 500, bg = "white")
ggsave("isolation_by_distance_part3.png", plot = ibd4, width = 5, height = 4, units = "in", dpi = 500, bg = "white")


### create high quality plot for PCAs
fig1 <- plot_grid(colony, region, nrow = 1, ncol = 2)
ggsave("landscape_genomics_draft_fig1.png", plot = fig1, width = 12, height = 5, units = "in", dpi=500, bg="white")


######################## Compute Mantel Statistics #########################

## all observations 
mantel_base <- mantel(formula = distance~ritland, data = ibd_combined |>
                        select(distance, ritland), mrank = TRUE, nperm = 9999)
mantel_base

## mantel r coefficient:-0.48394754
## pval1 (one-tailed, null r<= 0) -> p = 1 (fail to reject)
## pval2 (one-tailed, null r>= 0) -> p = 0.00010001 (reject)
## pval3(two-tail, null r = 0) -> p = 0.00010001 (reject)


### unknown error when running this test related to data frame size
# mantel(formula = distance~relatedness, data = ibd_diff_colony |>
#          select(distance, relatedness), mrank = TRUE, nperm = 9999)

## up-valley observations 
mantel(formula = distance~ritland, data = ibd_up |>
         select(distance, ritland), mrank = TRUE, nperm = 9999)

## mantel r coefficient: -0.50294206
## pval1 (one-tailed, null r<= 0) -> p = 1 (fail to reject)
## pval2 (one-tailed, null r>= 0) -> p = 0.00010001 (reject)
## pval3(two-tail, null r = 0) -> p = 0.00010001 (reject)

## down-valley observations
mantel(formula = distance~ritland, data = ibd_down |>
         select(distance, ritland), mrank = TRUE, nperm = 9999)

## mantel r coefficient: -0.38761603 
## pval1 (one-tailed, null r<= 0) -> p = 1 (fail to reject)
## pval2 (one-tailed, null r>= 0) -> p = 0.00010001 (reject)
## pval3(two-tail, null r = 0) -> p = 0.00010001 (reject)

### unknown error when running this test related to data frame size
# mantel(formula = distance~relatedness, data = ibd_up_down |>
#          select(distance, relatedness), mrank = TRUE, nperm = 9999)



## compare average Fst by colonies and do a chi-squared test?
## look at area vs. average relatedness within colonies?


