# AUTHOR / EMAIL: Aya Surheyao (as4789@princeton.edu)
# CREATION DATE: 03/16/2026
# MODIFIED DATE: 04/06/2026
# PURPOSE: remake PCAs of neutral/unlinked refiltered data


library(ggplot2) ## for plotting
library(sf) ## for presenting shapefile data
library(geosphere) ## for calculate geographic distances
library(R.utils)
library(ecodist) ## for performing mantel test
library(readxl) # for reading Excel files
library(dplyr) # for data manipulation
library(wesanderson) # to get some fun colors
library(cowplot) # for plotting multiple plots in one

# library("ggmap") 
# library("dartR")
# library("dartRverse")
# library("vcfR")
# library("vegan")

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
## just added southmound to be counted as river

marm_meta <- read.csv("data/radseq_pedigree_metadata.csv") |>
  select(uid, col_area)

length(unique(marm_meta$col_area))

## clean up column names for joining
marm_meta_cleaned <- marm_meta %>%
  mutate(col_area = gsub("picnic_(lower|upper|middle)$", "Picnic", col_area)) %>%
  mutate(col_area = gsub("mm_(aspen|maintalus)$", "MarmotMeadow", col_area)) %>%
  mutate(col_area = gsub("river_(rivermound|sagemound|southmound|middlemound)$", "River", col_area)) %>%
  mutate(col_area = gsub("northpk", "NorthPicnic", col_area)) %>%
  mutate(col_area = gsub("rvannex", "RiverAnnex", col_area)) %>%
  mutate(col_area = gsub("gothictown", "GothicTown", col_area)) %>%
  mutate(col_area = gsub("cliff_(lower|upper|middle)$", "MarmotMeadow", col_area)) %>%
  mutate(col_area = capitalize(col_area)) |>
  mutate(region = case_when(
    col_area %in% down_valley ~ "down",
    col_area %in% up_valley ~ "up")) |>
  mutate(region = as.factor(region)) |>
  filter(!is.na(region))

head(marm_meta_cleaned)
length(unique(marm_meta_cleaned$col_area))

# now MERGE new PCA data with new metadata so we can easily make updated plots!
merged_data <- merge(pcs, marm_meta_cleaned, by.x = "uid", by.y = "uid")
View(marm_meta_cleaned)

## check number of colonies -> 10 colonies (because we combined cliff in MarmotMeadow)
merged_data |>
  select(col_area) |>
  distinct()

#### PCA by colony  ####
colony <- ggplot(merged_data, aes(x = PC1, y = PC2)) + 
  geom_point(aes(fill=col_area), color = "black", size = 3, pch = 21) + 
  labs(title = "PCA Analysis by Colony Area", # rename and add PC %s
       x = "PC 1 (%)",
       y = "PC 2 (%)",
       fill = "Colony") +
  eeb_theme()
colony

# anovas: PCs x colony 
col_aov_pc1 <- aov(merged_data$PC1 ~ as.factor(merged_data$col_area)) 
summary(col_aov_pc1)  # F = 64.63, p = <2e-16 ***
col_aov_pc2 <- aov(merged_data$PC2 ~ as.factor(merged_data$col_area))  
summary(col_aov_pc2)  # F = 62.52, p = <2e-16 *** ***

## PCA by region (up vs. down) 
region <- ggplot(merged_data, aes(x = PC1, y = PC2)) + 
  geom_point(aes(fill=region), color = "black", size = 3, pch = 21) + 
  labs(title = "PCA Analysis by Region", # rename and add PC %s
       x = "PC 1 (%)",
       y = "PC 2 (%)",
       fill = "Colony") +
  eeb_theme()
region

# anovas: PCs x region
col_aov_pc1 <- aov(merged_data$PC1 ~ as.factor(merged_data$region)) 
summary(col_aov_pc1)  # F = 137.5, p = <2e-16 ***
col_aov_pc2 <- aov(merged_data$PC2 ~ as.factor(merged_data$region))  
summary(col_aov_pc2)  # F = 65.15, p = 4.8e-14 ***

################# Part I: Reading in and Cleaning Geospatial Data ##########################

## read in shape file, take out colonies we have no observations from
colony_sf <- st_read("part_3/marmot_polygons_wgs84.shp") |>
  select(Site, Area_sqm, geometry)

## get package for converting coordinate system into most common system
library("PROJ")

## get centroids
colony_centroids <- st_centroid(colony_sf) |>
  ## convert current UTM coordinates into longitude and latitude
  proj_trans(target_crs = "EPSG:4326", source_crs = "+proj=utm +zone=13 +datum=WGS84")

## separate longitude and latitude into separate coloumns
colony_centroids <- colony_centroids |>
  cbind(colony_centroids, st_coordinates(colony_centroids)) |>
  select(Site, Area_sqm, X,Y) 


write.csv(colony_centroids, "colony_locations.csv")

## overlay with google earth (come back to this) -> I need an API key which costs money :(

# str(colony_centroids)

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
  # geom_sf_label(aes(label = Site), size = 2, nudge_x = 0.01) +
  labs(title = "Geography of Colony Centers", x = "Longitude", y = "Latitude") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))

centers
## show on google maps 
# library("ggmap")
# my_map <- ggmap(get_googlemap(center = c( lon = colony_centroids$X[1], lat = colony_centroids$Y[1]), maptype = 'terrain'))


######################################################################################
#################### PART II: metadata and genomic information #######################

  
## combine metadata with colony coordinates -> assigns individuals to their colonies!!

## designate up and down-valley by individual
down_valley <- c("Avalanche", "Bench", "RiverAnnex", "River", 
                 "Horsemound", "GothicTown")
up_valley <- c("MarmotMeadow", "Picnic", "NorthPicnic", "Boulder")


## designate up-valley vs. down-valley sites
marm_by_colony <- merged_data |>
  left_join(colony_centroids, by = join_by(col_area == Site)) |>
  select(uid, col_area, X, Y, region)

## get list of IDs for down-valley
up_valleys <- merged_data |>
  filter(region == "up") |>
  select(uid, region)

up_ids <- up_valleys$uid

down_valleys <- merged_data |>
  filter(region == "down") |>
  select(uid, region)
down_ids <- down_valleys$uid

## plot  number of marmots by colony
#### NOTE: change to 218 individuals to match ANOVA and remaining data?? ##
m_distribution <- marm_by_colony |>
  group_by(col_area)|>
  mutate(individuals = n()) |>
  ungroup() |>
  select(col_area,individuals) |>
  ggplot(aes(x = col_area, fill = col_area)) +
  geom_bar(color = "black") +
  labs(title = "Individuals by Colony", subtitle = "Metadata from 236 Individuals", 
       x = "Collection Area", y = "Number of Marmots in Sample") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8), legend.position = "none")
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


# read in relatedness data
load("part_3/rel.RData")
re_vals <- get("rel")[1]$relatedness

# re_vals |>
#   arrange(desc(trioml)) |>
#   head()


## clean up indices on related dataframe
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

  # select(pair.no, ind1.id, ind2.id, trioml, group)
 
  # could add these ones two back, don't think they show up in the metadata
  # filter(group != "cURM")

  # unique(re_vals$group)

#### rename the columns of pair_distances ####

## names that we want
id_names = marm_by_colony$uid
f_id1 = unique(re_vals$ind1.id)
f_id2 = unique(re_vals$ind2.id)
f_id1
f_id2

## current column names in the distance pairwise matrix
colnames(pair_distances) = id_names
rownames(pair_distances) = id_names

### check that row names were added
colnames(pair_distances)
rownames(pair_distances)


######## combining data frame with relatedness info and geographic info ######

## numbers of rows in relatedness matrix
rel_length <- length(rownames(re_vals)) ## 23,635


## make new condensed pairwise distances
ibd_combined <- data.frame()

## iterate through the length of the relatedness data-frame row-wise
for (i in seq(rel_length)) {
  
  ## get id values to match to
  ind1 = re_vals$ind1.id[i]
  ind2 = re_vals$ind2.id[i]
  
  ## find the values corresponding in the distance matrix and record the index
  x_dist <- match(ind1, rownames(pair_distances), nomatch = 0)
  y_dist <- match(ind2, colnames(pair_distances), nomatch = 0)
  
  
  ### add new row to data frame containing the two indices, geographic distance, and relatedness
  row <- data.frame(id_1 = ind1, 
                    id_2 = ind2, 
                    distance = pair_distances[x_dist,y_dist],
                    relatedness = re_vals$trioml[i]) 
  ibd_combined <- rbind(ibd_combined, row)
}



### 
ibd_diff_colony <- ibd_combined |>
  filter(distance != 0)

ibd_up <- ibd_combined |>

ibd_down
ibd_up_down

######## repeat but only comparing values that are not from the same colony #########

# ## make new condensed pairwise distances
# ibd_diff_colony <- data.frame()
# 
# ## iterate through the length of the relatedness data-frame row-wise
# for (i in seq(rel_length)) {
#   
#   ## get id values to match to
#   ind1 = re_vals$ind1.id[i]
#   ind2 = re_vals$ind2.id[i]
#   
#   ## find the values corresponding in the distance matrix and record the index
#   x_dist <- match(ind1, rownames(pair_distances), nomatch = 0)
#   y_dist <- match(ind2, colnames(pair_distances), nomatch = 0)
#   
#   
#   ### add new row to data frame containing the two indices, geographic distance, and relatedness
#   ### if geographic distance is 0, will not be added
#   if (pair_distances[x_dist,y_dist] != 0) {
#     row <- data.frame(id_1 = ind1, 
#                     id_2 = ind2, 
#                     distance = pair_distances[x_dist,y_dist],
#                     relatedness = re_vals$trioml[i])
#     ibd_diff_colony <- rbind(ibd_combined, row)
#   }
# }


# ibd_combined <- ibd_combined |>
#   mutate(log_relatedness = case_when(
#     relatedness != 0 ~ log(relatedness),
#     ) 


## look now at up vs. down-colony relatedness
ibd_combined <- ibd_combined |>
  mutate(ind1_region = case_when(
    ind1.id %in% up_ids ~ as.factor("up"),
    ind1.id %in% down_ids ~ as.factor("down"))) |>
  mutate(ind2_region = case_when(
    ind2.id %in% up_ids ~ as.factor("up"), 
    ind2.id %in% down_ids ~ as.factor("down"))) |>
  filter(!is.na(ind2_region)) |>
  filter(!is.na(ind1_region))


## plot relatedness vs. geographic distance
ibd <- ibd_combined |>
  filter(!is.na(ind1) && !is.na(ind2)) |>
  distinct() |>
  ggplot(aes(x = distance, y = log_relatedness)) +
  geom_point(color = "seagreen") +
  labs(title = "Pairwise relatedness vs. geographic distance between 
       marmots from RMBL colonies", x = "Distance in Meters", 
       y = "Relatedness (trioml)", 
       subtitle = "Distance Measured Between Colony Centers") +
  theme_minimal()


ibd

fig1 <- plot_grid(centers, colony, ibd, rows = 2, cols = 2)


ggsave("landscape_genomics_draft_fig1.png", plot = fig1, width = 18, height = 14, units = "in", dpi=500, bg="white")


##### compute mantel value #############
mantel(formula = distance~relatedness, data = ibd_combined |>
         select(distance, relatedness), mrank = TRUE, nperm = 9999)



  



## mantel r coefficient: -0.486
## pval1 (one-tailed, null r<= 0) -> p = 1 (fail to reject)
## pval2 (one-tailed, null r>= 0) -> p = 0.00010001 (reject)
## pval3(two-tail, null r = 0) -> p = 0.00010001 (reject)


## do something logistic to look at up vs. down-valley? does it make sense to use F-statistics instead of relatedness?


## compare average Fst by colonies and do a chi-squared test


## figure out which metric is the one that we would need, probably:
## using trioml


