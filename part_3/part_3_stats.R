# AUTHOR / EMAIL: Aya Surheyao (as4789@princeton.edu)
# CREATION DATE: 04/06/2026
# MODIFIED DATE: 04/06/2026
# PURPOSE: remake PCAs of neutral/unlinked refiltered data

library(multcompView)
library(readxl) # for reading Excel files
library(dplyr) # for data manipulation
library(ggplot2) # for plotting
library(wesanderson) # to get some fun colors
library(cowplot) # for plotting multiple plots in one
library(agricolae)
library(broom)
library(gt)
library(tidyr)
## testing assumptions (parametric, etc.)



########################## Q-Q Plot ########################################

# Load necessary libraries
library(tidyverse)
library(metBrewer)

# Calculate slope and intercept for linear regression line
slope <- (cor(ibd_combined$distance, ibd_combined$log_relatedness) * 
            (sd(ibd_combined$distance)) / 
            sd(ibd_combined$log_relatedness))
intercept <- (mean(ibd_combined$log_relatedness) - slope * mean(ibd_combined$distance))

# Create scatter plot with regression line
ibd_graph <- ggplot(ibd_combined, aes(x = distance, y = log_relatedness)) +
  geom_point() +
  ggtitle("Isolation By Distance") +
  geom_abline(slope = slope, intercept = intercept, color = '#376795')



# Fit a linear model

l_model <- lm(relatedness ~ distance, data = ibd_combined)

# Generate Q-Q plot for residuals
qq_plot <- ggplot(data = data.frame(resid = residuals(l_model)), aes(sample = resid)) +
  stat_qq() +
  stat_qq_line(linetype = 'dashed', color = 'firebrick', size = 1) +
  labs(
    title = "Isolation by Distance in RMBL Marmot Colonies",
    subtitle = "Residual QQ Plot"
  ) +
  theme_minimal()

# Combine scatter plot and Q-Q plot using patchwork
library(patchwork)
car_prices_graph + qq_plot


qq_plot ## looks like a skewed distribution -> using spearman's

############################# Creating PCA (modifid code from precept) ###########################
setwd("C:/Users/sophi/OneDrive/Documents/EEB331") #set your wd

#### Read in marmot metadata for pedigree construction, newest metadata with 236 indidividuals ####
metadata <- read_excel("data/radseq_pedigree_metadata.xlsx")
View(metadata) # take a look

#### Now read in our newly filtered PCA files ####
pve = read.table("data/pve.txt", header=F)
head(pve) # record % variance explained on PC1 and PC2, control F to add updated values to all PCA plots!
pcs = read.table("data/pcs.txt", header=T)

# how many individuals (rows) do we have now?
dim(pcs) # 221 marmots

############ data cleaning ########################
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
## both rivermound and southmound are counted under river
metadata_cleaned <- metadata %>%
  mutate(col_area = gsub("_(lower|upper|middle)$", "", col_area)) %>%
  mutate(col_area = gsub("_(aspen|maintalus)$", "", col_area)) %>%
  mutate(col_area = gsub("_(rivermound|sagemound|southmound)$", "", col_area)) %>%
  mutate(col_area = gsub("cliff_(lower|upper|middle)$", "mm", col_area)) %>%
  mutate(col_area = gsub("cliff", "mm", col_area))

View(metadata_cleaned)

# now MERGE new PCA data with new metadata so we can easily make updated plots!
merged_data <- merge(pcs, metadata_cleaned, by.x = "uid", by.y = "uid") |>
  select(uid, PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10, col_area) |>
  drop_na()

View(merged_data)

## CHECK number of colonies - considering 10 distinct colonies
merged_data |>
  select(col_area) |>
  distinct() ## 10 distinct colonies

## plotting number of individuals per colony
m_distribution <- marm_meta_cleaned|>
  group_by(col_area)|>
  mutate(individuals = n()) |>
  ungroup() |>
  select(col_area,individuals) |>
  ggplot(aes(x = col_area, fill = col_area)) +
  geom_bar( color = "black") +
  labs(title = "Individuals by Colony", subtitle = "Metadata from 218 Individuals", 
       x = "Collection Area", y = "Number of Marmots in Sample") +
  eeb_theme() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 7))

## look at variation in the number of individuals from each colony
m_distribution


###### anova by colony ##############
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
summary(col_aov_pc1)  # F = 67.81, p = <2e-16 ***
col_aov_pc2 <- aov(merged_data$PC2 ~ as.factor(merged_data$col_area))  
summary(col_aov_pc2)  # F = 72.43, p = <2e-16 ***       

## tukey ad-hoc for PC1
tt_1 <- tidy(TukeyHSD(col_aov_pc1))

str(tt_1) ## 45 total comparisons
nrow(tt_1 |>
       filter(adj.p.value < 0.05)) ## 22 comparisons are significant

tt_1 |>
  ## look only a significant values
  filter(adj.p.value < 0.05) |>
  select(contrast, adj.p.value) |>
  arrange(adj.p.value) |>
  slice(1:10) |>
  ## format into a nice table to look at the most significant p-values
  gt() |>
  tab_footnote(footnote = "RMBL Marmot Dataset including 218 Individuals") |>
  tab_header(
    title = "Most Significant Colony-Colony Comparisons", subtitle = 
      "Tukey Ad-Hoc of PC1") |>
  tab_style(
    style = list(
      cell_fill(color = "grey")),
    locations = cells_body(columns = contrast)
  ) |>
  
  tab_style(
    style = list(
      cell_fill(color = "lightgreen")),
    locations = cells_body(columns = adj.p.value)
  )
  
 
## picnic-bench: 0 *** ***
## picnic-gothictown: 0 *** ***
## northpk-mm: 0 *** ***
## picnic-mm: 0 *** ***
## river-picnic: 0 *** ***


## Tukey test for PC2
tt_2 <- tidy(TukeyHSD(col_aov_pc2))

## will come back later and add with comparisons are up vs. down and which ones aren't
up_valley <- c("mm", "picnic", "northpk", "boulder")

## visualize for PC2
tt_2 |>
  ## look only a significant values
  filter(adj.p.value < 0.05) |>
  select(contrast, adj.p.value) |>
  arrange(adj.p.value) |>
  slice(1:10) |>
  ## format into a nice table to look at the most significant p-values
  gt() |>
  tab_footnote(footnote = "RMBL Marmot Dataset including 218 Individuals") |>
  tab_header(
    title = "Most Significant Colony-Colony Comparisons", subtitle = 
      "Tukey Ad-Hoc of PC2") |>
  tab_style(
    style = list(
      cell_fill(color = "grey")),
    locations = cells_body(columns = contrast)
  ) |>
  
  tab_style(
    style = list(
      cell_fill(color = "lightblue")),
    locations = cells_body(columns = adj.p.value)
  )

## gothictown - bench: 0 *** ***
## mm - bench: 0 *** ***
## picnic_bench: 0 *** ***
## river-mm: 0 *** ***
## northpk-mm: 0 *** ***

plot(tt_2)

## mantel testing


## testing different relatedness metrics