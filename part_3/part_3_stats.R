# AUTHOR / EMAIL: Aya Surheyao (as4789@princeton.edu)
# CREATION DATE: 04/06/2026
# MODIFIED DATE: 04/20/2026
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
library(MetBrewer)

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

qq_plot ## looks like a skewed distribution -> using spearman's


## new qq_plot using log relatedness
# Fit a linear model
model2 <- lm(relatedness ~ distance, data = ibd_up)

# Generate Q-Q plot for residuals
up_qq_plot <- ggplot(data = data.frame(resid = residuals(model2)), aes(sample = resid)) +
  stat_qq() +
  stat_qq_line(linetype = 'dashed', color = 'seagreen', size = 1) +
  labs(
    title = "Isolation by Distance in RMBL Marmot Colonies",
    subtitle = "Residual QQ Plot"
  ) +
  theme_minimal()

up_qq_plot


# Fit a linear model
model3 <- lm(relatedness ~ distance, data = ibd_down)

# Generate Q-Q plot for residuals
down_qq_plot <- ggplot(data = data.frame(resid = residuals(model3)), aes(sample = resid)) +
  stat_qq() +
  stat_qq_line(linetype = 'dashed', color = 'slateblue', size = 1) +
  labs(
    title = "Isolation by Distance in RMBL Marmot Colonies",
    subtitle = "Residual QQ Plot"
  ) +
  theme_minimal()

down_qq_plot


# Fit a linear model
model4 <- lm(relatedness ~ distance, data = ibd_diff_colony)

# Generate Q-Q plot for residuals
diff_qq_plot <- ggplot(data = data.frame(resid = residuals(model4)), aes(sample = resid)) +
  stat_qq() +
  stat_qq_line(linetype = 'dashed', color = 'salmon', size = 1) +
  labs(
    title = "Isolation by Distance in RMBL Marmot Colonies",
    subtitle = "Residual QQ Plot"
  ) +
  theme_minimal()

diff_qq_plot

## tukey ad-hoc for PC1
tt_1 <- tidy(TukeyHSD(col_aov_pc1))


str(tt_1) ## 45 total comparisons
colnames(tt_1)
nrow(tt_1 |>
       filter(adj.p.value < 0.001)) ## 22 comparisons are significant

ttable_1 <- tt_1 |>
  ## look only a significant values
  filter(adj.p.value < 0.001) |>
  select(contrast, adj.p.value) |>
  arrange(adj.p.value) |>
  ## format into a nice table to look at the most significant p-values
  gt() |>
  tab_footnote(footnote = "RMBL Marmot Dataset including 218 Individuals") |>
  tab_header(
    title = "Significant Colony-Colony Differences", subtitle = 
      "Post-Hoc Analysis of PC1 (p < 0.001)") |>
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

nrow(tt_2 |>
       filter(adj.p.value < 0.001))

## visualize for PC2
ttable_2 <- tt_2 |>
  ## look only a significant values
  filter(adj.p.value < 0.001) |>
  select(contrast, adj.p.value) |>
  arrange(adj.p.value) |>
  ## format into a nice table to look at the most significant p-values
  gt() |>
  tab_footnote(footnote = "RMBL Marmot Dataset including 218 Individuals") |>
  tab_header(
    title = "Significant Colony-Colony Differences", subtitle = 
      "Post-Hoc Analysis of PC2 (p < 0.001)") |>
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

gtsave(ttable_1, "tukey1.png")
gtsave(ttable_2, "tukey2.png")

## gothictown - bench: 0 *** ***
## mm - bench: 0 *** ***
## picnic_bench: 0 *** ***
## river-mm: 0 *** ***
## northpk-mm: 0 *** ***

plot(tt_2)

## mantel testing


## testing different relatedness metrics