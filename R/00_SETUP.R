### 01. PACKAGES ####
suppressPackageStartupMessages({
library(janitor)
library(lubridate)
library(readxl)
library(glmmTMB)
library(DHARMa)
library(performance)
library(emmeans)
library(vegan)
library(patchwork)
library(tidyverse)
library(mvabund)
library(dplyr)
library(beepr)
})
set.seed(42)


### 02. DIRECTORIES ####

analysis_date <- format(Sys.Date(), "%Y.%m.%d")

data_raw_dir <- "data/raw"
data_clean_dir <- "data/clean"
docs_dir <- "docs"

output_dir <- file.path("outputs", paste0("Analysis_", analysis_date))
fits_dir <- file.path(output_dir, "fits")
plots_dir <- file.path(output_dir, "plots")
tables_dir <- file.path(output_dir, "tables")
stats_dir <- file.path(output_dir, "stats")

dirs <- c(data_raw_dir, data_clean_dir, docs_dir,
          fits_dir, plots_dir, tables_dir, stats_dir)

for (dir in dirs) dir.create(dir, recursive = TRUE, showWarnings = FALSE)
rm(dirs, dir)


### 03. PLOT SETTINGS ####

#### Plot theme ####

theme_clean <- theme_minimal(base_family = "sans") +
  theme(
    legend.position = "bottom",
    panel.grid = element_blank(),
    plot.title = element_blank(),
    panel.background = element_rect(fill = "transparent", colour = NA),
    plot.background = element_rect(fill = "transparent", colour = NA)
  )

theme_set(theme_clean)


#### Colour palettes ####

stability_cols <- c(
  "Mobile"      = "#B9C7B5",
  "Semi-stable" = "#477F91",
  "Stable"      = "#1F4E6B"
)

degradation_cols <- c(
  "Low"      = "#E8D8B5",
  "Moderate" = "#B9C7B5",
  "High"     = "#477F91",
  "Severe"   = "#1F4E6B"
)

coral_cols <- c(
  "Recruits"     = "#1F4E6B",
  "Other corals" = "#7FA9A8"
)

benthic_cols <- c(
  "Macroalgae"    = "#477F91",
  "Hard coral"   = "#1F4E6B",
  "CCA"          = "#B9C7B5",
  "Sponge"       = "#6372EB",
  "Turf algae"   = "#7FA9A8",
  "Soft coral"   = "#90766C",
  "Zoanthids"    = "#D1C8B5",
  "Corallimorphs" = "#8C819D",
  "Tunicates"    = "#454545",
  "Sediment"     = "#E8D8B5",
  "Unknown"      = "#B0B0B0"
)

action_cols_traffic <- c(
  "Remove"           = "#C45C54",
  "Consider Removal" = "#E8C46C",
  "Retain"           = "#53866D"
)


action_cols_accessible <- c(
  "Remove"           = "#E2A76F",
  "Consider Removal" = "#D1C8B5",
  "Retain"           = "#1F4E6B"
)

heatmap_cols <- c("#F5F3EE", "#B9C7B5", "#477F91", "#1F4E6B")

### 04. SESSION INFO ####

# sessionInfo()

