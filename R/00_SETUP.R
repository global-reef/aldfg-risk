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
  "Mobile"      = "#8FA888",
  "Semi-stable" = "#477F91",
  "Stable"      = "#1F4E6B"
)
stability_shapes <- c(
  "Mobile"      = 17,
  "Semi-stable" = 16,
  "Stable"      = 15
)

degradation_cols <- c(
  "Low"      = "#E8D8B5",
  "Moderate" = "#8FA888",
  "High"     = "#477F91",
  "Severe"   = "#1F4E6B"
)

coral_cols <- c(
  "Recruits"     = "#1F4E6B",
  "Other corals" = "#7FA9A8"
)

benthic_cols <- c(
  "Macroalgae"    = "#B9C7B5",
  "Hard coral"   = "#477F91",
  "CCA"          = "#8C819D",
  "Sponge"       = "#7888C2",
  "Turf algae"   = "#7FA9A8",
  "Soft coral"   = "#1F4E6B",
  "Zoanthids"    = "#D1C8B5",
  "Corallimorphs" = "#B58A78",
  "Tunicates"    = "#454545",
  "Sediment"     = "#E8D8B5",
  "Unknown"      = "#B0B0B0"
)

action_cols <- c(
  "Prioritise removal"   = "#E2A76F",
  "Monitor"              = "#D1C8B5",
  "Low removal priority" = "#1F4E6B"
)

heatmap_cols <- c("#F5F3EE", "#B9C7B5", "#477F91", "#1F4E6B")

### 04. SESSION INFO AND CITATIONS ####

#### Session information ####

sessionInfo()
writeLines(capture.output(sessionInfo()), file.path(stats_dir, "session_info.txt"))

#### R and package citations ####

citation_packages <- c(
  "janitor", "lubridate", "readxl", "glmmTMB", "DHARMa",
  "performance", "emmeans", "vegan", "patchwork", "tidyverse",
  "mvabund", "dplyr", "sandwich", "lmtest"
)

writeLines(capture.output({
  print(citation(), style = "text")
  
  for (pkg in citation_packages) {
    cat("\n\n###", pkg, "- version", as.character(packageVersion(pkg)), "###\n")
    print(citation(pkg), style = "text")
  }
}), file.path(stats_dir, "citations.txt"))

rm(citation_packages, pkg)
