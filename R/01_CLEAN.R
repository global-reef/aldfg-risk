### 01. IMPORT DATA ####

dfg_raw <- read_csv(file.path(data_raw_dir, "raw_net_datasheet.csv"), show_col_types = FALSE)
benthos_raw <- read_csv(file.path(data_raw_dir, "raw_datapoint.csv"), show_col_types = FALSE)


### 02. CLEAN ALDFG DATA ####

# Exclude Net 5, classified as unusable in the original field datasheet.
dfg <- dfg_raw %>%
  clean_names() %>%
  mutate(across(where(is.character), str_squish)) %>%
  filter(!is.na(fishinggear_id), fishinggear_id != "Net 5") %>%
  transmute(
    dfg_id = fishinggear_id,
    date = dmy(if_else(date == "27203/2026", "27/03/2026", date)),
    site = dive_site,
    depth_m,
    substrate,
    length_m = lenght_m,
    width_m,
    surface_area_m2,
    mesh_diameter,
    stability = structural_stability,
    degradation = degree_of_degradation,
    dfg_type = gear_type,
    coral_abundance = total_number_of_corals,
    coral_recruits,
    genus_richness = species_richness,
    growth_form_richness = growth_form_diversity,
    living_tissue,
    partial_mortality,
    dead,
    sign_of_disease,
    bleaching,
    fauna_abundance = parse_number(cryptic_organisms),
    fauna_notes = cryptic_organisms
  ) %>%
  mutate(
    site = recode(site, "Red rock" = "Red Rock",
                  "Leam Thiaen South" = "Laem Thiaen South",
                  "Japanese gardens" = "Japanese Garden"),
    substrate = recode(substrate, "Corals" = "Coral",
                       "hard corals" = "Coral",
                       "Hard coral/Rock" = "Coral/Rock",
                       "Rock/ Hard coral" = "Coral/Rock"),
    dfg_type = str_to_lower(dfg_type),
    dfg_type = recode(dfg_type, "rope/wire" = "wire/rope",
                      "rope/wheel" = "rope/tire", "wheel/rope" = "rope/tire",
                      "tire (wheel)" = "tire", "cloths" = "cloth"),
    stability = factor(stability, levels = c("Mobile", "Semi-stable", "Stable")),
    degradation = factor(degradation, levels = c("Low", "Moderate", "High", "Severe")),
    dfg_group = case_when(
      dfg_type %in% c("wire", "steel", "metal") ~ "Metal",
      dfg_type %in% c("plastic", "rope", "nylon", "fishing line") ~ "Plastic",
      str_detect(dfg_type, "/") ~ "Mixed",
      dfg_type %in% c("cloth", "tire") ~ "Other"
    ),
    substrate_group = case_when(
      substrate %in% c("Coral", "Coral/Rock") ~ "Biotic",
      substrate %in% c("Rock", "Rock/sand", "Sand") ~ "Abiotic"
    )
  )


### 03. CLEAN REEFCLOUD DATA ####

benthos <- benthos_raw %>%
  clean_names() %>%
  transmute(
    dfg_id = str_squish(survey_title),
    survey_id, image_id, point_id, point_num, image_disabled,
    benthic_group = point_human_group_code,
    benthic_code = point_human_classification
  ) %>%
  filter(dfg_id %in% dfg$dfg_id, !image_disabled) %>%
  mutate(benthic_group = recode(benthic_group, "Coralline Algae" = "CCA",
                                "Macro" = "Macroalgae", "Turf Algae" = "Turf",
                                "Other" = "Unknown")) %>%
  select(-image_disabled)

sampling_effort <- benthos %>%
  group_by(dfg_id) %>%
  summarise(n_images = n_distinct(image_id), n_points = n(),
            n_unknown = sum(benthic_group == "Unknown"),
            n_identified = n_points - n_unknown, .groups = "drop")

# Cover is calculated relative to all annotation points, including unknown classifications.
cover <- benthos %>%
  count(dfg_id, benthic_group, name = "n") %>%
  group_by(dfg_id) %>%
  mutate(cover = 100 * n / sum(n)) %>%
  ungroup() %>%
  select(-n) %>%
  pivot_wider(names_from = benthic_group, values_from = cover,
              names_prefix = "cover_", values_fill = 0) %>%
  clean_names()


### 04. BUILD MASTER DATASET ####

dfg <- dfg %>%
  left_join(sampling_effort, by = "dfg_id") %>%
  left_join(cover, by = "dfg_id")

rm(sampling_effort, cover)


### 05. VALIDATE CLEANED DATA ####

stopifnot(
  nrow(dfg) == 106,
  !anyDuplicated(dfg$dfg_id),
  !anyDuplicated(benthos$point_id),
  setequal(dfg$dfg_id, benthos$dfg_id),
  !anyNA(dfg$date),
  !anyNA(dfg$stability),
  !anyNA(dfg$degradation),
  !anyNA(dfg$dfg_group),
  !anyNA(dfg$substrate_group),
  !anyNA(dfg$n_points),
  !anyNA(benthos$benthic_group),
  all(dfg$surface_area_m2 > 0),
  all(abs(dfg$surface_area_m2 - dfg$length_m * dfg$width_m) < 0.00001),
  all(abs(rowSums(select(dfg, starts_with("cover_"))) - 100) < 0.00001),
  all(dfg$n_points == dfg$n_images * 25),
  all(dfg$n_points == dfg$n_identified + dfg$n_unknown),
  all(dfg$coral_recruits <= dfg$coral_abundance)
)

glimpse(dfg)
glimpse(benthos)

dfg %>% count(stability)
dfg %>% count(degradation)
dfg %>% count(dfg_type)
dfg %>% count(dfg_group)
dfg %>% count(substrate)


### 06. SAVE CLEANED DATA ####

write_csv(dfg, file.path(data_clean_dir, "dfg_clean.csv"))
write_csv(benthos, file.path(data_clean_dir, "benthos_clean.csv"))

beep()