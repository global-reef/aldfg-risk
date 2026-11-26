### 01. SAMPLING DESIGN ####

dfg %>% count(site, sort = TRUE) %>% print(n = Inf)
dfg %>% count(dfg_group)
dfg %>% count(substrate_group)

table(dfg$stability, dfg$degradation)
table(dfg$site, dfg$stability)
table(dfg$site, dfg$degradation)
table(dfg$dfg_group, dfg$stability)
table(dfg$substrate_group, dfg$stability)

dfg %>%
  group_by(stability) %>%
  summarise(n = n(), n_sites = n_distinct(site),
            median_area = median(surface_area_m2),
            min_area = min(surface_area_m2),
            max_area = max(surface_area_m2))

dfg %>%
  group_by(degradation) %>%
  summarise(n = n(), n_sites = n_distinct(site),
            median_area = median(surface_area_m2),
            min_area = min(surface_area_m2),
            max_area = max(surface_area_m2))


### 02. DERIVE ECOLOGICAL METRICS ####

# Unknown classifications are excluded from biological group richness.

dfg <- dfg %>%
  mutate(
    total_biotic_cover = rowSums(pick(all_of(biotic_cols))),
    benthic_richness = rowSums(pick(all_of(biotic_cols)) > 0),
    recruit_density = coral_recruits / surface_area_m2,
    fauna_density = fauna_abundance / surface_area_m2
  )

count_vars <- c("coral_recruits", "genus_richness", "growth_form_richness",
                "fauna_abundance", "benthic_richness")

continuous_vars <- c("surface_area_m2", "depth_m", "mesh_diameter",
                     "cover_sediment", "n_images", "total_biotic_cover")


### 03. VARIABLE DISTRIBUTIONS ####

dfg %>%
  select(all_of(c(count_vars, continuous_vars))) %>%
  pivot_longer(everything(), names_to = "variable") %>%
  group_by(variable) %>%
  summarise(
    n = sum(!is.na(value)), missing = sum(is.na(value)),
    zeros = sum(value == 0, na.rm = TRUE),
    pct_zero = 100 * zeros / n,
    min = min(value, na.rm = TRUE),
    median = median(value, na.rm = TRUE),
    mean = mean(value, na.rm = TRUE),
    max = max(value, na.rm = TRUE),
    variance = var(value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(across(where(is.numeric), ~round(., 3))) %>%
  print(n = Inf, width = Inf)

p_distributions <- dfg %>%
  select(all_of(c(count_vars, continuous_vars))) %>%
  pivot_longer(everything(), names_to = "variable") %>%
  ggplot(aes(value)) +
  geom_histogram(bins = 20, fill = "#4393C3", colour = "white") +
  facet_wrap(~variable, scales = "free", ncol = 3) +
  labs(x = NULL, y = "Frequency")

p_distributions


### 04. GEAR CHARACTERISTICS AND SAMPLING EFFORT ####

p_area <- (
  ggplot(dfg, aes(stability, surface_area_m2, fill = stability)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.15, size = 1.5) +
    scale_y_log10() +
    scale_fill_manual(values = stability_cols) +
    labs(x = "Structural stability", y = expression("ALDFG area (m"^2*")"))
) + (
  ggplot(dfg, aes(degradation, surface_area_m2, fill = degradation)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.15, size = 1.5) +
    scale_y_log10() +
    scale_fill_manual(values = degradation_cols) +
    labs(x = "Degree of degradation", y = expression("ALDFG area (m"^2*")"))
) & theme(legend.position = "none")

p_area

p_effort <- (
  ggplot(dfg, aes(surface_area_m2, n_images)) +
    geom_point() +
    scale_x_log10() +
    labs(x = expression("ALDFG area (m"^2*")"), y = "Number of images")
) + (
  ggplot(dfg, aes(n_images, benthic_richness)) +
    geom_jitter(width = 0.15, height = 0.1) +
    labs(x = "Number of images", y = "Benthic richness")
)

p_effort

dfg %>%
  select(surface_area_m2, depth_m, mesh_diameter,
         cover_sediment, n_images, benthic_richness) %>%
  cor(method = "spearman", use = "pairwise.complete.obs") %>%
  round(2)


### 05. ECOLOGICAL RESPONSES ####

p_stability <- dfg %>%
  select(stability, all_of(count_vars)) %>%
  pivot_longer(-stability, names_to = "response") %>%
  ggplot(aes(stability, value, fill = stability)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.6) +
  scale_fill_manual(values = stability_cols) +
  facet_wrap(~response, scales = "free_y") +
  labs(x = "Structural stability", y = NULL) +
  theme(legend.position = "none")

p_stability

p_degradation <- dfg %>%
  select(degradation, all_of(count_vars)) %>%
  pivot_longer(-degradation, names_to = "response") %>%
  ggplot(aes(degradation, value, fill = degradation)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.6) +
  scale_fill_manual(values = degradation_cols) +
  facet_wrap(~response, scales = "free_y") +
  labs(x = "Degree of degradation", y = NULL) +
  theme(legend.position = "none")

p_degradation

p_sediment <- dfg %>%
  select(cover_sediment, all_of(count_vars)) %>%
  pivot_longer(-cover_sediment, names_to = "response") %>%
  ggplot(aes(cover_sediment, value)) +
  geom_point(alpha = 0.6) +
  facet_wrap(~response, scales = "free_y") +
  labs(x = "Sediment cover (%)", y = NULL)

p_sediment

p_counts_area <- dfg %>%
  select(surface_area_m2, stability, all_of(count_vars)) %>%
  pivot_longer(all_of(count_vars), names_to = "response") %>%
  ggplot(aes(surface_area_m2, value, colour = stability)) +
  geom_point(alpha = 0.7) +
  scale_x_log10() +
  scale_colour_manual(values = stability_cols) +
  facet_wrap(~response, scales = "free_y") +
  labs(x = expression("ALDFG area (m"^2*")"), y = NULL, colour = "Stability")

p_counts_area

dfg %>%
  select(surface_area_m2, n_images, cover_sediment, all_of(count_vars)) %>%
  cor(method = "spearman", use = "pairwise.complete.obs") %>%
  round(2)


### 06. BENTHIC COMMUNITY COMPOSITION ####

dfg %>%
  select(all_of(biotic_cols), cover_sediment, cover_unknown) %>%
  pivot_longer(everything(), names_to = "category", values_to = "cover") %>%
  group_by(category) %>%
  summarise(
    mean_cover = mean(cover),
    median_cover = median(cover),
    max_cover = max(cover),
    prevalence = 100 * mean(cover > 0),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_cover)) %>%
  mutate(across(where(is.numeric), ~round(., 2)))

# All-zero biological communities require consideration before Bray-Curtis analysis.
dfg %>%
  filter(total_biotic_cover == 0) %>%
  select(dfg_id, stability, degradation, n_points, cover_sediment, cover_unknown)

p_composition <- dfg %>%
  select(stability, degradation, all_of(biotic_cols), cover_sediment, cover_unknown) %>%
  pivot_longer(starts_with("cover_"), names_to = "category", values_to = "cover") %>%
  group_by(stability, degradation, category) %>%
  summarise(mean_cover = mean(cover), .groups = "drop") %>%
  ggplot(aes(degradation, mean_cover, fill = category)) +
  geom_col() +
  facet_wrap(~stability) +
  labs(x = "Degree of degradation", y = "Mean cover (%)", fill = "Category")

p_composition

p_biotic <- ggplot(dfg, aes(cover_sediment, total_biotic_cover,
                            colour = stability)) +
  geom_point() +
  geom_abline(intercept = 100, slope = -1, linetype = "dashed") +
  scale_colour_manual(values = stability_cols) +
  labs(x = "Sediment cover (%)", y = "Total biotic cover (%)",
       colour = "Stability")

p_biotic


### 07. SAVE EXPLORATORY OUTPUTS ####

explore_dir <- file.path(plots_dir, "exploratory")
dir.create(explore_dir, recursive = TRUE, showWarnings = FALSE)

ggsave(file.path(explore_dir, "distributions.png"), p_distributions, width = 10, height = 9, dpi = 300)
ggsave(file.path(explore_dir, "gear_area.png"), p_area, width = 9, height = 4, dpi = 300)
ggsave(file.path(explore_dir, "sampling_effort.png"), p_effort, width = 9, height = 4, dpi = 300)
ggsave(file.path(explore_dir, "responses_stability.png"), p_stability, width = 10, height = 6, dpi = 300)
ggsave(file.path(explore_dir, "responses_degradation.png"), p_degradation, width = 10, height = 6, dpi = 300)
ggsave(file.path(explore_dir, "responses_sediment.png"), p_sediment, width = 10, height = 6, dpi = 300)
ggsave(file.path(explore_dir, "responses_area.png"), p_counts_area, width = 10, height = 6, dpi = 300)
ggsave(file.path(explore_dir, "benthic_composition.png"), p_composition, width = 10, height = 5, dpi = 300)
ggsave(file.path(explore_dir, "biotic_sediment.png"), p_biotic, width = 6, height = 5, dpi = 300)

beep()