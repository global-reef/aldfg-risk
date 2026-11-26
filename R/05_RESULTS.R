### 01. SAMPLE CHARACTERISTICS ####

results_dir <- file.path(tables_dir, "manuscript")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

sample_summary <- dfg %>%
  summarise(
    n_dfg = n(),
    n_sites = n_distinct(site),
    n_images = sum(n_images),
    n_points = sum(n_points),
    total_area_m2 = sum(surface_area_m2),
    across(c(surface_area_m2, depth_m, coral_abundance, coral_recruits,
             genus_richness, growth_form_richness, fauna_abundance,
             benthic_richness, total_biotic_cover, all_of(biotic_cols),
             cover_sediment),
           list(mean = mean, sd = sd, median = median, min = min, max = max))
  ) %>%
  pivot_longer(everything(), names_to = "statistic", values_to = "value")

sample_summary %>% print(n = Inf)

dfg %>% count(dfg_type)
dfg %>% count(stability, degradation)
dfg %>% count(substrate_group)


#### B. Benthic category occurrence and cover ####

# Occurrence is detection in photoquadrats; positive cover summaries exclude zeros.
cover_long <- dfg %>%
  select(dfg_id, all_of(biotic_cols)) %>%
  pivot_longer(all_of(biotic_cols), names_to = "category", values_to = "cover") %>%
  mutate(category = recode(category, !!!cover_labels))

cover_summary <- cover_long %>%
  group_by(category) %>%
  summarise(n_total = n(), n_present = sum(cover > 0),
            occurrence = 100 * n_present / n_total,
            mean_cover = mean(cover), median_present = median(cover[cover > 0]),
            q25_present = quantile(cover[cover > 0], 0.25),
            q75_present = quantile(cover[cover > 0], 0.75),
            max_cover = max(cover), .groups = "drop") %>%
  arrange(desc(n_present))

print(cover_summary, n = Inf, width = Inf)

### 02. PRIMARY MODEL RESULTS ####

#### A. Overall predictor tests ####

extract_tests <- function(model) {
  drop1(model, test = "Chisq") %>%
    as.data.frame() %>%
    rownames_to_column("predictor") %>%
    filter(predictor != "<none>") %>%
    transmute(predictor, df = Df, statistic = LRT, p = `Pr(>Chi)`) %>%
    as_tibble()
}

primary_tests <- map_dfr(
  list(O1A_biotic_cover = m1a, O1C_benthic_richness = m1c,
       O2B_recruits = m2b, O2C_genus_richness = m2c),
  extract_tests, .id = "model"
)

primary_tests %>% print(n = Inf)


#### B. O1B: Community composition ####

m1b_test
RsquareAdj(m1b)
nmds1b$stress

primary_tests <- bind_rows(
  primary_tests,
  as.data.frame(m1b_test) %>%
    rownames_to_column("predictor") %>%
    filter(!predictor %in% c("Residual", "Total")) %>%
    transmute(model = "O1B_composition", predictor,
              df = Df, statistic = F, p = `Pr(>F)`)
)


#### C. O2A: Coral and associated fauna abundance ####

m2a_stability
m2a_degradation
m2a_sediment
m2a_area

# Retain joint tests and adjusted response-specific inference.
writeLines(
  capture.output({
    print(m2a_stability)
    print(m2a_degradation)
    print(m2a_sediment)
    print(m2a_area)
  }),
  file.path(results_dir, "O2A_multivariate_tests.txt")
)


### 03. KEY EFFECT SIZES AND POST-HOC COMPARISONS ####

#### A. Continuous predictor effects ####

extract_effects <- function(model, terms) {
  summary(model)$coefficients$cond %>%
    as.data.frame() %>%
    rownames_to_column("term") %>%
    filter(term %in% terms) %>%
    transmute(
      term, estimate = Estimate, SE = `Std. Error`,
      ratio = exp(estimate),
      lower_CL = exp(estimate - qnorm(0.975) * SE),
      upper_CL = exp(estimate + qnorm(0.975) * SE)
    )
}

key_effect_sizes <- bind_rows(
  extract_effects(m1c, c("cover_sediment", "log_effort")) %>%
    mutate(model = "O1C_benthic_richness"),
  extract_effects(m2b, "log_area") %>%
    mutate(model = "O2B_recruits"),
  extract_effects(m2c, "log_area") %>%
    mutate(model = "O2C_genus_richness")
) %>%
  relocate(model)

key_effect_sizes


#### B. Estimated marginal means ####

# Stability effects supported by the overall predictor tests.
emm2b_stability <- emmeans(m2b, ~ stability, type = "response")
emm2c_stability <- emmeans(m2c, ~ stability, type = "response")

summary(emm2b_stability, infer = c(TRUE, FALSE))
summary(emm2c_stability, infer = c(TRUE, FALSE))


#### C. Pairwise stability comparisons ####

pairs2b_stability <- pairs(emm2b_stability, adjust = "tukey")
pairs2c_stability <- pairs(emm2c_stability, adjust = "tukey")

pairwise_comparisons <- bind_rows(
  as.data.frame(summary(pairs2b_stability, infer = TRUE, type = "response")) %>%
    mutate(model = "O2B_recruits"),
  as.data.frame(summary(pairs2c_stability, infer = TRUE, type = "response")) %>%
    mutate(model = "O2C_genus_richness")
) %>%
  relocate(model)

pairwise_comparisons


### 04. OBJECTIVE 3: REMOVAL FRAMEWORK ####

#### A. Risk and ecological value ####

framework %>%
  summarise(
    mean_risk = mean(risk_scale),
    mean_value = mean(ecological_value),
    min_value = min(ecological_value), max_value = max(ecological_value),
    sd_value = sd(ecological_value),
    percent_zero = 100 * mean(ecological_value == 0)
  )

risk_ecology
lmtest::coeftest(m3, vcov. = m3_vcov)
lmtest::coefci(m3, vcov. = m3_vcov)

# Expected ecological scores and robust 95% confidence intervals.
risk_predictions <- tibble(risk_scale = seq(1, 5, length.out = 100))

risk_design <- model.matrix(~ risk_scale, data = risk_predictions)
risk_link <- predict(m3, newdata = risk_predictions, type = "link")
risk_se <- sqrt(diag(risk_design %*% m3_vcov %*% t(risk_design)))

risk_predictions <- risk_predictions %>%
  mutate(
    predicted = 8 * plogis(risk_link),
    lower_CL = 8 * plogis(risk_link - qnorm(0.975) * risk_se),
    upper_CL = 8 * plogis(risk_link + qnorm(0.975) * risk_se)
  )

rm(risk_design, risk_link, risk_se)


#### B. Management recommendations ####

action_summary <- framework %>%
  count(action, .drop = FALSE) %>%
  mutate(percent = 100 * n / sum(n))

action_summary

table(framework$risk_scale, framework$action)
table(framework$ecological_value, framework$action)


#### C. Illustrative management cases ####

example_ids <- paste("Net", c(86, 13, 45, 90, 20, 21))

framework_examples <- framework %>%
  filter(dfg_id %in% example_ids) %>%
  select(
    dfg_id, site, dfg_type, depth_m, surface_area_m2,
    stability, degradation, risk_scale, ecological_value,
    coral_abundance, coral_recruits, genus_richness,
    growth_form_richness, combined_benthic_richness, fauna_abundance, fauna_notes,
    starts_with("cover_"), total_biotic_cover, action
  ) %>%
  arrange(action, desc(risk_scale))

stopifnot(nrow(framework_examples) == 6)

framework_examples %>% print(n = Inf, width = Inf)


### 05. SENSITIVITY ANALYSES ####

#### A. Site effects and growth-form richness ####

sensitivity_summary <- map_dfr(
  list(O1A_no_site = m1a_beta, O1C_no_site = m1c_comp,
       O2A_coral_site = m2a_coral_site, O2A_fauna_site = m2a_fauna_site,
       O2B_site = m2b_site, O2C_site = m2c_site,
       O2C_growth_form = m2c_gf),
  extract_tests, .id = "model"
) %>%
  select(model, predictor, p)

sensitivity_summary %>% print(n = Inf)

m1b_unrestricted
m1b_restricted


#### B. Earlier coral-only framework and stability bonus ####

risk_ecology_current
lmtest::coeftest(m3_current, vcov. = m3_current_vcov)

framework_comparison <- framework %>%
  count(action_current, action, .drop = FALSE)
framework_comparison

risk_ecology_bonus
lmtest::coeftest(m3_bonus, vcov. = m3_bonus_vcov)

lmtest::coeftest(m3, vcov. = m3_site_vcov)


### 06. SAVE MANUSCRIPT RESULTS ####

write_csv(cover_summary, file.path(results_dir, "benthic_occurrence_cover.csv"))
write_csv(framework_comparison, file.path(results_dir, "framework_comparison.csv"))
write_csv(sample_summary, file.path(results_dir, "sample_summary.csv"))
write_csv(primary_tests, file.path(results_dir, "primary_tests.csv"))
write_csv(key_effect_sizes, file.path(results_dir, "key_effect_sizes.csv"))
write_csv(pairwise_comparisons, file.path(results_dir, "pairwise_comparisons.csv"))
write_csv(sensitivity_summary, file.path(results_dir, "sensitivity_summary.csv"))
write_csv(action_summary, file.path(results_dir, "framework_summary.csv"))
write_csv(framework_examples, file.path(results_dir, "framework_examples.csv"))

# Supplementary ALDFG-level management table.
framework %>%
  select(
    dfg_id, site, dfg_type, depth_m, surface_area_m2,
    stability, degradation, risk_scale, coral_abundance,
    coral_recruits, genus_richness, growth_form_richness,
    fauna_abundance, cover_hard_coral, combined_benthic_richness,
    score_recruits, score_hc_cover, score_combined_richness, score_fauna,
    ecological_value, action, ecological_value_current, action_current
  ) %>%
  arrange(action, desc(risk_scale), ecological_value) %>%
  write_csv(file.path(results_dir, "framework_gear.csv"))

# Retain estimated means and predictions for figure preparation.
saveRDS(
  list(
    emm2b_stability = emm2b_stability,
    emm2c_stability = emm2c_stability,
    risk_predictions = risk_predictions,
    action_summary = action_summary,
    framework_examples = framework_examples
  ),
  file.path(fits_dir, "manuscript_results.rds")
)

rm(extract_tests, extract_effects)

beep()

# other small results 
# range of "other substrates" 
dfg %>%
  group_by(stability) %>%
  summarise(mean = mean(total_biotic_cover),
            min = min(total_biotic_cover),
            max = max(total_biotic_cover))

dfg %>%
  group_by(degradation) %>%
  summarise(mean = mean(total_biotic_cover),
            min = min(total_biotic_cover),
            max = max(total_biotic_cover))