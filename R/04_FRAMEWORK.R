### 01. CALCULATE RISK AND ECOLOGICAL VALUE ####

framework <- dfg %>%
  mutate(
    risk_scale = case_when(
      stability == "Stable" ~ 1L,
      stability == "Semi-stable" & degradation %in% c("Low", "Moderate") ~ 2L,
      stability == "Semi-stable" & degradation %in% c("High", "Severe") ~ 3L,
      stability == "Mobile" & degradation %in% c("Low", "Moderate") ~ 4L,
      stability == "Mobile" & degradation %in% c("High", "Severe") ~ 5L
    ),
    score_recruits = case_when(coral_recruits >= 6 ~ 2L,
                               coral_recruits >= 1 ~ 1L, TRUE ~ 0L),
    score_genus = case_when(genus_richness >= 3 ~ 2L,
                            genus_richness == 2 ~ 1L, TRUE ~ 0L),
    score_growth = case_when(growth_form_richness >= 3 ~ 2L,
                             growth_form_richness == 2 ~ 1L, TRUE ~ 0L),
    score_stability = as.integer(risk_scale == 1),
    ecological_value = score_recruits + score_genus + score_growth + score_stability,
    ecological_no_bonus = ecological_value - score_stability
  )

stopifnot(
  nrow(framework) == 106,
  !anyDuplicated(framework$dfg_id),
  !anyNA(framework %>% select(coral_recruits, genus_richness, growth_form_richness)),
  all(framework$risk_scale %in% 1:5),
  all(framework$ecological_value %in% 0:7),
  all(framework$ecological_no_bonus %in% 0:6)
)

framework %>% count(risk_scale)
framework %>% count(ecological_value)

framework %>%
  summarise(
    mean_risk = mean(risk_scale),
    mean_value = mean(ecological_value),
    sd_value = sd(ecological_value),
    percent_zero = 100 * mean(ecological_value == 0)
  )


### 02. REMOVAL DECISION FRAMEWORK ####

action_levels <- c("Remove", "Consider Removal", "Retain")

classify_action <- function(risk, value) {
  factor(
    case_when(
      risk == 1 ~ "Retain",
      value <= 2 ~ "Remove",
      risk %in% c(2, 3) & value >= 6 ~ "Retain",
      TRUE ~ "Consider Removal"
    ),
    levels = action_levels
  )
}

framework <- framework %>%
  mutate(action = classify_action(risk_scale, ecological_value))

decision_matrix <- expand_grid(risk_scale = 1:5, ecological_value = 0:7) %>%
  mutate(action = classify_action(risk_scale, ecological_value)) %>%
  left_join(framework %>% count(risk_scale, ecological_value),
            by = c("risk_scale", "ecological_value")) %>%
  mutate(n = replace_na(n, 0L))

stopifnot(
  !anyNA(framework$action),
  sum(decision_matrix$n) == nrow(framework),
  all(framework$action[framework$risk_scale == 1] == "Retain"),
  !any(framework$action == "Retain" & framework$risk_scale >= 4)
)

framework %>%
  count(action, .drop = FALSE) %>%
  mutate(percent = 100 * n / sum(n))

table(framework$risk_scale, framework$action)
table(framework$ecological_value, framework$action)


### 03. RISK AND ECOLOGICAL VALUE TRADE-OFF ####

#### A. Primary analysis ####

risk_ecology <- cor.test(
  framework$risk_scale, framework$ecological_value,
  method = "spearman", exact = FALSE
)

risk_ecology

# Fractional logistic regression accommodates bounded composite scores,
# including exact zeros and ones, without assuming binomial trials.

m3 <- glm(
  I(ecological_value / 7) ~ risk_scale,
  family = quasibinomial(link = "logit"), data = framework
)

# Robust covariance accounts for variance misspecification.
m3_vcov <- sandwich::vcovHC(m3, type = "HC3")

lmtest::coeftest(m3, vcov. = m3_vcov)
lmtest::coefci(m3, vcov. = m3_vcov)


#### B. Stability bonus sensitivity ####

# Remove the component directly derived from the risk score.
# Rescale using the revised theoretical maximum of six.

risk_ecology_no_bonus <- cor.test(
  framework$risk_scale, framework$ecological_no_bonus,
  method = "spearman", exact = FALSE
)

risk_ecology_no_bonus

m3_no_bonus <- glm(
  I(ecological_no_bonus / 6) ~ risk_scale,
  family = quasibinomial(link = "logit"), data = framework
)

m3_no_bonus_vcov <- sandwich::vcovHC(m3_no_bonus, type = "HC3")

lmtest::coeftest(m3_no_bonus, vcov. = m3_no_bonus_vcov)
lmtest::coefci(m3_no_bonus, vcov. = m3_no_bonus_vcov)


### 04. MODEL ADEQUACY ####

stopifnot(m3$converged, m3_no_bonus$converged)

# Compare observed and fitted means across the ordinal risk gradient
# to assess the assumed logit-linear relationship.

framework %>%
  mutate(
    fitted_value = fitted(m3) * 7,
    fitted_no_bonus = fitted(m3_no_bonus) * 6
  ) %>%
  group_by(risk_scale) %>%
  summarise(
    n = n(),
    observed = mean(ecological_value),
    predicted = mean(fitted_value),
    observed_no_bonus = mean(ecological_no_bonus),
    predicted_no_bonus = mean(fitted_no_bonus),
    .groups = "drop"
  )

# Assess sensitivity of regression inference to site-level clustering.
lmtest::coeftest(m3, vcov. = sandwich::vcovCL(
  m3, cluster = framework$site, type = "HC3"
))


### 05. SAVE FRAMEWORK OUTPUTS ####

saveRDS(
  list(
    framework = framework,
    decision_matrix = decision_matrix,
    risk_ecology = risk_ecology,
    risk_ecology_no_bonus = risk_ecology_no_bonus,
    m3 = m3,
    m3_vcov = m3_vcov,
    m3_no_bonus = m3_no_bonus,
    m3_no_bonus_vcov = m3_no_bonus_vcov
  ),
  file.path(fits_dir, "removal_framework.rds")
)


### FRAMEWORK CASE STUDIES ####

framework %>%
  select(
    dfg_id, site, dfg_type, depth_m, surface_area_m2,
    stability, degradation, risk_scale, ecological_value,
    coral_abundance, coral_recruits, genus_richness,
    growth_form_richness, fauna_abundance,
    cover_hard_coral, cover_macroalgae, cover_sediment,
    total_biotic_cover, action
  ) %>%
  arrange(action, desc(risk_scale), desc(ecological_value)) %>%
  print(n = Inf, width = Inf)

#### Illustrative management cases

example_ids <- paste("Net", c(86, 13, 45, 90, 20, 21))

framework %>%
  filter(dfg_id %in% example_ids) %>%
  select(
    dfg_id, site, dfg_type, depth_m, surface_area_m2,
    stability, degradation, risk_scale, ecological_value,
    coral_abundance, coral_recruits, genus_richness,
    growth_form_richness, fauna_abundance,
    starts_with("cover_"), total_biotic_cover, action
  ) %>%
  arrange(action, desc(risk_scale)) %>%
  print(n = Inf, width = Inf)
