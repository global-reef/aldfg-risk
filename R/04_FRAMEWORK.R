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
    score_growth = case_when(growth_form_richness >= 4 ~ 2L,
                             growth_form_richness >= 2 ~ 1L, TRUE ~ 0L),
    ecological_value = score_recruits + score_genus + score_growth,
    
    # Original stability bonus retained for sensitivity analysis only.
    score_stability = as.integer(risk_scale == 1),
    ecological_with_bonus = ecological_value + score_stability
  )

stopifnot(
  nrow(framework) == 106,
  !anyDuplicated(framework$dfg_id),
  !anyNA(framework %>% select(coral_recruits, genus_richness, growth_form_richness)),
  all(framework$risk_scale %in% 1:5),
  all(framework$ecological_value %in% 0:6),
  all(framework$ecological_with_bonus %in% 0:7)
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


### 02. REMOVAL PRIORITISATION FRAMEWORK ####

action_levels <- c("Prioritise removal", "Monitor", "Low removal priority")

classify_action <- function(risk, value) {
  factor(
    case_when(
      risk == 1 ~ "Low removal priority",
      value <= 2 ~ "Prioritise removal",
      risk %in% c(2, 3) & value >= 6 ~ "Low removal priority",
      TRUE ~ "Monitor"
    ),
    levels = action_levels
  )
}

framework <- framework %>%
  mutate(action = classify_action(risk_scale, ecological_value))

decision_matrix <- expand_grid(risk_scale = 1:5, ecological_value = 0:6) %>%
  mutate(action = classify_action(risk_scale, ecological_value)) %>%
  left_join(framework %>% count(risk_scale, ecological_value),
            by = c("risk_scale", "ecological_value")) %>%
  mutate(n = replace_na(n, 0L))

# Removing the stability bonus must not alter management classifications.
stopifnot(
  !anyNA(framework$action),
  sum(decision_matrix$n) == nrow(framework),
  all(framework$action ==
        classify_action(framework$risk_scale, framework$ecological_with_bonus))
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

# Fractional logistic regression for the bounded ecological index.
m3 <- glm(
  I(ecological_value / 6) ~ risk_scale,
  family = quasibinomial(link = "logit"), data = framework
)

m3_vcov <- sandwich::vcovHC(m3, type = "HC3")

lmtest::coeftest(m3, vcov. = m3_vcov)
lmtest::coefci(m3, vcov. = m3_vcov)


#### B. Stability bonus sensitivity ####

# Evaluate the influence of including stability directly in both indices.
risk_ecology_bonus <- cor.test(
  framework$risk_scale, framework$ecological_with_bonus,
  method = "spearman", exact = FALSE
)

risk_ecology_bonus

m3_bonus <- glm(
  I(ecological_with_bonus / 7) ~ risk_scale,
  family = quasibinomial(link = "logit"), data = framework
)

m3_bonus_vcov <- sandwich::vcovHC(m3_bonus, type = "HC3")

lmtest::coeftest(m3_bonus, vcov. = m3_bonus_vcov)
lmtest::coefci(m3_bonus, vcov. = m3_bonus_vcov)


### 04. MODEL ADEQUACY ####

stopifnot(m3$converged, m3_bonus$converged)

framework %>%
  mutate(
    fitted_value = fitted(m3) * 6,
    fitted_bonus = fitted(m3_bonus) * 7
  ) %>%
  group_by(risk_scale) %>%
  summarise(
    n = n(),
    observed = mean(ecological_value),
    predicted = mean(fitted_value),
    observed_bonus = mean(ecological_with_bonus),
    predicted_bonus = mean(fitted_bonus),
    .groups = "drop"
  )

# Assess sensitivity of inference to within-site dependence.
lmtest::coeftest(m3, vcov. = sandwich::vcovCL(
  m3, cluster = framework$site, type = "HC3"
))


### 05. SAVE FRAMEWORK OUTPUTS ####

saveRDS(
  list(
    framework = framework,
    decision_matrix = decision_matrix,
    risk_ecology = risk_ecology,
    risk_ecology_bonus = risk_ecology_bonus,
    m3 = m3,
    m3_vcov = m3_vcov,
    m3_bonus = m3_bonus,
    m3_bonus_vcov = m3_bonus_vcov
  ),
  file.path(fits_dir, "removal_framework.rds")
)

beep()