### 01. CALCULATE RISK AND ECOLOGICAL VALUE ####

action_levels <- c("Prioritise removal", "Monitor", "Low removal priority")

classify_action_current <- function(risk, value) {
  factor(
    case_when(
      risk == 1 ~ "Low removal priority",
      value <= 2 ~ "Prioritise removal",
      risk %in% c(2L, 3L) & value >= 6 ~ "Low removal priority",
      TRUE ~ "Monitor"
    ),
    levels = action_levels
  )
}

classify_action <- function(risk, value) {
  factor(case_when(
    risk == 1 ~ "Low removal priority",
    value <= 3 ~ "Prioritise removal",
    risk %in% c(2, 3) & value >= 6 ~ "Low removal priority",
    TRUE ~ "Monitor"
  ), levels = action_levels)
}


risk_matrix <- expand_grid(
  stability = factor(
    c("Mobile", "Semi-stable", "Stable"),
    levels = c("Mobile", "Semi-stable", "Stable")
  ),
  degradation = factor(
    c("Low", "Moderate", "High", "Severe"),
    levels = c("Low", "Moderate", "High", "Severe")
  )
) %>%
  mutate(
    risk_scale = case_when(
      stability == "Stable" ~ 1L,
      stability == "Semi-stable" & degradation %in% c("Low", "Moderate") ~ 2L,
      stability == "Semi-stable" & degradation %in% c("High", "Severe") ~ 3L,
      stability == "Mobile" & degradation %in% c("Low", "Moderate") ~ 4L,
      stability == "Mobile" & degradation %in% c("High", "Severe") ~ 5L
    )
  )

# Validate inputs before scoring; missing values must not be treated as zero.
stopifnot(!anyNA(dfg %>% select(stability, degradation, coral_recruits,
                                genus_richness, growth_form_richness, cover_hard_coral, benthic_richness, fauna_abundance)))

framework <- dfg %>%
  left_join(risk_matrix, by = c("stability", "degradation")) %>%
  mutate(
    
    # Earlier coral-only index retained for sensitivity analysis.
    score_recruits = case_when(coral_recruits >= 6 ~ 2L,
                               coral_recruits >= 1 ~ 1L, TRUE ~ 0L),
    score_genus_current = case_when(
      genus_richness >= 3 ~ 2L,
      genus_richness == 2 ~ 1L,
      TRUE ~ 0L
    ),
    score_growth_current = case_when(
      growth_form_richness >= 4 ~ 2L,
      growth_form_richness >= 2 ~ 1L,
      TRUE ~ 0L
    ),
    ecological_value_current =
      score_recruits + score_genus_current + score_growth_current,
    action_current = classify_action_current(risk_scale, ecological_value_current),
    
    # Expanded primary index
    hard_coral_present = as.integer(cover_hard_coral > 0),
    
    # Replace the single hard-coral benthic category with genus richness,
    # while retaining at least one coral contribution when hard coral cover
    # is present but no genus was recorded in situ.
    combined_benthic_richness =
      benthic_richness - hard_coral_present +
      pmax(genus_richness, hard_coral_present),
    
    score_hc_cover = case_when(
      cover_hard_coral > 30 ~ 2L,
      cover_hard_coral >= 10 ~ 1L, TRUE ~ 0L
    ),
    
    score_fauna = case_when(
      fauna_abundance >= 6 ~ 2L,
      fauna_abundance >= 1 ~ 1L, TRUE ~ 0L
    ),
    
    score_combined_richness = case_when(
      combined_benthic_richness >= 6 ~ 2L,
      combined_benthic_richness >= 3 ~ 1L,
      TRUE ~ 0L
    ),
    
    ecological_value =
      score_recruits + score_hc_cover + score_combined_richness + score_fauna,
    
    action = classify_action(risk_scale, ecological_value),
    ecological_with_bonus = ecological_value_current + as.integer(risk_scale == 1)
  )

stopifnot(nrow(framework) == 106, !anyDuplicated(framework$dfg_id),
          all(framework$risk_scale %in% 1:5),
          all(framework$ecological_value %in% 0:8),
          all(framework$ecological_value_current %in% 0:6),
          all(framework$ecological_with_bonus %in% 0:7), !anyNA(framework$action),
          all(framework$action_current ==
                classify_action_current(framework$risk_scale, framework$ecological_with_bonus)))

### 02. REMOVAL PRIORITISATION FRAMEWORK ####

decision_matrix <- expand_grid(risk_scale = 1:5, ecological_value = 0:8) %>%
  mutate(action = classify_action(risk_scale, ecological_value)) %>%
  left_join(framework %>% count(risk_scale, ecological_value),
            by = c("risk_scale", "ecological_value")) %>%
  mutate(n = replace_na(n, 0L))

stopifnot(sum(decision_matrix$n) == nrow(framework))
framework %>% count(action, .drop = FALSE) %>% mutate(percent = 100 * n / sum(n))
table(framework$risk_scale, framework$action)

### 03. RISK AND ECOLOGICAL VALUE TRADE-OFF ####

#### A. Primary analysis ####

risk_ecology <- cor.test(framework$risk_scale, framework$ecological_value,
                         method = "spearman", exact = FALSE)
# Fractional logistic regression for the bounded ecological index.
m3 <- glm(I(ecological_value / 8) ~ risk_scale,
          family = quasibinomial(link = "logit"), data = framework)
m3_vcov <- sandwich::vcovHC(m3, type = "HC3")
m3_site_vcov <- sandwich::vcovCL(m3, cluster = framework$site, type = "HC3")
risk_ecology
lmtest::coeftest(m3, vcov. = m3_vcov)
lmtest::coefci(m3, vcov. = m3_vcov)
lmtest::coeftest(m3, vcov. = m3_site_vcov)

#### B. Earlier coral-only index sensitivity ####

risk_ecology_current <- cor.test(framework$risk_scale, framework$ecological_value_current,
                                 method = "spearman", exact = FALSE)
m3_current <- glm(I(ecological_value_current / 6) ~ risk_scale,
                  family = quasibinomial(link = "logit"), data = framework)
m3_current_vcov <- sandwich::vcovHC(m3_current, type = "HC3")
risk_ecology_current
lmtest::coeftest(m3_current, vcov. = m3_current_vcov)

# Historical stability bonus comparison, not part of the expanded index.
risk_ecology_bonus <- cor.test(framework$risk_scale, framework$ecological_with_bonus,
                               method = "spearman", exact = FALSE)
m3_bonus <- glm(I(ecological_with_bonus / 7) ~ risk_scale,
                family = quasibinomial(link = "logit"), data = framework)
m3_bonus_vcov <- sandwich::vcovHC(m3_bonus, type = "HC3")
risk_ecology_bonus
lmtest::coeftest(m3_bonus, vcov. = m3_bonus_vcov)

#### C. Index comparison ####

framework %>%
  select(coral_recruits, cover_hard_coral, combined_benthic_richness, fauna_abundance,
         ecological_value_current, ecological_value) %>% cor(method = "spearman") %>% round(2)
framework %>%
  select(score_recruits, score_hc_cover, score_combined_richness, score_fauna) %>%
  cor(method = "spearman") %>% round(2)
table(Current = framework$action_current, Expanded = framework$action)

### 04. MODEL ADEQUACY ####

stopifnot(m3$converged, m3_current$converged, m3_bonus$converged)
framework %>% mutate(predicted = fitted(m3) * 8) %>%
  group_by(risk_scale) %>%
  summarise(n = n(), observed = mean(ecological_value), predicted = mean(predicted), .groups = "drop")

### 05. SAVE FRAMEWORK OUTPUTS ####

saveRDS(list(framework = framework, decision_matrix = decision_matrix, risk_matrix = risk_matrix,
             risk_ecology = risk_ecology, m3 = m3, m3_vcov = m3_vcov, m3_site_vcov = m3_site_vcov,
             risk_ecology_current = risk_ecology_current, m3_current = m3_current,
             m3_current_vcov = m3_current_vcov, risk_ecology_bonus = risk_ecology_bonus,
             m3_bonus = m3_bonus, m3_bonus_vcov = m3_bonus_vcov),
        file.path(fits_dir, "removal_framework.rds"))

beep()
