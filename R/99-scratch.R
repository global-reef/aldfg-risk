### 01. EXPANDED ECOLOGICAL INDEX: SETUP ####

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

classify_action_expanded <- function(risk, value) {
  factor(case_when(
    risk == 1 ~ "Low removal priority",
    value <= 3 ~ "Prioritise removal",
    risk %in% c(2, 3) & value >= 6 ~ "Low removal priority",
    TRUE ~ "Monitor"
  ), levels = action_levels)
}


### 02. BUILD WORKING EXPANDED INDEX ####

eco_candidate <- dfg %>%
  mutate(
    risk_scale = case_when(
      stability == "Stable" ~ 1L,
      stability == "Semi-stable" & degradation %in% c("Low", "Moderate") ~ 2L,
      stability == "Semi-stable" & degradation %in% c("High", "Severe") ~ 3L,
      stability == "Mobile" & degradation %in% c("Low", "Moderate") ~ 4L,
      stability == "Mobile" & degradation %in% c("High", "Severe") ~ 5L
    ),
    
    ## current framework
    score_recruits_current = case_when(
      coral_recruits >= 6 ~ 2L,
      coral_recruits >= 1 ~ 1L,
      TRUE ~ 0L
    ),
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
      score_recruits_current + score_genus_current + score_growth_current,
    action_current = classify_action_current(risk_scale, ecological_value_current),
    
    ## expanded candidate metrics
    hard_coral_present = as.integer(cover_hard_coral > 0),
    
    # Replace the single hard-coral benthic category with genus richness,
    # while retaining at least one coral contribution when hard coral cover
    # is present but no genus was recorded in situ.
    combined_benthic_richness =
      benthic_richness - hard_coral_present +
      pmax(genus_richness, hard_coral_present),
    
    score_recruits = case_when(
      coral_recruits >= 6 ~ 2L,
      coral_recruits >= 1 ~ 1L,
      TRUE ~ 0L
    ),
    
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
    
    ecological_value_expanded =
      score_recruits + score_hc_cover + score_combined_richness + score_fauna,
    
    action_expanded = classify_action_expanded(risk_scale, ecological_value_expanded)
  )

stopifnot(
  all(eco_candidate$ecological_value_current >= 0),
  all(eco_candidate$ecological_value_current <= 6),
  all(eco_candidate$ecological_value_expanded >= 0),
  all(eco_candidate$ecological_value_expanded <= 8)
)


### 03. QUICK CHECKS ####

eco_candidate %>%
  summarise(
    current_mean = mean(ecological_value_current),
    expanded_mean = mean(ecological_value_expanded),
    current_sd = sd(ecological_value_current),
    expanded_sd = sd(ecological_value_expanded)
  )

eco_candidate %>%
  select(
    coral_recruits, cover_hard_coral, combined_benthic_richness,
    fauna_abundance, ecological_value_current, ecological_value_expanded
  ) %>%
  cor(method = "spearman") %>%
  round(2)

eco_candidate %>%
  select(score_recruits, score_hc_cover, score_combined_richness, score_fauna) %>%
  cor(method = "spearman") %>%
  round(2)

eco_candidate %>%
  select(
    Current = action_current,
    Expanded = action_expanded
  ) %>%
  pivot_longer(everything(), names_to = "framework", values_to = "priority") %>%
  count(framework, priority, .drop = FALSE) %>%
  group_by(framework) %>%
  mutate(percent = round(100 * n / sum(n), 1)) %>%
  print(n = Inf)

eco_candidate %>%
  filter(action_current != action_expanded) %>%
  select(
    dfg_id, risk_scale,
    coral_recruits, cover_hard_coral, combined_benthic_richness, fauna_abundance,
    ecological_value_current, ecological_value_expanded,
    action_current, action_expanded
  ) %>%
  arrange(desc(risk_scale), dfg_id) %>%
  print(n = Inf, width = Inf)


### 04. TRADE-OFF CHECK ####

cor.test(
  eco_candidate$risk_scale,
  eco_candidate$ecological_value_expanded,
  method = "spearman",
  exact = FALSE
)

m3_expanded <- glm(
  I(ecological_value_expanded / 8) ~ risk_scale,
  family = quasibinomial(link = "logit"),
  data = eco_candidate
)

lmtest::coeftest(m3_expanded, vcov. = sandwich::vcovHC(m3_expanded, type = "HC3"))

lmtest::coeftest(
  m3_expanded,
  vcov. = sandwich::vcovCL(m3_expanded, cluster = eco_candidate$site, type = "HC3")
)


### 05. DECISION MATRIX FOR EXPANDED INDEX ####

decision_matrix_expanded <- expand_grid(
  risk_scale = 1:5,
  ecological_value_expanded = 0:8
) %>%
  mutate(
    action = classify_action_expanded(risk_scale, ecological_value_expanded)
  )

matrix_plot_data_expanded <- decision_matrix_expanded %>%
  mutate(
    risk_group = case_when(
      risk_scale == 1 ~ "Low (1)",
      risk_scale <= 3 ~ "Medium (2-3)",
      TRUE ~ "High (4-5)"
    ),    
    risk_group = factor(
      risk_group,
      levels = c("Low (1)", "Medium (2-3)", "High (4-5)")
    ),
    value_group = case_when(
      ecological_value_expanded <= 3 ~ "Low (0–3)",
      ecological_value_expanded <= 5 ~ "Medium (4–5)",
      TRUE ~ "High (6–8)"
    ),
    value_group = factor(
      value_group,
      levels = c("High (6–8)", "Medium (4–5)", "Low (0–3)")
    ),
    action = factor(action, levels = action_levels)
  ) %>%
  distinct(risk_group, value_group, action)

stopifnot(nrow(matrix_plot_data_expanded) == 9)


### 06. CONCEPTUAL FIGURE: EXPANDED FRAMEWORK ####

theme_matrix <- theme_figure +
  theme(
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_text(colour = "grey20"),
    axis.title = element_text(face = "bold"),
    legend.position = "none",
    plot.title = element_blank(),
    plot.subtitle = element_blank()
  )

score_cols <- c(
  "0" = "#F5F3EE",
  "1" = "#B9C7B5",
  "2" = "#477F91"
)


#### A. Entanglement risk ####

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

p4a_expanded <- ggplot(risk_matrix, aes(degradation, stability)) +
  geom_tile(aes(fill = risk_scale), colour = "white", linewidth = 0.8) +
  geom_text(aes(label = risk_scale, colour = risk_scale >= 4), size = 4.5) +
  scale_fill_gradientn(colours = heatmap_cols, limits = c(1, 5)) +
  scale_colour_manual(values = c("FALSE" = "grey20", "TRUE" = "white")) +
  scale_x_discrete(position = "top", expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(
    x = "Degradation",
    y = "Physical\nstability"
  ) +
  theme_matrix


#### B. Expanded ecological value ####

ecological_matrix_expanded <- tribble(
  ~metric,                     ~`0`,  ~`1`,    ~`2`,
  "Recruit abundance",         "0",   "1–5",   "≥6",
  "Hard coral cover (%)",      "<10", "10–30", ">30",
  "Combined benthic richness", "0–2", "3–5",   "≥6",
  "Associated fauna",          "0",   "1–5",   "≥6"
) %>%
  pivot_longer(-metric, names_to = "score", values_to = "threshold") %>%
  mutate(
    score = factor(score, levels = c("0", "1", "2")),
    metric = factor(metric, levels = rev(c(
      "Recruit abundance",
      "Hard coral cover (%)",
      "Combined benthic richness",
      "Associated fauna"
    )))
  )

p4b_expanded <- ggplot(ecological_matrix_expanded, aes(score, metric)) +
  geom_tile(aes(fill = score), colour = "white", linewidth = 0.8) +
  geom_text(aes(label = threshold, colour = score == "2"), size = 4.2) +
  scale_fill_manual(values = score_cols, drop = FALSE) +
  scale_colour_manual(values = c("FALSE" = "grey20", "TRUE" = "white")) +
  scale_x_discrete(
    position = "top",
    expand = c(0, 0),
    labels = c("0" = "0\npoints", "1" = "1\npoint", "2" = "2\npoints")
  ) +
  scale_y_discrete(
    labels = function(x) stringr::str_wrap(x, width = 16),
    expand = c(0, 0)
  ) +
  labs(
    x = "Component\nscore",
    y = "Ecological\nmetric"
  ) +
  theme_matrix


#### C. Prioritisation matrix ####

p4c_expanded <- ggplot(matrix_plot_data_expanded, aes(risk_group, value_group)) +
  geom_tile(aes(fill = action), colour = "white", linewidth = 0.8) +
  geom_text(
    aes(
      label = recode(
        as.character(action),
        "Prioritise removal" = "Priority",
        "Monitor" = "Monitor",
        "Low removal priority" = "Low"
      ),
      colour = action == "Low removal priority"
    ),
    size = 4.5
  ) +
  scale_fill_manual(values = action_cols, drop = FALSE) +
  scale_colour_manual(values = c("FALSE" = "grey20", "TRUE" = "white")) +
  scale_x_discrete(position = "top", expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(
    x = "Entanglement\nrisk",
    y = "Ecological\nvalue"
  ) +
  theme_matrix


#### D. Assemble ####

fig4_alt_expanded <- (p4a_expanded | p4b_expanded) /
  p4c_expanded +
  plot_layout(heights = c(1, 1.15)) +
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")")

print(fig4_alt_expanded)

save_figure(fig4_alt_expanded, "Fig4_expanded_framework_draft", "half", dir = figures_dir)
#### Management classification summary ####

eco_candidate %>%
  select(Current = action_current, Expanded = action_expanded) %>%
  pivot_longer(everything(), names_to = "framework", values_to = "priority") %>%
  count(framework, priority, .drop = FALSE) %>%
  group_by(framework) %>%
  mutate(percent = round(100 * n / sum(n), 1)) %>%
  ungroup() %>%
  print(n = Inf)

table(Current = eco_candidate$action_current,
      Expanded = eco_candidate$action_expanded)