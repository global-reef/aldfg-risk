### 01. FIGURE SETUP ####
#### Publication dimensions ####

figure_width <- 7.2
half_page_height <- 4.5
full_page_height <- 8
figure_dpi <- 600

#### Figure export helper ####
save_figure <- function(plot,
                        filename,
                        size = c("half", "full"),
                        width = figure_width,
                        dir = plots_dir,
                        dpi = figure_dpi,
                        bg = "white") {
  
  size <- match.arg(size)
  height <- if (size == "half") half_page_height else full_page_height
  
  png_path <- file.path(dir, paste0(filename, ".png"))
  pdf_path <- file.path(dir, paste0(filename, ".pdf"))
  
  ggsave(filename = png_path,plot = plot,
         width = width,height = height,
         units = "in", dpi = dpi, bg = bg)
  ggsave(filename = pdf_path, plot = plot, width = width,
         height = height, units = "in", bg = bg, device = "pdf")
}



# Match plotting labels to the cleaned ReefCloud cover variables.
cover_labels <- c(
  cover_hard_coral = "Hard coral",
  cover_soft_coral = "Soft coral",
  cover_cca = "CCA",
  cover_sponge = "Sponge",
  cover_tunicates = "Tunicates",
  cover_zoanthids = "Zoanthids",
  cover_corallimorphs = "Corallimorphs",
  cover_macroalgae = "Macroalgae",
  cover_turf = "Turf algae"
)

theme_figure <- theme_clean +
  theme(
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 9, colour = "grey20"),
    strip.text = element_text(size = 10),
    legend.text = element_text(size = 9),
    legend.title = element_text(size = 9),
    plot.tag = element_text(size = 11, face = "bold"),
    axis.line = element_line(colour = "grey30", linewidth = 0.3),
    axis.ticks = element_line(colour = "grey30", linewidth = 0.3)
  )

theme_set(theme_figure)

theme_composite <- theme(plot.tag = element_text(face = "bold", size = 11))


### 02. FIGURE 1: STUDY AREA ####

# Study area and survey locations prepared separately in QGIS.


### 03. FIGURE 2: BENTHIC INTEGRATION ####

#### A. Biological community composition ####

nmds_data <- as.data.frame(vegan::scores(nmds1b, display = "sites")) %>%
  bind_cols(dfg %>% select(dfg_id, stability, degradation))

nmds_vectors <- as.data.frame(vegan::scores(taxa1b, display = "vectors")) %>%
  rownames_to_column("category") %>%
  mutate(category = recode(category, !!!cover_labels)) 

p2a <- ggplot(nmds_data, aes(NMDS1, NMDS2, colour = stability)) +
  stat_ellipse(aes(fill = stability), geom = "polygon",
               level = 0.95, alpha = 0.12, colour = NA) +
  stat_ellipse(level = 0.95, linewidth = 0.5) +
  geom_point(aes(shape = stability), size = 1.5, alpha = 0.95) +
  scale_colour_manual(values = stability_cols) +
  scale_fill_manual(values = stability_cols) +
  scale_shape_manual(values = stability_shapes) +
  annotate("text", x = Inf, y = Inf,
           label = paste0("Stress = ", round(nmds1b$stress, 3)),
           hjust = 1.1, vjust = 1.5, size = 3) +
  coord_equal() +
  labs(x = "NMDS1", y = "NMDS2",
       colour = NULL, fill = NULL, shape = NULL)


#### B. Total biotic cover ####

# Match observed cover to the identified-point denominator used in O1A.
cover_emm <- emmeans(m1a, ~ stability, type = "response") %>%
  summary(infer = c(TRUE, FALSE)) %>%
  as.data.frame()

p2b <- ggplot(dfg, aes(stability)) +
  geom_jitter(aes(y = total_biotic_cover * n_points / n_identified,
                  colour = stability),
              width = 0.14, height = 0, alpha = 0.6, size = 1.4) +
  geom_pointrange(
    data = cover_emm,
    aes(y = prob * 100, ymin = asymp.LCL * 100, ymax = asymp.UCL * 100),
    colour = "black", linewidth = 0.5, size = 0.7
  ) +
  scale_colour_manual(values = stability_cols, guide = "none") +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 25)) +
  labs(x = "Structural stability", y = "Biotic cover (%)")


#### C. Benthic category richness ####

# Predictions at mean photographic effort and the most common physical
# condition: semi-stable gear with high degradation.
richness_predictions <- tibble(
  cover_sediment = seq(min(dfg$cover_sediment),
                       max(dfg$cover_sediment), length.out = 100),
  stability = factor("Semi-stable", levels = levels(dfg$stability)),
  degradation = factor("High", levels = levels(dfg$degradation)),
  log_effort = mean(dfg$log_effort),
  site = dfg$site[1]
)

richness_pred <- predict(m1c, newdata = richness_predictions,
                         type = "link", re.form = NA, se.fit = TRUE)

richness_predictions <- richness_predictions %>%
  mutate(
    predicted = exp(richness_pred$fit),
    lower_CL = exp(richness_pred$fit - qnorm(0.975) * richness_pred$se.fit),
    upper_CL = exp(richness_pred$fit + qnorm(0.975) * richness_pred$se.fit)
  )

p2c <- ggplot(dfg, aes(cover_sediment, benthic_richness)) +
  geom_ribbon(
    data = richness_predictions,
    aes(x = cover_sediment, ymin = lower_CL, ymax = upper_CL),
    inherit.aes = FALSE, fill = "#477F91", alpha = 0.2
  ) +
  geom_point(colour = "#477F91", alpha = 0.4, size = 1.7,
             position = position_jitter(height = 0.08, width = 0)) +
  geom_line(data = richness_predictions,
            aes(y = predicted), colour = "#1F4E6B", linewidth = 0.9) +
  scale_y_continuous(breaks = seq(0, 9, 3), limits = c(0, 9)) +
  labs(x = "Sediment cover (%)", y = "Benthic richness")

rm(richness_pred)


#### D. Assemble Figure 2 ####

fig2 <- (p2a | (p2b / p2c)) +
  plot_layout(widths = c(1.25, 1)) +
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") &
  theme_composite

print(fig2)


### 04. FIGURE 3: CORAL AND ASSOCIATED FAUNA ABUNDANCE ####

#### A. Coral abundance and recruitment ####

# Total coral abundance includes recruits.
coral_summary <- dfg %>%
  group_by(stability) %>%
  summarise(
    Recruits = mean(coral_recruits),
    `Other corals` = mean(coral_abundance - coral_recruits),
    .groups = "drop"
  ) %>%
  pivot_longer(-stability, names_to = "coral_group",
               values_to = "mean_abundance") %>%
  mutate(coral_group = factor(coral_group,
                              levels = c("Other corals", "Recruits")))

p3a <- ggplot(coral_summary,
              aes(stability, mean_abundance, fill = coral_group)) +
  geom_col(width = 0.65, colour = "white", linewidth = 0.3) +
  scale_fill_manual(values = coral_cols,
                    breaks = c("Recruits", "Other corals")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  labs(x = "Structural stability",
       y = "Mean coral abundance per gear item", fill = NULL)


#### B. Associated fauna by stability ####

p3b <- ggplot(dfg, aes(stability, fauna_abundance)) +
  geom_jitter(aes(colour = stability), width = 0.16, height = 0,
              size = 1.8, alpha = 0.5) +
  stat_summary(fun = mean, geom = "point",
               shape = 23, size = 3, fill = "black", colour = "white") +
  scale_colour_manual(values = stability_cols, guide = "none") +
  scale_y_continuous(limits = c(0, NA)) +
  labs(x = "Structural stability", y = "Associated fauna")


#### C. Associated fauna by degradation ####

p3c <- ggplot(dfg, aes(degradation, fauna_abundance)) +
  geom_jitter(aes(colour = degradation), width = 0.16, height = 0,
              size = 1.8, alpha = 0.5) +
  stat_summary(fun = mean, geom = "point",
               shape = 23, size = 3, fill = "black", colour = "white") +
  scale_colour_manual(values = degradation_cols, guide = "none") +
  scale_y_continuous(limits = c(0, NA)) +
  labs(x = "Degree of degradation", y = "Associated fauna")


#### D. Assemble Figure 3 ####

fig3 <- (p3a | (p3b / p3c)) +
  plot_layout(widths = c(1.15, 1)) +
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") &
  theme_composite

print(fig3)


### 05a. FIGURE 4: REMOVAL PRIORITISATION FRAMEWORK ####

#### A. Decision matrix ####

matrix_plot_data <- decision_matrix %>%
  mutate(
    risk_group = case_when(
      risk_scale == 1 ~ "Low (1)",
      risk_scale <= 3 ~ "Medium (2–3)",
      TRUE ~ "High (4–5)"
    ),
    value_group = case_when(
      ecological_value <= 2 ~ "Low (0–2)",
      ecological_value <= 5 ~ "Medium (3–5)",
      TRUE ~ "High (6)"
    ),
    risk_group = factor(risk_group, levels = c("Low (1)", "Medium (2–3)", "High (4–5)")),
    value_group = factor(value_group, levels = c("Low (0–2)", "Medium (3–5)", "High (6)"))
  ) %>%
  distinct(risk_group, value_group, action)

p4a <- ggplot(matrix_plot_data, aes(risk_group, value_group)) +
  geom_tile(aes(fill = action), colour = "white", linewidth = 0.8) +
  scale_fill_manual(values = action_cols) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(x = "Entanglement risk", y = "Ecological value", fill = NULL) +
  guides(fill = guide_legend(ncol = 1)) +
  theme(
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "right",
    legend.direction = "vertical",
    legend.justification = "center"
  )

#### B. Risk-value trade-off ####

p4b <- ggplot(framework, aes(risk_scale, ecological_value)) +
  geom_ribbon(
    data = risk_predictions,
    aes(x = risk_scale, ymin = lower_CL, ymax = upper_CL),
    inherit.aes = FALSE, fill = "#477F91", alpha = 0.22
  ) +
  geom_jitter(width = 0.12, height = 0.08,
              colour = "#477F91", size = 1.7, alpha = 0.35) +
  geom_line(data = risk_predictions, aes(y = predicted),
            colour = "#1F4E6B", linewidth = 0.9) +
  scale_x_continuous(breaks = 1:5) +
  scale_y_continuous(breaks = 0:6) +
  coord_cartesian(ylim = c(-0.2, 6.2)) +
  labs(x = "Entanglement risk score", y = "Ecological value score")


#### C. Management priorities ####

p4c <- action_summary %>%
  mutate(action = factor(action, levels = rev(action_levels))) %>%
  ggplot(aes(n, action, fill = action)) +
  geom_col(width = 0.65) +
  geom_text(aes(label = sprintf("%d (%.1f%%)", n, percent)),
            hjust = -0.12, size = 3) +
  scale_fill_manual(values = action_cols, guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.55))) +
  labs(x = "Number of gear items", y = NULL) +
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank())


#### D. Assemble Figure 4 ####

fig4 <- (
  p4b | ((p4a + theme(legend.position = "none")) / p4c)
) +
  plot_layout(widths = c(1, 1.2)) +
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") &
  theme_composite

print(fig4)


save_figure(fig4, "Fig4_removal_framework", "half", dir = figures_dir)


### 05b. ALTERNATIVE FIGURE 4: CONCEPTUAL FRAMEWORK ####

#### A. Shared matrix theme ####

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


#### B. Entanglement risk assessment ####

risk_matrix <- expand_grid(
  stability = factor(c("Mobile", "Semi-stable", "Stable"),
                     levels = c("Mobile", "Semi-stable", "Stable")),
  degradation = factor(c("Low", "Moderate", "High", "Severe"),
                       levels = c("Low", "Moderate", "High", "Severe"))
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

p4a_alt <- ggplot(risk_matrix, aes(degradation, stability)) +
  geom_tile(aes(fill = risk_scale), colour = "white", linewidth = 0.8) +
  geom_text(aes(label = risk_scale, colour = risk_scale >= 4), size = 4.5) +
  scale_fill_gradientn(colours = heatmap_cols, limits = c(1, 5)) +
  scale_colour_manual(values = c("FALSE" = "grey20", "TRUE" = "white")) +
  scale_x_discrete(position = "top", expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(x = "Degradation", y = "Physical stability") +
  theme_matrix


#### C. Ecological value assessment ####

ecological_matrix <- tribble(
  ~metric,                ~`0`,  ~`1`,  ~`2`,
  "Recruit abundance",    "0",   "1–5", "≥6",
  "Genus richness",       "0–1", "2",   "≥3",
  "Growth form richness", "0–1", "2–3", "≥4"
) %>%
  pivot_longer(-metric, names_to = "score", values_to = "threshold") %>%
  mutate(
    score = as.integer(score),
    metric = factor(metric, levels = rev(c(
      "Recruit abundance", "Genus richness", "Growth form richness"
    )))
  )

p4b_alt <- ggplot(ecological_matrix, aes(factor(score), metric)) +
  geom_tile(aes(fill = score), colour = "white", linewidth = 0.8) +
  geom_text(aes(label = threshold, colour = score == 2), size = 4.5) +
  scale_fill_gradientn(colours = c("#F5F3EE", "#B9C7B5", "#477F91"),
                       limits = c(0, 2)) +
  scale_colour_manual(values = c("FALSE" = "grey20", "TRUE" = "white")) +
  scale_x_discrete(position = "top", expand = c(0, 0),
                   labels = c("0 points", "1 point", "2 points")) +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 12),
                   expand = c(0, 0)) +
  labs(x = "Component score", y = "Ecological metric") +
  theme_matrix


#### D. Removal prioritisation matrix ####

matrix_plot_data <- decision_matrix %>%
  mutate(
    risk_group = case_when(
      risk_scale == 1 ~ "Low (1)",
      risk_scale <= 3 ~ "Medium (2–3)",
      TRUE ~ "High (4–5)"
    ),
    value_group = case_when(
      ecological_value <= 2 ~ "Low (0–2)",
      ecological_value <= 5 ~ "Medium (3–5)",
      TRUE ~ "High (6)"
    ),
    risk_group = factor(risk_group, levels = c("Low (1)", "Medium (2–3)", "High (4–5)")),
    value_group = factor(value_group, levels = c("High (6)", "Medium (3–5)", "Low (0–2)"))
  ) %>%
  distinct(risk_group, value_group, action)

stopifnot(nrow(matrix_plot_data) == 9)

p4c_alt <- ggplot(matrix_plot_data, aes(risk_group, value_group)) +
  geom_tile(aes(fill = action), colour = "white", linewidth = 0.8) +
  geom_text(
    aes(
      label = recode(as.character(action),
                     "Prioritise removal" = "Priority",
                     "Monitor" = "Monitor",
                     "Low removal priority" = "Low"),
      colour = action == "Low removal priority"
    ),
    size = 4.5
  ) +
  scale_fill_manual(values = action_cols) +
  scale_colour_manual(values = c("FALSE" = "grey20", "TRUE" = "white")) +
  scale_x_discrete(position = "top", expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(x = "Entanglement risk", y = "Ecological value") +
  theme_matrix


#### E. Assemble conceptual framework ####

fig4_alt <- wrap_plots(
  A = p4a_alt,
  B = p4b_alt,
  C = p4c_alt,
  design = "AB
            CC",
  heights = c(1, 1)
) +
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")")

print(fig4_alt)


### 06. FIGURE 5: ILLUSTRATIVE MANAGEMENT CASES ####

# Photographic placeholders pending selection of field images.
# Examples illustrate contrasting pathways rather than representative samples.

case_data <- framework_examples %>%
  group_by(action) %>%
  mutate(case_number = row_number()) %>%
  ungroup() %>%
  mutate(
    physical_label = paste0(
      site, " | ", str_to_sentence(dfg_type), "\n",
      stability, " | ", degradation, " degradation"
    ),
    score_label = paste0(
      "Risk: ", risk_scale, "/5    |    Value: ", ecological_value, "/7"
    ),
    biological_label = paste0(
      "Corals: ", coral_abundance, " (", coral_recruits, " recruits)\n",
      "Genera: ", genus_richness, "    |    Growth forms: ", growth_form_richness,
      "\nAssociated fauna: ", fauna_abundance
    )
  )

fig5 <- ggplot(case_data) +
  geom_rect(aes(xmin = 0.04, xmax = 0.96, ymin = 0.45, ymax = 0.89),
            fill = "grey94", colour = "grey80", linewidth = 0.3) +
  geom_text(aes(x = 0.5, y = 0.95, label = dfg_id),
            fontface = "bold", size = 3.8) +
  annotate("text", x = 0.5, y = 0.67,
           label = "Field photograph", colour = "grey55", size = 3.2) +
  geom_text(aes(x = 0.5, y = 0.37, label = physical_label),
            size = 3, lineheight = 1) +
  geom_label(aes(x = 0.5, y = 0.25, label = score_label, fill = action),
             colour = "grey10", size = 3, linewidth = 0) +
  geom_text(aes(x = 0.5, y = 0.10, label = biological_label),
            size = 2.9, lineheight = 1.1) +
  facet_grid(case_number ~ action) +
  scale_fill_manual(values = action_cols_traffic, guide = "none") +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  theme_void(base_family = "sans") +
  theme(
    strip.text.x = element_text(face = "bold", size = 11,
                                margin = margin(b = 6)),
    strip.text.y = element_blank(),
    panel.spacing = unit(0.5, "lines")
  )

print(fig5)


### 07. FIGURE S1: INDIVIDUAL BENTHIC COMPOSITION ####

#### A. Prepare composition data ####

# Retain absolute percentage cover to represent both colonisation
# extent and biological composition.

composition_plot_data <- dfg %>%
  arrange(stability, total_biotic_cover) %>%
  mutate(dfg_id = factor(dfg_id, levels = unique(dfg_id))) %>%
  select(dfg_id, stability, all_of(biotic_cols)) %>%
  pivot_longer(all_of(biotic_cols), names_to = "category", values_to = "cover") %>%
  mutate(
    category = recode(category, !!!cover_labels),
    category = factor(category, levels = names(benthic_cols))
  )


#### B. Plot stability groups ####

plot_composition <- function(groups) {
  composition_plot_data %>%
    filter(stability %in% groups) %>%
    ggplot(aes(dfg_id, cover, fill = category)) +
    geom_col(width = 0.9) +
    facet_grid(. ~ stability, scales = "free_x", space = "free_x") +
    scale_fill_manual(values = benthic_cols, drop = FALSE,
                      breaks = names(benthic_cols)[1:9]) +
    scale_x_discrete(labels = function(x) str_remove(x, "^Net "),
                     expand = expansion(add = 0.5)) +
    scale_y_continuous(
      breaks = seq(0, 100, 25),
      labels = scales::label_number(suffix = "%"),
      expand = expansion(mult = c(0, 0.02))
    ) +
    coord_cartesian(ylim = c(0, 100)) +
    guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
    labs(x = "Fishing gear ID", y = "Biological cover", fill = NULL) +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 7),
      axis.ticks.x = element_blank(),
      panel.spacing.x = unit(1, "lines"),
      panel.grid.major.y = element_line(colour = "grey90", linewidth = 0.3),
      strip.text = element_text(face = "bold")
    )
}


#### C. Assemble Figure S1 ####

figS1 <- (
  (plot_composition(c("Mobile", "Stable")) + labs(x = NULL)) /
    plot_composition("Semi-stable")
) +
  plot_layout(guides = "collect", heights = c(1, 1)) &
  theme(legend.position = "bottom")

print(figS1)


### 08. FIGURE S2: BENTHIC CATEGORY COVER HEATMAP ####

#### A. Mean observed percentage cover ####

heatmap_data <- bind_rows(
  dfg %>%
    group_by(stability) %>%
    summarise(across(c(total_biotic_cover, all_of(biotic_cols),
                       cover_sediment), mean), .groups = "drop") %>%
    rename(condition = stability) %>%
    mutate(predictor = "Structural stability"),
  
  dfg %>%
    group_by(degradation) %>%
    summarise(across(c(total_biotic_cover, all_of(biotic_cols),
                       cover_sediment), mean), .groups = "drop") %>%
    rename(condition = degradation) %>%
    mutate(predictor = "Degree of degradation")
) %>%
  pivot_longer(-c(condition, predictor), names_to = "category",
               values_to = "mean_cover") %>%
  mutate(
    category = recode(category, !!!cover_labels,
                      total_biotic_cover = "Total biotic cover",
                      cover_sediment = "Sediment"),
    category = factor(category, levels = rev(c(
      "Total biotic cover", names(benthic_cols)[1:9], "Sediment"
    ))),
    condition = factor(condition, levels = c(
      "Mobile", "Semi-stable", "Stable", "Low", "Moderate", "High", "Severe"
    )),
    predictor = factor(predictor, levels = c(
      "Structural stability", "Degree of degradation"
    )),
    cover_label = case_when(
      mean_cover == 0 ~ "0",
      mean_cover < 0.1 ~ "<0.1",
      TRUE ~ sprintf("%.1f", mean_cover)
    )
  )


#### B. Assemble heatmap ####

figS2 <- ggplot(heatmap_data, aes(condition, category, fill = mean_cover)) +
  geom_tile(colour = "white", linewidth = 0.7) +
  geom_text(aes(label = cover_label,
                colour = mean_cover >= 60), size = 3.3) +
  facet_grid(. ~ predictor, scales = "free_x", space = "free_x") +
  scale_fill_gradientn(
    colours = heatmap_cols, limits = c(0, 100),
    breaks = seq(0, 100, 25), name = "Mean cover (%)"
  ) +
  scale_colour_manual(values = c("FALSE" = "grey15", "TRUE" = "white"),
                      guide = "none") +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(x = NULL, y = NULL) +
  theme(
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text.x = element_text(size = 9),
    axis.text.y = element_text(size = 9),
    strip.text = element_text(face = "bold"),
    panel.spacing.x = unit(0.7, "lines"),
    legend.key.width = unit(1.2, "cm")
  )

print(figS2)



### 09. SAVE ALL FIGURES ####

figures_dir <- file.path(plots_dir, "manuscript")
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

save_figure(fig2, "Fig2_benthic_integration", "half", dir = figures_dir)
save_figure(fig3, "Fig3_ecological_value", "half", dir = figures_dir)
save_figure(fig4_alt, "Fig4_removal_framework", "half", dir = figures_dir)
# save_figure(fig5, "Fig5_management_cases_draft", "full", dir = figures_dir)
save_figure(figS1, "FigS1_benthic_composition", "full", dir = figures_dir)
save_figure(figS2, "FigS2_benthic_heatmap", "half", dir = figures_dir)

beep()