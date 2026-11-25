### 01. MODEL DATA PREPARATION ####

biotic_cols <- c("cover_hard_coral", "cover_soft_coral", "cover_cca",
                 "cover_sponge", "cover_tunicates", "cover_zoanthids",
                 "cover_corallimorphs", "cover_macroalgae", "cover_turf")

dfg <- dfg %>%
  mutate(
    log_area = log(surface_area_m2),
    log_effort = log(n_images),
    total_biotic_cover = rowSums(pick(all_of(biotic_cols))),
    benthic_richness = rowSums(pick(all_of(biotic_cols)) > 0)
  )

# Aggregate annotations by image, excluding unknown classifications from
# the denominator of the biotic cover model.
benthic_images <- benthos %>%
  group_by(dfg_id, image_id) %>%
  summarise(
    n_biotic = sum(!benthic_group %in% c("Sediment", "Unknown")),
    n_sediment = sum(benthic_group == "Sediment"),
    n_unknown = sum(benthic_group == "Unknown"),
    n_points = n(),
    .groups = "drop"
  ) %>%
  left_join(dfg %>% select(dfg_id, site, stability, degradation), by = "dfg_id")

# Standardise biological cover independently of total colonisation.
biotic_matrix <- dfg %>%
  select(all_of(biotic_cols)) %>%
  decostand(method = "total")

biotic_dist <- vegdist(sqrt(biotic_matrix), method = "bray")

stopifnot(
  all(dfg$coral_recruits <= dfg$coral_abundance),
  all(dfg$total_biotic_cover > 0),
  all(abs(rowSums(biotic_matrix) - 1) < 1e-8),
  all(benthic_images$n_points == 25),
  all(benthic_images$n_biotic + benthic_images$n_sediment > 0)
)


### 02. OBJECTIVE 1A: TOTAL BIOTIC COVER ####

#### A. Candidate models ####

m1a_bin <- glmmTMB(
  cbind(n_biotic, n_sediment) ~ stability + degradation + (1 | dfg_id),
  family = binomial, data = benthic_images
)

m1a_beta <- update(m1a_bin, family = betabinomial())
m1a_site <- update(m1a_beta, . ~ . + (1 | site))

AIC(m1a_bin, m1a_beta, m1a_site)

#### B. Primary model ####

# Beta-binomial accounts for extra-binomial variation within images.
# Site and gear random intercepts account for hierarchical sampling.

m1a <- m1a_site

summary(m1a)
drop1(m1a, test = "Chisq")

#### C. Sensitivity analysis ####

# Assess inference without accounting for site-level variation.
drop1(m1a_beta, test = "Chisq")


### 03. OBJECTIVE 1B: BIOLOGICAL COMMUNITY COMPOSITION ####

#### A. Multivariate dispersion ####

# Lingoes correction accommodates negative eigenvalues in Bray-Curtis distances.
disp_stability <- betadisper(biotic_dist, dfg$stability,
                             bias.adjust = TRUE, add = "lingoes")

disp_degradation <- betadisper(biotic_dist, dfg$degradation,
                               bias.adjust = TRUE, add = "lingoes")

permutest(disp_stability, permutations = 999)
permutest(disp_degradation, permutations = 999)

#### B. Primary model ####

# Condition on site to isolate within-site relationships between
# gear characteristics and relative biological composition.

m1b <- dbrda(
  biotic_dist ~ stability + degradation + cover_sediment + Condition(site),
  data = dfg
)

m1b_test <- anova(m1b, by = "margin", permutations = 999, strata = dfg$site)

m1b_test
RsquareAdj(m1b)

#### C. Sensitivity analyses ####

# Compare unrestricted and site-restricted PERMANOVA without
# explicitly removing between-site variation.

m1b_unrestricted <- adonis2(
  biotic_dist ~ stability + degradation + cover_sediment,
  data = dfg, permutations = 999, by = "margin"
)

m1b_restricted <- adonis2(
  biotic_dist ~ stability + degradation + cover_sediment,
  data = dfg, permutations = 999, by = "margin", strata = dfg$site
)

m1b_unrestricted
m1b_restricted

#### D. Ordination ####

# Unconstrained NMDS describes overall community patterns.
# Environmental and benthic category fits are descriptive and
# are not substitutes for the site-adjusted primary analysis.

nmds1b <- metaMDS(biotic_dist, k = 2, trymax = 200, trace = FALSE)

env1b <- envfit(
  nmds1b,
  dfg %>% select(stability, degradation, cover_sediment),
  permutations = 999
)

taxa1b <- envfit(nmds1b, biotic_matrix, permutations = 999)

nmds1b
env1b
taxa1b


### 04. OBJECTIVE 1C: BENTHIC CATEGORY RICHNESS ####

#### A. Candidate models ####

# Photographic effort is included because category detection
# depends on the number of images analysed.

m1c_pois <- glmmTMB(
  benthic_richness ~ stability + degradation + cover_sediment + log_effort,
  family = poisson, data = dfg
)

m1c_comp <- update(m1c_pois, family = compois())
m1c_site <- update(m1c_comp, . ~ . + (1 | site))

AIC(m1c_pois, m1c_comp, m1c_site)

#### B. Primary model ####

# COM-Poisson accounts for underdispersion in category richness.
m1c <- m1c_site

summary(m1c)
drop1(m1c, test = "Chisq")

#### C. Sensitivity analysis ####

drop1(m1c_comp, test = "Chisq")


### 05. OBJECTIVE 2A: CORAL AND ASSOCIATED FAUNA ABUNDANCE ####

#### A. Primary model ####

# Model coral and associated fauna abundance jointly without
# combining counts. Estimate the area-abundance relationship
# rather than assuming proportional scaling through an offset.

m2a <- manyglm(
  cbind(coral_abundance, fauna_abundance) ~
    stability + degradation + cover_sediment + log_area,
  family = "negative.binomial", data = dfg
)

summary(m2a)

# Nested model comparisons provide marginal predictor tests.
# Response-specific P-values are adjusted across both responses.

m2a_stability <- anova(
  update(m2a, . ~ . - stability), m2a,
  test = "LR", resamp = "pit.trap", p.uni = "adjusted", nBoot = 999
)

m2a_degradation <- anova(
  update(m2a, . ~ . - degradation), m2a,
  test = "LR", resamp = "pit.trap", p.uni = "adjusted", nBoot = 999
)

m2a_sediment <- anova(
  update(m2a, . ~ . - cover_sediment), m2a,
  test = "LR", resamp = "pit.trap", p.uni = "adjusted", nBoot = 999
)

m2a_area <- anova(
  update(m2a, . ~ . - log_area), m2a,
  test = "LR", resamp = "pit.trap", p.uni = "adjusted", nBoot = 999
)

m2a_stability
m2a_degradation
m2a_sediment
m2a_area

#### B. Site-level sensitivity analyses ####

# Response-specific mixed models assess robustness to site-level
# clustering because manyglm does not support random intercepts.

m2a_coral_site <- glmmTMB(
  coral_abundance ~ stability + degradation + cover_sediment +
    log_area + (1 | site),
  family = nbinom2(), data = dfg
)

m2a_fauna_site <- glmmTMB(
  fauna_abundance ~ stability + degradation + cover_sediment +
    log_area + (1 | site),
  family = nbinom2(), data = dfg
)

summary(m2a_coral_site)
summary(m2a_fauna_site)

drop1(m2a_coral_site, test = "Chisq")
drop1(m2a_fauna_site, test = "Chisq")


### 06. OBJECTIVE 2B: CORAL RECRUIT ABUNDANCE ####

#### A. Candidate models ####

m2b_pois <- glmmTMB(
  coral_recruits ~ stability + degradation + cover_sediment + log_area,
  family = poisson, data = dfg
)

m2b_nb <- update(m2b_pois, family = nbinom2())
m2b_site <- update(m2b_nb, . ~ . + (1 | site))

AIC(m2b_pois, m2b_nb, m2b_site)

#### B. Primary model ####

# Negative binomial accounts for overdispersion.
# Site adjustment does not substantially improve model fit.

m2b <- m2b_nb

summary(m2b)
drop1(m2b, test = "Chisq")

#### C. Site-level sensitivity analysis ####

summary(m2b_site)
drop1(m2b_site, test = "Chisq")


### 07. OBJECTIVE 2C: CORAL GENUS RICHNESS ####

#### A. Candidate models ####

m2c_pois <- glmmTMB(
  genus_richness ~ stability + degradation + cover_sediment + log_area,
  family = poisson, data = dfg
)

m2c_nb <- update(m2c_pois, family = nbinom2())
m2c_site <- update(m2c_nb, . ~ . + (1 | site))

AIC(m2c_pois, m2c_nb, m2c_site)

#### B. Primary model ####

# Negative binomial resolves overdispersion.
# Estimated site-level variance is negligible.

m2c <- m2c_nb

summary(m2c)
drop1(m2c, test = "Chisq")

#### C. Site-level sensitivity analysis ####

summary(m2c_site)
drop1(m2c_site, test = "Chisq")

#### D. Growth-form richness sensitivity analysis ####

# Assess consistency between taxonomic and morphological richness.
m2c_gf_pois <- update(m2c, growth_form_richness ~ ., family = poisson)
m2c_gf_nb <- update(m2c_gf_pois, family = nbinom2())
m2c_gf_comp <- update(m2c_gf_pois, family = compois())

AIC(m2c_gf_pois, m2c_gf_nb, m2c_gf_comp)

m2c_gf <- m2c_gf_comp

summary(m2c_gf)
drop1(m2c_gf, test = "Chisq")


### 08. MODEL DIAGNOSTICS ####

#### A. Convergence and collinearity ####

primary_glms <- list(m1a = m1a, m1c = m1c, m2b = m2b, m2c = m2c)

sensitivity_glms <- list(
  m1a_beta = m1a_beta, m1c_comp = m1c_comp,
  m2a_coral_site = m2a_coral_site, m2a_fauna_site = m2a_fauna_site,
  m2b_site = m2b_site, m2c_site = m2c_site, m2c_gf = m2c_gf
)

map_dfr(c(primary_glms, sensitivity_glms), function(x) {
  tibble(convergence = x$fit$convergence, positive_hessian = x$sdr$pdHess)
}, .id = "model")

map(primary_glms, check_collinearity)

#### B. Simulated residual diagnostics ####

candidate_glms <- list(
  m1a_bin = m1a_bin, m1c_pois = m1c_pois,
  m2b_pois = m2b_pois, m2c_pois = m2c_pois,
  m2c_gf_pois = m2c_gf_pois, m2c_gf_nb = m2c_gf_nb
)

residuals_primary <- map(primary_glms, simulateResiduals, n = 1000)
residuals_sensitivity <- map(sensitivity_glms, simulateResiduals, n = 1000)
residuals_candidates <- map(candidate_glms, simulateResiduals, n = 1000)

# Numerical diagnostics without plotting.
for (model_name in names(c(primary_glms, sensitivity_glms, candidate_glms))) {
  cat("\n", model_name, "\n")
  res <- c(residuals_primary, residuals_sensitivity, residuals_candidates)[[model_name]]
  print(DHARMa::testDispersion(res, plot = FALSE))
  print(DHARMa::testUniformity(res, plot = FALSE))
}

DHARMa::testZeroInflation(residuals_primary$m2b, plot = FALSE)
DHARMa::testZeroInflation(residuals_primary$m2c, plot = FALSE)
DHARMa::testZeroInflation(residuals_sensitivity$m2a_coral_site, plot = FALSE)
DHARMa::testZeroInflation(residuals_sensitivity$m2a_fauna_site, plot = FALSE)

DHARMa::testOutliers(residuals_primary$m1a, type = "bootstrap", plot = FALSE)

rm(res, model_name)

#### C. Benthic richness prediction bounds ####

# COM-Poisson is unbounded; assess predictions against the maximum
# number of available biological categories.

range(predict(m1c, type = "response"))

mean(residuals_primary$m1c$simulatedResponse > length(biotic_cols))
mean(residuals_primary$m1c$simulatedResponse == 0)


### 09. SAVE MODEL OUTPUTS ####

saveRDS(
  list(
    m1a = m1a, m1b = m1b, m1b_test = m1b_test, m1c = m1c,
    m2a = m2a, m2b = m2b, m2c = m2c,
    m2a_stability = m2a_stability, m2a_degradation = m2a_degradation,
    m2a_sediment = m2a_sediment, m2a_area = m2a_area,
    nmds1b = nmds1b, env1b = env1b, taxa1b = taxa1b,
    disp_stability = disp_stability, disp_degradation = disp_degradation
  ),
  file.path(fits_dir, "primary_models.rds")
)

saveRDS(
  c(sensitivity_glms, list(
    m1b_unrestricted = m1b_unrestricted,
    m1b_restricted = m1b_restricted
  )),
  file.path(fits_dir, "sensitivity_models.rds")
)

# Save diagnostic plots directly to PDF without displaying them in RStudio.
diagnostics_dir <- file.path(plots_dir, "diagnostics")
dir.create(diagnostics_dir, recursive = TRUE, showWarnings = FALSE)

pdf(file.path(diagnostics_dir, "primary_model_diagnostics.pdf"),
    width = 10, height = 6)

for (model_name in names(residuals_primary)) {
  plot(residuals_primary[[model_name]], title = model_name)
}

plot(m2a)
stressplot(nmds1b)

plotResiduals(residuals_primary$m1a, benthic_images$stability)
plotResiduals(residuals_primary$m1a, benthic_images$degradation)
plotResiduals(residuals_primary$m1c, dfg$cover_sediment)
plotResiduals(residuals_primary$m1c, dfg$log_effort)
plotResiduals(residuals_primary$m2b, dfg$log_area)
plotResiduals(residuals_primary$m2b, dfg$cover_sediment)
plotResiduals(residuals_primary$m2c, dfg$log_area)
plotResiduals(residuals_primary$m2c, dfg$cover_sediment)

dev.off()

pdf(file.path(diagnostics_dir, "sensitivity_model_diagnostics.pdf"),
    width = 10, height = 6)

for (model_name in names(residuals_sensitivity)) {
  plot(residuals_sensitivity[[model_name]], title = model_name)
}

dev.off()

writeLines(capture.output(sessionInfo()), file.path(stats_dir, "session_info.txt"))

# Remove temporary diagnostic collections after saving.
rm(primary_glms, sensitivity_glms, candidate_glms,
   residuals_primary, residuals_sensitivity, residuals_candidates, model_name)

beep()