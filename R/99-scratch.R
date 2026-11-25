### 08. SELECT PRIMARY MODELS ####

m1a <- m1a_beta
m1c <- m1c_comp
m2b <- m2b_nb
m2c <- m2c_nb

# Check convergence and positive-definite Hessians.
sapply(list(m1a = m1a, m1c = m1c, m2b = m2b, m2c = m2c),
       function(x) x$sdr$pdHess)

sapply(list(m1a = m1a, m1c = m1c, m2b = m2b, m2c = m2c),
       function(x) x$fit$convergence)


### 09. FINAL MODEL DIAGNOSTICS ####

# Investigate residual non-uniformity in the beta-binomial model.
plotResiduals(res1a_beta, benthic_images$stability)
plotResiduals(res1a_beta, benthic_images$degradation)
testOutliers(res1a_beta, type = "bootstrap")

# Assess residual patterns against continuous predictors.
plotResiduals(res1c_comp, dfg$cover_sediment)
plotResiduals(res1c_comp, dfg$log_effort)

plotResiduals(res2b_nb, dfg$log_area)
plotResiduals(res2b_nb, dfg$cover_sediment)

plotResiduals(res2c_nb, dfg$log_area)
plotResiduals(res2c_nb, dfg$cover_sediment)

# Assess whether simulated richness exceeds the nine available categories
# or generates biologically impossible negative counts.
range(predict(m1c, type = "response"))

mean(res1c_comp$simulatedResponse > length(biotic_cols))
mean(res1c_comp$simulatedResponse == 0)

# Check whether Euclidean correction affects dispersion tests.
permutest(
  betadisper(biotic_dist, dfg$stability, bias.adjust = TRUE, add = "lingoes"),
  permutations = 999
)

permutest(
  betadisper(biotic_dist, dfg$degradation, bias.adjust = TRUE, add = "lingoes"),
  permutations = 999
)

# Assess collinearity among model predictors.
check_collinearity(m1a)
check_collinearity(m1c)
check_collinearity(m2b)
check_collinearity(m2c)

### 10. OVERALL PREDICTOR TESTS ####

drop1(m1a, test = "Chisq")
drop1(m1c, test = "Chisq")
drop1(m2b, test = "Chisq")
drop1(m2c, test = "Chisq")

### 11. SITE-LEVEL SENSITIVITY ANALYSES ####

# Assess whether accounting for spatial clustering changes inference.

m1a_site <- update(m1a, . ~ . + (1 | site))
m1c_site <- update(m1c, . ~ . + (1 | site))
m2b_site <- update(m2b, . ~ . + (1 | site))
m2c_site <- update(m2c, . ~ . + (1 | site))

AIC(m1a, m1a_site)
AIC(m1c, m1c_site)
AIC(m2b, m2b_site)
AIC(m2c, m2c_site)

summary(m1a_site)
summary(m1c_site)
summary(m2b_site)
summary(m2c_site)

drop1(m1a_site, test = "Chisq")
drop1(m1c_site, test = "Chisq")
drop1(m2b_site, test = "Chisq")
drop1(m2c_site, test = "Chisq")

# Restrict permutations within sites for community composition.
m1b_site <- adonis2(
  biotic_dist ~ stability + degradation + cover_sediment,
  data = dfg, permutations = 999, by = "margin", strata = dfg$site
)

m1b_site

# Response-specific mixed models assess site-level sensitivity because
# manyglm does not support random intercepts.

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


### 12. GROWTH-FORM RICHNESS SENSITIVITY ####

# Assess consistency between taxonomic and morphological richness responses.

m2c_gf_pois <- update(m2c, growth_form_richness ~ ., family = poisson)
m2c_gf_nb <- update(m2c_gf_pois, family = nbinom2())

AIC(m2c_gf_pois, m2c_gf_nb)

summary(m2c_gf_pois)
summary(m2c_gf_nb)

res2c_gf_pois <- simulateResiduals(m2c_gf_pois, n = 1000)
res2c_gf_nb <- simulateResiduals(m2c_gf_nb, n = 1000)

plot(res2c_gf_pois)
plot(res2c_gf_nb)

testDispersion(res2c_gf_pois)
testDispersion(res2c_gf_nb)

### 13. FINAL MODEL REFINEMENTS ####

#### Site-adjusted benthic models ####

res1a_site <- simulateResiduals(m1a_site, n = 1000)
res1c_site <- simulateResiduals(m1c_site, n = 1000)

testDispersion(res1a_site)
testUniformity(res1a_site)
testDispersion(res1c_site)
testUniformity(res1c_site)

#### Site-adjusted community composition ####

# Separate within-site associations from between-site variation.
m1b_partial <- dbrda(
  biotic_dist ~ stability + degradation + cover_sediment + Condition(site),
  data = dfg
)

anova(m1b_partial, by = "margin", permutations = 999, strata = dfg$site)

#### Growth-form richness ####

m2c_gf_comp <- update(m2c_gf_pois, family = compois())

AIC(m2c_gf_pois, m2c_gf_nb, m2c_gf_comp)

summary(m2c_gf_comp)
drop1(m2c_gf_comp, test = "Chisq")

res2c_gf_comp <- simulateResiduals(m2c_gf_comp, n = 1000)

testDispersion(res2c_gf_comp)
testUniformity(res2c_gf_comp)


### 01. COMMUNITY COMPOSITION AUDIT ####

# Reproduce the original cover-based approach using corrected data.
cover_dist <- vegdist(sqrt(dfg %>% select(all_of(biotic_cols))),
                      method = "bray")

adonis2(cover_dist ~ stability + degradation + cover_sediment,
        data = dfg, permutations = 999, by = "margin")

# Assess the overall degradation association without adjusting
# for other gear characteristics, but controlling for site.
adonis2(biotic_dist ~ degradation + Condition(site),
        data = dfg, permutations = 999, strata = dfg$site)