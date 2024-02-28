### Statistical Analysis ###

packages <- c('dplyr', 'tidyr', 'ggplot2', 'ggpubr', 'ggeffects',
              'Hmisc', 'lme4', 'MuMIn', 'visreg', "nlme")
lapply(packages, require, character.only = T)

setwd("~/Desktop/Research/FireBiodiversity/Data/ProcessedData")
responses_herbs_all <- read.csv("herbaceous_responses_all.csv")
responses_grass_all <- read.csv("grass_responses_all.csv")
responses_nongrass_all <- read.csv("nongrass_responses_all.csv")
responses_woody_all <- read.csv("woody_responses_all.csv")

# functions for figure generation and formatting


pct_chng <- function(x) {
  return (paste0(round( (100*(exp(x) -1)), 2), "%"))
} 

supp_figure <- function(df_RR, lower_lim) {
  plot(ggplot(df_RR, aes(x = FireFreq, y = logRR2, color = Vegetation)) +
         geom_point(size = I(3)) +
         ylab('% difference in species richness from unburned') +
         xlab('fire frequency (fires per year)') +
         theme_classic() + 
         theme(axis.ticks.length = unit(.3, 'cm'), 
               axis.title = element_text(size = 16), axis.text = element_text(size = 12),
               plot.title = element_text(hjust = 0.5),
               legend.position = "none") +
         scale_y_continuous(labels = function(i) pct_chng(i), expand = c(0, 0), 
                            limits = c(-log(lower_lim), log(5)), breaks = c(-log(lower_lim),log(1), log(2), log(3), log(4), log(5))) +
         scale_x_continuous(expand = c(0, 0), limits = c(0,1.01), breaks = seq(0, 1, by = .2)) +
         geom_hline(yintercept=0, linetype='dashed') +
         geom_errorbar(aes(ymin = logRR2 - logRR2_v, ymax = logRR2 + logRR2_v), width = 0.02) +
         coord_cartesian(clip= "off") +
         scale_color_manual(values=c("green4", "gold")))
}

supp_figure(responses_herbs_all, 5)
supp_figure(responses_grass_all, 5)
supp_figure(responses_nongrass_all, 5)
supp_figure(responses_woody_all, 40)


figure <- function(gg) {
  plot(gg +
         ylab('% difference in species richness from unburned') +
         xlab('fire frequency (fires per year)') +
         theme_classic() + 
         theme(axis.ticks.length = unit(.3, 'cm'), 
               axis.title = element_text(size = 16), axis.text = element_text(size = 12, color = "black"), 
               axis.ticks = element_line(color = "black"),
               plot.title = element_text(hjust = 0.5),
               legend.position = "none") +
         scale_y_continuous(labels = function(i) pct_chng(i), expand = c(0, 0), 
                            limits = c(-log(5), log(4)), breaks = c(-log(5),log(1), log(2), log(3), log(4))) +
         scale_x_continuous(expand = c(0, 0), limits = c(0,1.015), breaks = seq(0, 1, by = .2)) +
         geom_hline(yintercept=0, linetype='dashed') +
         scale_color_manual(values=c("green4", "gold"))) +
    coord_cartesian(clip = "off")
}

# library('mgcv') 
# library('gratia')

# # GAM
# rcorr(as.matrix(responses_herbs_all[,c(11,9:10,12)]),type="pearson")
# responses_herbs_all$Vegetation <- as.factor(responses_herbs_all$Vegetation)
# responses_herbs_all$Site <- as.factor(responses_herbs_all$Site)
# herbs_gam_all <- gam(logRR2 ~ s(scale(FireFreq), by = Vegetation) + scale(Aridity) + scale(MAT) + scale(Duration) + s(Site, bs = "re"), method = "REML", data = responses_herbs_all, select = T)
# summary(herbs_gam_all)
# appraise(herbs_gam_all)
# draw(herbs_gam_all, residuals = TRUE, scales = "fixed")
# 
# rcorr(as.matrix(responses_grass_all[,c(11,9:10,12)]),type="pearson")
# responses_grass_all$Vegetation <- as.factor(responses_grass_all$Vegetation)
# responses_grass_all$Site <- as.factor(responses_grass_all$Site)
# grass_gam_all <- gam(logRR2 ~ s(scale(FireFreq), by = Vegetation) + scale(Aridity) + scale(MAT) + scale(Duration) + s(Site, bs = "re"), method = "REML", data = responses_grass_all, select = T)
# summary(grass_gam_all)
# appraise(grass_gam_all)
# draw(grass_gam_all, residuals = TRUE, scales = "fixed")
# 
# 
# rcorr(as.matrix(responses_nongrass_all[,c(11,9:10,12)]),type="pearson")
# responses_nongrass_all$Vegetation <- as.factor(responses_nongrass_all$Vegetation)
# responses_nongrass_all$Site <- as.factor(responses_nongrass_all$Site)
# nongrass_gam_all <- gam(logRR2 ~ s(scale(FireFreq), by = Vegetation) + scale(Aridity) + scale(MAT) + scale(Duration) + s(Site, bs = "re"), method = "REML", data = responses_nongrass_all, select = T)
# summary(nongrass_gam_all)
# appraise(nongrass_gam_all)
# draw(nongrass_gam_all, residuals = TRUE, scales = "fixed")
# 
# rcorr(as.matrix(responses_woody_all[,c(11,9:10,12)]),type="pearson")
# responses_woody_all$Vegetation <- as.factor(responses_woody_all$Vegetation)
# responses_woody_all$Site <- as.factor(responses_woody_all$Site)
# woody_gam_all <- gam(logRR2 ~ s(scale(FireFreq), by = Vegetation) + scale(Aridity) + scale(MAT) + scale(Duration), method = "REML", data = responses_woody_all, select = T)
# summary(woody_gam_all)
# appraise(woody_gam_all)
# draw(woody_gam_all, residuals = TRUE, scales = "fixed")


### Herbaceous species richness for all sites ###

# Correlation matrix
rcorr(as.matrix(responses_herbs_all[,c(11,9:10,12)]),type="pearson")

# full model and model selection
herbs_mod_all <- lmer(logRR2 ~ (scale(FireFreq) + scale(FireFreq^2) + scale(Aridity) + scale(MAT) + scale(Duration))*Vegetation + (1|Site), data = responses_herbs_all, REML = F, na.action = "na.fail")
options(max.print=1000000)
dredge(herbs_mod_all)

# best model and summary
best_herbs_all <- lme(logRR2 ~ scale(FireFreq) + scale(FireFreq^2) + scale(MAT)*Vegetation, random = ~1|Site, data = responses_herbs_all)
best_herbs_all_pctchng <- lme(pct_chng ~ scale(FireFreq) + scale(FireFreq^2) + scale(MAT)*Vegetation, random = ~1|Site, data = responses_herbs_all)

# Model diagnostics
plot(best_herbs_all)
qqnorm(resid(best_herbs_all))
qqline(resid(best_herbs_all))

summary(best_herbs_all)
r.squaredGLMM(best_herbs_all)

# model with unscaled predictors
responses_herbs_all$Vegetation <- as.factor()
best_herbs_all <- lme(logRR2 ~ FireFreq + I(FireFreq^2) + MAT*Vegetation, random = ~1|Site, data = responses_herbs_all)
summary(best_herbs_all)

# Examples for in-text results
examples <- data.frame(FireFreq= c(0.25, 0.5, 1), 
                       MAT = rep(6.267937, 3),
                       Vegetation = rep("Savanna", 3),
                       Site = rep("Cedar Creek", 3),
                       preds = rep(NA, 3),
                       pct_chg = rep(NA, 3))

examples$preds <- predict(best_herbs_all, examples)
examples$pct_chg <- 100*(exp(predict(best_herbs_all, examples))-1)

# visualization
gg_herbaceous <- visreg(lmer(logRR2 ~ scale(FireFreq) + scale(FireFreq^2) + scale(MAT) * Vegetation + (1|Site), data = responses_herbs_all, na.action = "na.fail"),
                        "FireFreq", gg = T, line = list(col = "black")) + 
  ggtitle("Herbaceous vegetation") +
  geom_point(aes(color=responses_herbs_all$Vegetation, size = I(3)))
figure(gg_herbaceous)


## Grass species richness for all sites ##

# correlation matrix
rcorr(as.matrix(responses_grass_all[,c(11,9:10,12)]),type="pearson")

# full model and model selection
grass_mod_all <- lmer(logRR2 ~ (scale(FireFreq) + scale(FireFreq^2) + scale(Aridity) + scale(MAT) + scale(Duration))*Vegetation + (1|Site), data = responses_grass_all, REML = F, na.action = "na.fail")
dredge(grass_mod_all)

# best model and summary
best_grass_all <- lme(logRR2 ~ scale(FireFreq) + scale(FireFreq^2) * Vegetation + scale(MAT), random = ~1|Site, data = responses_grass_all)
summary(best_grass_all)
r.squaredGLMM(best_grass_all)

# model with unscaled predictors
best_grass_all <- lme(logRR2 ~ FireFreq + I(FireFreq^2) * Vegetation + MAT, random = ~1|Site, data = responses_grass_all)
summary(best_grass_all)

examples <- data.frame(FireFreq= c(0.25, 0.5, 1, 0.25, 0.5, 1), 
                       MAT = rep(6.267937, 6),
                       Vegetation = as.factor(c("Savanna", "Savanna", "Savanna", "Forest", "Forest", "Forest")),
                       Site = rep("Cedar Creek", 6),
                       preds = rep(NA, 6),
                       pct_chg = rep(NA, 6))

examples$preds <- predict(best_herbs_all, examples)
examples$pct_chg <- 100*(exp(predict(best_herbs_all, examples))-1)

# visualization
gg_grass <- visreg(lmer(logRR2 ~ scale(FireFreq) + scale(FireFreq^2) * Vegetation + scale(MAT) + (1|Site), data = responses_grass_all, na.action = "na.fail"),
                   "FireFreq", by = "Vegetation", line = "Vegetation", gg = T, overlay = T) +  
  ggtitle("Grasses") +
  geom_point(aes(color=Vegetation, size = I(3)))
figure(gg_grass)



## All sites nongrass species richness

# correlation matrix
rcorr(as.matrix(responses_nongrass_all[,c(11,9:10,12)]),type="pearson")

# full model and model selection
nongrass_mod_all <- lmer(logRR2 ~ (scale(FireFreq) + scale(FireFreq^2) + scale(MAT) + scale(Aridity) + scale(Duration))*Vegetation + (1|Site), data = responses_nongrass_all, REML = F, na.action = "na.fail")
dredge(nongrass_mod_all)

# best model and summary
best_nongrass_all <- lme(logRR2 ~ scale(FireFreq) + scale(FireFreq^2) + (scale(MAT) + scale(Aridity) + scale(Duration))*Vegetation, random = ~1|Site, data = responses_nongrass_all)
summary(best_nongrass_all)
r.squaredGLMM(best_nongrass_all)

# model with unscaled predictors
best_nongrass_all <- lme(logRR2 ~ FireFreq + I(FireFreq^2) + (MAT + Aridity + Duration)*Vegetation, random = ~1|Site, data = responses_nongrass_all)
summary(best_nongrass_all)


gg_nongrass <- visreg(lmer(logRR2 ~ scale(FireFreq) + scale(FireFreq^2) + (scale(MAT) + scale(Aridity) + scale(Duration))*Vegetation + (1|Site), data = responses_nongrass_all, na.action = "na.fail"),
                      "FireFreq", overlay = T, gg = T, line = list(col = "black")) +  
  ggtitle("Non-grass herbaceous vegetation") +
  geom_point(aes(color=responses_nongrass_all$Vegetation, size = I(3)))
figure(gg_nongrass)


## All sites woody species richness ##

# correlation matrix
rcorr(as.matrix(responses_woody_all[,c(11,9:10,12)]),type="pearson")

# full model and model selection

woody_mod_all <- lmer(logRR2 ~ (scale(FireFreq) + scale(FireFreq^2) + scale(Aridity) + scale(MAT) + scale(Duration))*Vegetation + (1|Site), data = responses_woody_all, REML = F, na.action = "na.fail")
dredge(woody_mod_all)

# best model and summary
best_woody_all <- lme(logRR2 ~ scale(FireFreq) + (scale(FireFreq^2)), random = ~1|Site, data = responses_woody_all)
summary(best_woody_all)
r.squaredGLMM(best_woody_all)

# model with unscaled predictors
best_woody_all <- lme(logRR2 ~ FireFreq + I(FireFreq^2), random = ~1|Site, data = responses_woody_all)
summary(best_woody_all)

examples <- data.frame(FireFreq= c(0.25, 0.5, 1, 0.25, 0.5, 1), 
                       Seasonality = rep(median(responses_woody_all$Seasonality), 6),
                       Vegetation = as.factor(c("Savanna", "Savanna", "Savanna", "Forest", "Forest", "Forest")),
                       Site = rep("Cedar Creek"), 6)


gg_woody <- visreg(lmer(logRR2 ~ scale(FireFreq^2) + scale(FireFreq) + (1|Site), data = responses_woody_all, na.action = "na.fail"),
                   "FireFreq", overlay = T, gg = T, line = list(col = "black")) +  
  ggtitle("Woody vegetation") +
  geom_point(aes(color=responses_woody_all$Vegetation, size = I(3)))
figure(gg_woody)

visreg(lmer(logRR2 ~ scale(FireFreq^2) + (scale(FireFreq) + scale(Seasonality)) * Vegetation + (1|Site), data = responses_woody_all, na.action = "na.fail"),
       "Seasonality", by = "Vegetation", gg = T)



## Plot-level herbaceous models ##

figure2 <- function(gg) {
  plot(gg +
         ylab('log response ratio') +
         xlab('fire frequency (fires per year)') +
         theme_classic() + 
         theme(axis.ticks.length = unit(.3, 'cm'), 
               axis.title = element_text(size = 16), axis.text = element_text(size = 12, color = "black"), 
               axis.ticks = element_line(color = "black"),
               plot.title = element_text(hjust = 0.5),
               legend.position = "none") +
         scale_y_continuous(expand = c(0, 0), limits = c(-1,1), breaks = seq(-1, 1, by = .5)) +
         scale_x_continuous(expand = c(0, 0), limits = c(0,1.015), breaks = seq(0, 1, by = .2)) +
         geom_hline(yintercept=0, linetype='dashed') +
         scale_color_manual(values=c("green4", "gold"))) +
    coord_cartesian(clip = "off")
}





responses_herbaceous2 <- read.csv("herbaceous_responses_plotlevel.csv")

# correlation matrix 
rcorr(as.matrix(responses_herbaceous2[,c(29,27:28,30)]),type="pearson")

# Species richness # 
herbs_mod_richness <- lmer(logRR2 ~ (scale(FireFreq) + scale(FireFreq^2)) * Vegetation + scale(Aridity) + scale(MAT) + scale(Duration) + (1|Site), data = responses_herbaceous2, REML = F, na.action = "na.fail")
dredge(herbs_mod_richness)

best_mod_richness <- lme(logRR2 ~ scale(FireFreq) + scale(FireFreq^2), random = ~1|Site, data = responses_herbaceous2)
summary(best_mod_richness)
r.squaredGLMM(best_mod_richness)

best_mod_richness <- lme(logRR2 ~ FireFreq + I(FireFreq^2), random = ~1|Site, data = responses_herbaceous2)
summary(best_mod_richness)

gg_richness <- visreg(lmer(logRR2 ~ scale(FireFreq) + scale(FireFreq^2) + (1|Site), data = responses_herbaceous2, REML = F, na.action = "na.fail"), 
                      "FireFreq", gg = T, line = list(col = "black")) +
  ggtitle("Species richness") +
  geom_point(aes(color=responses_herbaceous2$Vegetation), size = I(3)) 
figure2(gg_richness)


# Shannon diversity # 

herbs_mod_shannon <- lmer(logRR_H2 ~ (scale(FireFreq) + scale(FireFreq^2)) * Vegetation + scale(Aridity) + scale(MAT) + scale(Duration) + (1|Site), data = responses_herbaceous2, REML = F, na.action = "na.fail")
dredge(herbs_mod_shannon)

best_mod_shannon <- lme(logRR_H2 ~ scale(FireFreq) + scale(FireFreq^2), random = ~1|Site, data = responses_herbaceous2)
summary(best_mod_shannon)
r.squaredGLMM(best_mod_shannon)

best_mod_shannon <- lme(logRR_H2 ~ FireFreq + I(FireFreq^2), random = ~1|Site, data = responses_herbaceous2)
summary(best_mod_shannon)

gg_shannon <- visreg(lmer(logRR_H2 ~ scale(FireFreq) + scale(FireFreq^2) + (1|Site), data = responses_herbaceous2, REML = F, na.action = "na.fail"),
                     "FireFreq", gg = T, line = list(col = "black")) +  
  ggtitle ("Shannon diversity") +
  geom_point(aes(color=responses_herbaceous2$Vegetation), size = I(3))
figure2(gg_shannon)

# Evenness

herbs_mod_evenness <- lmer(logRR_J2 ~ (scale(FireFreq) + scale(FireFreq^2)) * Vegetation + scale(Aridity) + scale(MAT) + scale(Duration) + (1|Site), data = responses_herbaceous2, REML = F, na.action = "na.fail")
dredge(herbs_mod_evenness)


best_mod_evenness <- lme(logRR_J2 ~ scale(FireFreq) + scale(FireFreq^2), random = ~1|Site, data = responses_herbaceous2)
summary(best_mod_evenness)
r.squaredGLMM(best_mod_evenness)

best_mod_evenness <- lme(logRR_J2 ~ FireFreq + I(FireFreq^2), random = ~1|Site, data = responses_herbaceous2)
summary(best_mod_evenness)

gg_evenness <- visreg(lmer(logRR_J2 ~ scale(FireFreq) + scale(FireFreq^2) + scale(MAP) + (1|Site), data = responses_herbaceous2, REML = F, na.action = "na.fail"), "FireFreq", gg = T, line = list(col = "black")) +  
  ggtitle("Evenness") +
  geom_point(aes(color=responses_herbaceous2$Vegetation), size = I(3))
figure2(gg_evenness)


# Community dissimilarity herbaceous
herbaceous_dissimilarity <- read.csv("herbaceous_dissimilarity.csv")

ggplot(herbaceous_dissimilarity, aes(x = FireFreq, y = BC_dissimilarity, color = Vegetation)) +
       geom_point() + theme_classic() +
       scale_y_continuous(expand = c(0, 0), limits = c(0.01, 1.01)) +
       scale_x_continuous(expand = c(0, 0), limits = c(-0.01, 1.01), breaks = seq(0,1,.2)) +
       scale_color_manual(values=c("green4", "gold"))

herb_dissimilarity_mod <- lmer(BC_dissimilarity ~ (scale(FireFreq) + scale(FireFreq^2) + scale(Aridity) + scale(MAT) + scale(Duration))*Vegetation + (1|Site), data = herbaceous_dissimilarity, REML = F, na.action = "na.fail")
dredge(herb_dissimilarity_mod)

# rcorr(as.matrix(herbaceous_dissimilarity[,c(6:9)]),type="pearson")
# 
# library(glmmTMB)
# herbaceous_dissimilarity$Site <- as.factor(herbaceous_dissimilarity$Site)
# BC_model_herbaceous <- glmmTMB(BC_dissimilarity ~  scale(FireFreq)*Vegetation + scale(FireFreq^2)*Vegetation + scale(MAP) + scale(MAT) + scale(Duration) + (1|Site), herbaceous_dissimilarity, family=beta_family(link="logit"), na.action = "na.fail")
# head(dredge(BC_model_herbaceous))
# 
# BC_herbaceous_best <- glmmTMB(BC_dissimilarity ~  scale(FireFreq) + (1|Site), herbaceous_dissimilarity, family=beta_family(link="logit"), na.action = "na.fail")
# summary(BC_herbaceous_best)
# 
# plot(ggpredict(BC_herbaceous_best, terms = "FireFreq"), residuals = T) # this goes beyond the domain, so have to fix this somehow
# 
# BC_herbaceous_best <- glmmTMB(BC_dissimilarity ~  FireFreq + (1|Site), herbaceous_dissimilarity, family=beta_family(link="logit"), na.action = "na.fail")
# summary(BC_herbaceous_best)
# 
# # Dissimilarity grass

grass_dissimilarity <- read.csv("grass_dissimilarity.csv")

ggplot(grass_dissimilarity, aes(x = FireFreq, y = BC_dissimilarity, color = Vegetation)) +
  geom_point() + theme_classic() +
  scale_y_continuous(expand = c(0, 0), limits = c(0.01, 1.01)) +
  scale_x_continuous(expand = c(0, 0), limits = c(-0.01, 1.01), breaks = seq(0,1,.2)) +
  scale_color_manual(values=c("green4", "gold"))
# 
# rcorr(as.matrix(grass_dissimilarity[,c(6:9)]),type="pearson")
# 
# grass_dissimilarity$Site <- as.factor(grass_dissimilarity$Site)
# BC_model_grass <- glmmTMB(BC_dissimilarity ~  scale(FireFreq)*Vegetation + scale(FireFreq^2)*Vegetation + scale(MAP) + scale(MAT) + scale(Duration) + (1|Site), grass_dissimilarity, family=beta_family(link="logit"), na.action = "na.fail")
# head(dredge(BC_model_grass))
# BC_grass_best <- glmmTMB(BC_dissimilarity ~  scale(FireFreq) + (1|Site), grass_dissimilarity, family=beta_family(link="logit"), na.action = "na.fail")
# summary(BC_grass_best)
# plot(ggpredict(BC_grass_best, terms = "FireFreq"), residuals = T) # this goes beyond the domain, so have to fix this somehow
# BC_grass_best <- glmmTMB(BC_dissimilarity ~  FireFreq + (1|Site), grass_dissimilarity, family=beta_family(link="logit"), na.action = "na.fail")
# summary(BC_grass_best)


nongrass_dissimilarity <- read.csv("nongrass_dissimilarity.csv")

ggplot(nongrass_dissimilarity, aes(x = FireFreq, y = BC_dissimilarity, color = Vegetation)) +
  geom_point() + theme_classic() +
  scale_y_continuous(expand = c(0, 0), limits = c(0.01, 1.01)) +
  scale_x_continuous(expand = c(0, 0), limits = c(-0.01, 1.01), breaks = seq(0,1,.2)) +
  scale_color_manual(values=c("green4", "gold"))
# 
# 
# # Dissimilarity trees

woody_dissimilarity <- read.csv("woody_dissimilarity.csv")

ggplot(woody_dissimilarity, aes(x = FireFreq, y = BC_dissimilarity, color = Vegetation)) +
  geom_point() + theme_classic() +
  scale_y_continuous(expand = c(0, 0), limits = c(0.01, 1.01)) +
  scale_x_continuous(expand = c(0, 0), limits = c(-0.01, 1.01), breaks = seq(0,1,.2)) +
  scale_color_manual(values=c("green4", "gold"))
# 
# rcorr(as.matrix(ba_dissimilarity[,c(6:9)]),type="pearson")
# 
# ba_dissimilarity <- ba_dissimilarity %>% filter(BC_dissimilarity != 1)
# ba_dissimilarity$Site <- as.factor(ba_dissimilarity$Site)
# BC_model_ba <- glmmTMB(BC_dissimilarity ~  scale(FireFreq)*Vegetation + scale(FireFreq^2)*Vegetation + scale(MAP) + scale(MAT) + (1|Site), ba_dissimilarity, family=beta_family(link="logit"), na.action = "na.fail")
# head(dredge(BC_model_ba), 6)
# BC_ba_best <- glmmTMB(BC_dissimilarity ~  scale(FireFreq)*Vegetation + scale(FireFreq^2) + (1|Site), ba_dissimilarity, family=beta_family(link="logit"), na.action = "na.fail")
# summary(BC_ba_best)
# plot(ggpredict(BC_ba_best, terms = c("FireFreq", "Vegetation")), residuals = T) # this goes beyond the domain, so have to fix this somehow
# 
# BC_ba_best <- glmmTMB(BC_dissimilarity ~  FireFreq*Vegetation + I(FireFreq^2) + (1|Site), ba_dissimilarity, family=beta_family(link="logit"), na.action = "na.fail")
# summary(BC_ba_best)






### Rarefaction

ggplot(responses_herbaceous2, aes(x = cover_adj, y = logRR2, color = Vegetation)) + geom_point() +
  ylab ('log response ratio') + xlab("percent cover") +
  theme_classic() + 
  theme(axis.ticks.length = unit(.3, 'cm'), 
        axis.title = element_text(size = 16), axis.text = element_text(size = 12, color = "black"), 
        axis.ticks = element_line(color = "black"),
        plot.title = element_text(hjust = 0.5),
        legend.position = "none") +
  scale_y_continuous(expand = c(0, 0), limits = c(-1,1), breaks = seq(-1,1,0.5)) +
  scale_x_continuous(expand = c(0, 0), limits = c(0,103), breaks = seq(0, 100, by = 25)) +
  geom_hline(yintercept=0, linetype='dashed') +
  scale_color_manual(values=c("green4", "gold")) +
  coord_cartesian(clip = "off")


rarefaction <- lme(logRR2 ~ cover_adj*Vegetation, random = ~1|Site, data = responses_herbaceous2)
summary(rarefaction)
r.squaredGLMM(rarefaction)

gg_rarefaction <- visreg(lmer(logRR2 ~ cover_adj * Vegetation + (1|Site), data = responses_herbaceous2, na.action = "na.fail"),
                   "cover_adj", by = "Vegetation", line = "Vegetation", gg = T, overlay = T) +  
  geom_point(aes(color=Vegetation, size = I(3)))

gg_rarefaction +   
  ylab ('log response ratio') + xlab("percent cover") +
  theme_classic() + 
  theme(axis.ticks.length = unit(.3, 'cm'), 
        axis.title = element_text(size = 16), axis.text = element_text(size = 12, color = "black"), 
        axis.ticks = element_line(color = "black"),
        plot.title = element_text(hjust = 0.5),
        legend.position = "none") +
  scale_y_continuous(expand = c(0, 0), limits = c(-1.5,1.5), breaks = seq(-1.5,1.5,0.5)) +
  scale_x_continuous(expand = c(0, 0), limits = c(0,103), breaks = seq(0, 100, by = 25)) +
  geom_hline(yintercept=0, linetype='dashed') +
  scale_color_manual(values=c("green4", "gold")) +
  coord_cartesian(clip = "off")




### Herbaceous traits ###

responses_herbaceoustraits <- read.csv("herbaceous_traits.csv")

trait_figure <- function(gg) {
  plot(gg +
         ylab('log response ratio') +
         xlab('fire frequency (fires per year)') +
         theme_classic() + 
         theme(axis.ticks.length = unit(.3, 'cm'), 
               axis.title = element_text(size = 16), axis.text = element_text(size = 12, color = "black"), 
               axis.ticks = element_line(color = "black"),
               plot.title = element_text(hjust = 0.5),
               legend.position = "none") +
         scale_y_continuous(expand = c(0, 0), limits = c(-3,1), breaks = seq(-3,1,0.5)) +
         scale_x_continuous(expand = c(0, 0), limits = c(0,1.015), breaks = seq(0, 1, by = .2)) +
         geom_hline(yintercept=0, linetype='dashed') +
         scale_color_manual(values=c("green4", "gold")) +
         coord_cartesian(clip = "off"))
}

# Seed stats
options(na.action = "na.fail")
herbs_mod_seed <- lmer(logRR_Seed ~ (scale(FireFreq) + scale(FireFreq^2)) * Vegetation + scale(Aridity) + scale(MAT) + scale(Duration) + (1 | Site), REML = F, data = responses_herbaceoustraits)
dredge(herbs_mod_seed)

best_mod_seed <- lme(logRR_Seed ~ scale(FireFreq) + scale(Aridity), random = ~1|Site, data = responses_herbaceoustraits)
summary(best_mod_seed)
r.squaredGLMM(best_mod_seed)
best_mod_seed <- lme(logRR_Seed ~ FireFreq + Aridity, random = ~1|Site, data = responses_herbaceoustraits)
summary(best_mod_seed)

gg_seed <- visreg(lmer(logRR_Seed ~ scale(FireFreq) + scale(Aridity) + (1 | Site), data = responses_herbaceoustraits), "FireFreq", gg = T, line = list(col = "black"), band = T) +
  ggtitle("Seed mass") +
  geom_point(aes(color=responses_herbaceoustraits$Vegetation), size = I(3))
trait_figure(gg_seed)

# Height stats

herbs_mod_height <- lmer(logRR_Height ~ scale(FireFreq) * Vegetation + scale(FireFreq^2) * Vegetation + scale(Aridity) + scale(MAT) + scale(Duration) + (1 | Site), REML = F, data = responses_herbaceoustraits)
dredge(herbs_mod_height)

best_mod_herbheight <- lme(logRR_Height ~ scale(FireFreq) + scale(Aridity), random = ~1|Site, data = responses_herbaceoustraits)
summary(best_mod_herbheight)
r.squaredGLMM(best_mod_herbheight)
best_mod_herbheight <- lme(logRR_Height ~ FireFreq + Aridity, random = ~1|Site, data = responses_herbaceoustraits)
summary(best_mod_herbheight)

gg_heightherb <- visreg(lmer(logRR_Height ~ scale(FireFreq) + scale(Aridity) + (1 | Site), data = responses_herbaceoustraits), "FireFreq", gg = T, line = list(col = "black"), band = T) +
  ggtitle("Maximum height") +
  geom_point(aes(color=responses_herbaceoustraits$Vegetation), size = I(3))
trait_figure(gg_heightherb)

# N stats
herbs_mod_N <- lmer(logRR_N ~ scale(FireFreq) * Vegetation + scale(FireFreq^2) * Vegetation + scale(Aridity) + scale(MAT) + scale(Duration) + (1 | Site), REML = F, data = responses_herbaceoustraits)
dredge(herbs_mod_N)

best_mod_herbN <- lme(logRR_N ~ scale(Aridity), random = ~1|Site, data = responses_herbaceoustraits)
summary(best_mod_herbN)
r.squaredGLMM(best_mod_herbN)
best_mod_herbN <- lme(logRR_N ~ Aridity, random = ~1|Site, data = responses_herbaceoustraits)
summary(best_mod_herbN)

gg_Nherb <- visreg(lmer(logRR_N ~ scale(Aridity) + (1 | Site), data = responses_herbaceoustraits), "Aridity", gg = T, line = list(col = "black"), band = T) +
  ggtitle("Foliar nitrogen") +
  geom_point(aes(color=responses_herbaceoustraits$Vegetation), size = I(3))

gg_Nherb +
  ylab('log response ratio') +
  xlab('aridity index') +
  theme_classic() + 
  theme(axis.ticks.length = unit(.3, 'cm'), 
        axis.title = element_text(size = 16), axis.text = element_text(size = 12, color = "black"), 
        axis.ticks = element_line(color = "black"),
        plot.title = element_text(hjust = 0.5),
        legend.position = "none") +
  scale_y_continuous(expand = c(0, 0), limits = c(-3,1), breaks = seq(-3,1,0.5)) +
  scale_x_continuous(expand = c(0, 0), limits = c(2000,12700), breaks = seq(2000, 12000, by = 2000)) +
  geom_hline(yintercept=0, linetype='dashed') +
  scale_color_manual(values=c("green4", "gold")) +
  coord_cartesian(clip = "off")

# P stats
herbs_mod_P <- lmer(logRR_P ~ scale(FireFreq) * Vegetation + scale(FireFreq^2) * Vegetation + scale(Aridity) + scale(MAT) + scale(Duration) + (1 | Site), REML = F, data = responses_herbaceoustraits)
dredge(herbs_mod_P)

best_mod_herbP <- lme(logRR_P ~ scale(Aridity), random = ~1|Site, data = responses_herbaceoustraits)
summary(best_mod_herbP)
r.squaredGLMM(best_mod_herbP)
best_mod_herbP <- lme(logRR_P ~ Aridity, random = ~1|Site, data = responses_herbaceoustraits)
summary(best_mod_herbP)

gg_Pherb <- visreg(lmer(logRR_P ~ scale(Aridity) + (1 | Site), data = responses_herbaceoustraits), "Aridity", gg = T, line = list(col = "black"), band = T) +
  ggtitle("Foliar phosphorus") +
  geom_point(aes(color=responses_herbaceoustraits$Vegetation), size = I(3))

gg_Pherb +
  ylab('log response ratio') +
  xlab('aridity index') +
  theme_classic() + 
  theme(axis.ticks.length = unit(.3, 'cm'), 
        axis.title = element_text(size = 16), axis.text = element_text(size = 12, color = "black"), 
        axis.ticks = element_line(color = "black"),
        plot.title = element_text(hjust = 0.5),
        legend.position = "none") +
  scale_y_continuous(expand = c(0, 0), limits = c(-3,1), breaks = seq(-3,1,0.5)) +
  scale_x_continuous(expand = c(0, 0), limits = c(2000,12700), breaks = seq(2000, 12000, by = 2000)) +
  geom_hline(yintercept=0, linetype='dashed') +
  scale_color_manual(values=c("green4", "gold")) +
  coord_cartesian(clip = "off")

# Stem stats
responses_herbaceoustraits2 <- responses_herbaceoustraits %>%   filter(!is.na(logRR_Stem))
herbs_mod_stem <- lmer(logRR_Stem ~ scale(FireFreq) * Vegetation + scale(FireFreq^2) * Vegetation + scale(Aridity) + scale(MAT) + scale(Duration) + (1 | Site), REML = F, data = responses_herbaceoustraits2)
dredge(herbs_mod_stem)

ggplot(responses_herbaceoustraits2, aes(x = FireFreq, y = logRR_Stem, color = Vegetation)) +
  geom_point(size = I(3)) +
  ylab('log response ratio') +
  xlab('fire frequency (fires per year)') +
  theme_classic() + 
  theme(axis.ticks.length = unit(.3, 'cm'), 
        axis.title = element_text(size = 16), axis.text = element_text(size = 12),
        plot.title = element_text(hjust = 0.5),
        legend.position = "none") +
  scale_y_continuous(expand = c(0, 0), limits = c(-3,1), breaks = seq(-3,1,0.5)) +
  scale_x_continuous(expand = c(0, 0), limits = c(0,1.01), breaks = seq(0, 1, by = .2)) +
  geom_hline(yintercept=0, linetype='dashed') +
  coord_cartesian(clip= "off") +
  scale_color_manual(values=c("green4", "gold")) +
  ggtitle("Stem diameter")


### Woody Traits ###

responses_woodytraits <- read.csv("woody_traits.csv")

trait_figure_woody <- function(gg) {
  plot(gg +
         ylab('log response ratio') +
         xlab('fire frequency (fires per year)') +
         theme_classic() + 
         theme(axis.ticks.length = unit(.3, 'cm'), 
               axis.title = element_text(size = 16), axis.text = element_text(size = 12, color = "black"), 
               axis.ticks = element_line(color = "black"),
               plot.title = element_text(hjust = 0.5),
               legend.position = "none") +
         scale_y_continuous(expand = c(0, 0), limits = c(-1.5,1.5), breaks = seq(-1.5,1.5,0.5)) +
         scale_x_continuous(expand = c(0, 0), limits = c(0,1.015), breaks = seq(0, 1, by = .2)) +
         geom_hline(yintercept=0, linetype='dashed') +
         scale_color_manual(values=c("green4", "gold"))) +
    coord_cartesian(clip = "off")
}

options(na.action = "na.fail")
rcorr(as.matrix(responses_woodytraits[,c(33, 31:32, 34)]),type="pearson")

# Seed
woody_mod_seed <- lmer(logRR_Seed ~ scale(FireFreq) * Vegetation + scale(FireFreq^2) * Vegetation + scale(Aridity) + scale(MAT) + (1 | Site), REML = F, data = responses_woodytraits)
dredge(woody_mod_seed)

best_mod_woodyseed <- lme(logRR_Seed ~ scale(FireFreq) + scale(FireFreq^2), random = ~1|Site, data = responses_woodytraits)
summary(best_mod_woodyseed)
r.squaredGLMM(best_mod_woodyseed)
best_mod_woodyseed <- lme(logRR_Seed ~ FireFreq + I(FireFreq^2), random = ~1|Site, data = responses_woodytraits)
summary(best_mod_woodyseed)

gg_woodyseed <- visreg(lmer(logRR_Seed ~ scale(FireFreq) + scale(FireFreq^2) + (1 | Site), data = responses_woodytraits), "FireFreq", gg = T, line = list(col = "black"), band = T) +
  ggtitle("Seed mass") +
  geom_point(aes(color=responses_woodytraits$Vegetation), size = I(3))
trait_figure_woody(gg_woodyseed)

# Height
woody_mod_height <- lmer(logRR_Height ~ scale(FireFreq) * Vegetation + scale(FireFreq^2) * Vegetation + scale(Aridity) + scale(MAT) + (1 | Site), REML = F, data = responses_woodytraits)
dredge(woody_mod_height)
best_mod_woodyheight <- lme(logRR_Height ~ scale(FireFreq) + scale(FireFreq^2) + scale(MAT), random = ~1|Site, data = responses_woodytraits)
summary(best_mod_woodyheight)
r.squaredGLMM(best_mod_woodyheight)
best_mod_woodyheight <- lme(logRR_Height ~ FireFreq + I(FireFreq^2) + MAT, random = ~1|Site, data = responses_woodytraits)
summary(best_mod_woodyheight)

gg_woodyheight <- visreg(lmer(logRR_Height ~ scale(FireFreq) + scale(FireFreq^2) + scale(MAT) + (1 | Site), data = responses_woodytraits), "FireFreq", gg = T, line = list(col = "black"), band = T) +
  ggtitle("Maximum height") +
  geom_point(aes(color=responses_woodytraits$Vegetation), size = I(3))
trait_figure_woody(gg_woodyheight)

# N

responses_woodytraits2 <- responses_woodytraits %>% filter(!is.na(logRR_N))
woody_mod_N <- lmer(logRR_N ~ scale(FireFreq) * Vegetation + scale(FireFreq^2) * Vegetation + scale(Aridity) + scale(MAT) + (1 | Site), REML = F, data = responses_woodytraits2)
dredge(woody_mod_N)

best_mod_woodyN <- lme(logRR_N ~ scale(FireFreq) + scale(FireFreq^2) + scale(MAT) + Vegetation, random = ~1|Site, data = responses_woodytraits2)
summary(best_mod_woodyN)
r.squaredGLMM(best_mod_woodyN)
best_mod_woodyN <- lme(logRR_N ~ FireFreq + I(FireFreq^2) + MAT + Vegetation, random = ~1|Site, data = responses_woodytraits2)
summary(best_mod_woodyN)

gg_woodyN <- visreg(lmer(logRR_N ~ scale(FireFreq) + scale(FireFreq^2) + scale(MAT) + Vegetation + (1 | Site), data = responses_woodytraits2), "FireFreq", gg = T, line = list(col = "black"), band = T) +
  ggtitle("Foliar Nitrogen") +
  geom_point(aes(color=responses_woodytraits2$Vegetation), size = I(3))
trait_figure_woody(gg_woodyN)

# LMA
responses_woodytraits3 <- responses_woodytraits %>% filter(!is.na(logRR_LMA))
# removed MAT bc singular
woody_mod_LMA <- lmer(logRR_LMA ~ scale(FireFreq) * Vegetation + scale(FireFreq^2) * Vegetation + scale(Aridity) + (1 | Site), REML = F, data = responses_woodytraits3)
dredge(woody_mod_LMA)

best_mod_LMAwoody <- lme(logRR_LMA ~ scale(FireFreq) * Vegetation + scale(FireFreq^2), random = ~1|Site, data = responses_woodytraits3)
summary(best_mod_LMAwoody)
r.squaredGLMM(best_mod_LMAwoody)

best_mod_LMAwoody <- lme(logRR_LMA ~ FireFreq * Vegetation + I(FireFreq^2), random = ~1|Site, data = responses_woodytraits3)
summary(best_mod_LMAwoody)

gg_woodyLMA <- visreg(lmer(logRR_LMA ~ scale(FireFreq) * Vegetation + scale(FireFreq^2) + scale(MAT) + (1 | Site), data = responses_woodytraits3), "FireFreq", by = "Vegetation", overlay = T, gg = T, band = T) +
  ggtitle("LMA") +
  geom_point(aes(color=Vegetation), size = I(3))
trait_figure_woody(gg_woodyLMA)


# Wood Density

responses_woodytraits4 <- responses_woodytraits %>% filter(!is.na(logRR_Density))
woody_mod_density <- lmer(logRR_Density ~ scale(FireFreq) * Vegetation + scale(FireFreq^2) * Vegetation + scale(Aridity) + scale(MAT) + (1 | Site), REML = F, data = responses_woodytraits4)
dredge(woody_mod_density)

best_mod_woodydensity <- lme(logRR_Density ~ scale(FireFreq) * Vegetation + scale(FireFreq^2) + scale(MAT), random = ~1|Site, data = responses_woodytraits4)
summary(best_mod_woodydensity)
r.squaredGLMM(best_mod_woodydensity)
best_mod_woodydensity <- lme(logRR_Density ~ FireFreq * Vegetation + I(FireFreq^2) + MAT, random = ~1|Site, data = responses_woodytraits4)
summary(best_mod_woodydensity)

gg_woodydensity <- visreg(lmer(logRR_Density ~ scale(FireFreq) * Vegetation + scale(FireFreq^2) + scale(MAT) + (1 | Site), data = responses_woodytraits4), "FireFreq", by = "Vegetation", overlay = T, gg = T, band = T) +
  ggtitle("Wood density") +
  geom_point(aes(color=Vegetation), size = I(3))
trait_figure_woody(gg_woodydensity)
