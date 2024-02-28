setwd("~/Desktop/Research/FireBiodiversity/Data")

packages <- c('dplyr', 'tidyr', 'ggplot2', 'ggpubr', 'ggeffects',
              'Hmisc', 'lme4', 'MuMIn', 'visreg', "nlme")
lapply(packages, require, character.only = T)


### Data processing functions ###

# Calculate plot-level richness (S), Shannon diversity (H), Simpson diversity (D), 
# evenness (J), % represented by top 3 species (dom), number of species that comprise 70% of total (num70)
calcs <- function(cover) {
  df <- cover %>% group_by(Site, Plot, FireFreq, FireStart, FireLast, FireSeason, SurveyYear) %>% 
    arrange(desc(Abundance), .by_group = TRUE) %>% 
    mutate(N = sum(Abundance), p = Abundance/N) %>% 
    summarise(S = n_distinct(GenusSpecies),                                      # species richness
              H = -sum(p*log(p)),                                                # Shannon diversity
              # D = 1 - (sum(Abundance*(Abundance - 1))/(sum(Abundance)*(sum(Abundance)-1))), # Simpson diversity
              J = H/log(S),                                                      # evenness
              # dom = sum(first(p), nth(p, 2), nth(p, 3)),                        # percent represented by 3 most dominant sp
              # num70 = sum(cumsum(p) <= .7)                                      # number of species that comprise 70% of total
              cover = sum(Abundance),
              cover_adjusted = min(cover, 100)
              )                                      
  
  # for plots surveyed multiple years, choose the most recent
  df <- df %>% arrange(Site, Plot, FireFreq, desc(SurveyYear)) %>% 
    distinct(Site, Plot, FireFreq, .keep_all = TRUE)
  
  return(df)
}

# Calculate responses for each plot
# RR <- function(df) {
#   responses <- df %>% filter(FireFreq != 0)
#   controls <- df %>% filter(FireFreq == 0)
#   for (i in 1:nrow(responses)) {
#     for (j in 1:nrow(controls)) {
#       if (responses$Site[i] == controls$Site[j]) {
#         responses$RR_S[i] <- responses$S[i]/controls$S[j]
#         responses$RR_H[i] <- responses$H[i]/controls$H[j]
#         responses$RR_D[i] <- responses$D[i]/controls$D[j]
#         responses$RR_J[i] <- responses$J[i]/controls$J[j]
#         responses$RR_dom[i] <- responses$dom[i]/controls$dom[j]
#         responses$RR_num70[i] <- responses$num70[i]/controls$num70[j]
#       }
#     }
#   }
#   return(responses)
# }

RR2 <- function(df) { 
  responses <- df %>% filter(FireFreq != 0) %>% group_by(Site, FireFreq) %>% 
    summarise(meanS = mean(S), meanH = mean(H), meanJ = mean(J), sdS = sd(S), sdH = sd(H), sdJ = sd(J), n = n(), cover = mean(cover), cover_adj = mean(cover_adjusted)) %>% 
    mutate(CV = sdS/meanS, CV_H = sdH/meanH, CV_J = sdJ/meanJ)
  controls <- df %>% filter(FireFreq == 0) %>% group_by(Site) %>% 
    summarise(meanS = mean(S), meanH = mean(H), meanJ = mean(J), sdS = sd(S), sdH = sd(H), sdJ = sd(J), n = n(), cover = mean(cover), cover_adj = mean(cover_adjusted)) %>% 
    mutate(CV = sdS/meanS, CV_H = sdH/meanH, CV_J = sdJ/meanJ)
  for (i in 1:nrow(responses)) {
    for (j in 1:nrow(controls)) {
      if (responses$Site[i] == controls$Site[j]) {
        responses$RR_S[i] <- responses$meanS[i]/controls$meanS[j]
        responses$logRR[i] <- log(responses$meanS[i]/controls$meanS[j])
        responses$logRR2[i] <- log(responses$meanS[i]/controls$meanS[j]) + 0.5 * ( (((responses$CV[i])^2)/controls$n[j]) - (((controls$CV[j])^2)/responses$n[i]) )
        responses$logRR2_v[i] <- ((responses$CV[i])^2)/responses$n[i] + ((controls$CV[j])^2)/controls$n[j] + ((responses$CV[i])^4)/(2*(responses$n[i])^2) + ((controls$CV[j])^4)/(2*(controls$n[j])^2)
        
        responses$logRR_H[i] <- log(responses$meanH[i]/controls$meanH[j])
        responses$logRR_H2[i] <- log(responses$meanH[i]/controls$meanH[j]) + 0.5 * ( (((responses$CV_H[i])^2)/controls$n[j]) - (((controls$CV_H[j])^2)/responses$n[i]) )
        responses$logRR_H2_v[i] <- ((responses$CV_H[i])^2)/responses$n[i] + ((controls$CV_H[j])^2)/controls$n[j] + ((responses$CV_H[i])^4)/(2*(responses$n[i])^2) + ((controls$CV_H[j])^4)/(2*(controls$n[j])^2)
        
        
        responses$logRR_J[i] <- log(responses$meanJ[i]/controls$meanJ[j])
        responses$logRR_J2[i] <- log(responses$meanJ[i]/controls$meanJ[j]) + 0.5 * ( (((responses$CV_J[i])^2)/controls$n[j]) - (((controls$CV_J[j])^2)/responses$n[i]) )
        responses$logRR_J2_v[i] <- ((responses$CV_J[i])^2)/responses$n[i] + ((controls$CV_J[j])^2)/controls$n[j] + ((responses$CV_J[i])^4)/(2*(responses$n[i])^2) + ((controls$CV_J[j])^4)/(2*(controls$n[j])^2)
        
      }
    }
  }
  return(list(responses, controls))
}

# Calculate response ratios for each fire frequency at each site, and scale based on CV and sample size
RR3 <- function(df) { 
 responses <- df %>% filter(FireFreq != 0) %>% group_by(Site, FireFreq) %>% 
   summarise(meanS = mean(S), sdS = sd(S), n = n()) %>% 
   mutate(CV = sdS/meanS)
  controls <- df %>% filter(FireFreq == 0) %>% group_by(Site) %>% 
    summarise(meanS = mean(S), sdS = sd(S), n = n()) %>% 
    mutate(CV = sdS/meanS)
  for (i in 1:nrow(responses)) {
    for (j in 1:nrow(controls)) {
      if (responses$Site[i] == controls$Site[j]) {
        responses$RR_S[i] <- responses$meanS[i]/controls$meanS[j]
        responses$logRR[i] <- log(responses$meanS[i]/controls$meanS[j])
        responses$logRR2[i] <- log(responses$meanS[i]/controls$meanS[j]) + 0.5 * ( (((responses$CV[i])^2)/controls$n[j]) - (((controls$CV[j])^2)/responses$n[i]) )
        responses$logRR2_v[i] <- ((responses$CV[i])^2)/responses$n[i] + ((controls$CV[j])^2)/controls$n[j] + ((responses$CV[i])^4)/(2*(responses$n[i])^2) + ((controls$CV[j])^4)/(2*(controls$n[j])^2)
      }
    }
  }
  return(list(responses, controls))
}

# Calculate response ratios for sites that do not have plot-level data
RR4 <- function(df) {
  responses <- df %>% filter(FireFreq != 0)
  controls <- df %>% filter(FireFreq == 0)
  for (i in 1:nrow(responses)) {
    for (j in 1:nrow(controls)) {
      responses$logRR[i] <- log(responses$S[i]/controls$S[j])
      responses$logRR2[i] <- log(responses$S[i]/controls$S[j]) + 0.5 * ( (((responses$CV[i])^2)/controls$n[i]) - (((controls$CV[i])^2)/responses$n[i]) )
      responses$logRR2_v[i] <- ((responses$CV[i])^2)/responses$n[i] + ((controls$CV[j])^2)/controls$n[j] + ((responses$CV[i])^4)/(2*(responses$n[i])^2) + ((controls$CV[j])^4)/(2*(controls$n[j])^2)
    }
  }
  return(list(responses, controls))
}

# Match site and climate data data
sitedata <- read.csv("sitedata_full.csv")
match_sitedata <- function(df) {
  df$Vegetation <- as.factor(sitedata$Biome[match(df$Site, sitedata$Site)])
  df$MAP <- as.numeric(sitedata$MAP[match(df$Site, sitedata$Site)])
  df$MAT <- as.numeric(sitedata$MAT[match(df$Site, sitedata$Site)])
  df$Seasonality <- as.numeric(sitedata$Seasonality[match(df$Site, sitedata$Site)])
  # df$WetQ <- as.numeric(sitedata$WetQ[match(df$Site, sitedata$Site)])
  # df$DryQ <- as.numeric(sitedata$DryQ[match(df$Site, sitedata$Site)])
  df$Aridity <- as.numeric(sitedata$aridity[match(df$Site, sitedata$Site)])
  df$Duration <- as.numeric(sitedata$Duration[match(df$Site, sitedata$Site)])
  df$Continent <- as.factor(sitedata$Continent[match(df$Site, sitedata$Site)])
  return(df)
}

# Plot exploratory figures
# figs <- function(responses) {
#   for (var in colnames(responses[,grepl("RR", colnames(responses))])) {
#     assign(paste0("plot", var), 
#            ggplot(responses, aes(x = FireFreq, y = log(.data[[var]]), color = Vegetation)) + 
#              geom_point() + geom_hline(yintercept=0, linetype="dashed") + theme_classic() +
#              scale_y_continuous(expand = c(0, 0), limits = c(-3, 3)) +
#              scale_x_continuous(expand = c(0, 0), limits = c(0, 1.01), breaks = seq(0,1,.2)) +
#              scale_color_manual(values=c("green4", "gold")))
#   }
#   ggarrange(plotRR_S, plotRR_H, plotRR_D, plotRR_J, plotRR_dom, plotRR_num70, nrow = 2, ncol = 3, common.legend = T)
# }

figure <- function(gg) {
  plot(gg +
         ylab('log response ratio') +
         xlab('fire frequency (fires per year)') +
         theme_classic() + 
         theme(axis.ticks.length = unit(.3, 'cm'), 
               axis.title = element_text(size = 16), axis.text = element_text(size = 12, color = "black"), 
               axis.ticks = element_line(color = "black"),
               plot.title = element_text(hjust = 0.5),
               legend.position = "none") +
         scale_y_continuous(expand = c(0, 0), limits = c(-1,1), breaks = seq(-1,1,0.5)) +
         scale_x_continuous(expand = c(0, 0), limits = c(0,1.015), breaks = seq(0, 1, by = .2)) +
         geom_hline(yintercept=0, linetype='dashed') +
         scale_color_manual(values=c("green4", "gold"))) +
        coord_cartesian(clip = "off")
}

supp_figure <- function(df_RR) {
  plot(ggplot(df_RR, aes(x = FireFreq, y = logRR2, color = Vegetation)) +
         geom_point(size = I(3)) +
         ylab('log response ratio') +
         xlab('fire frequency (fires per year)') +
         theme_classic() + 
         theme(axis.ticks.length = unit(.3, 'cm'), 
               axis.title = element_text(size = 16), axis.text = element_text(size = 12),
               plot.title = element_text(hjust = 0.5),
               legend.position = "none") +
         scale_y_continuous(expand = c(0, 0), limits = c(-2,2), breaks = seq(-2,2,0.5)) +
         scale_x_continuous(expand = c(0, 0), limits = c(0,1.01), breaks = seq(0, 1, by = .2)) +
         geom_hline(yintercept=0, linetype='dashed') +
    geom_errorbar(aes(ymin = logRR2 - logRR2_v, ymax = logRR2 + logRR2_v), width = 0.02) +
    coord_cartesian(clip= "off") +
    scale_color_manual(values=c("green4", "gold")))
}

## All species cover ##

cover <- read.csv("percentCoverSurveys.csv")

# cover$SurveyYear[is.na(cover$SurveyYear)] <- "N/A"
# cover$FireStart[is.na(cover$FireStart)] <- "N/A"
# cover$FireLast[is.na(cover$FireLast)] <- "N/A"
# cover$FireSeason[is.na(cover$FireSeason)] <- "N/A"
cover <- cover %>% filter(Site != "JosephJones") %>% filter(!is.na(GenusSpecies)) # remove sites with no plot data and remove all unknown species

exclude <- c("ChimneySpring", "Limestone", "Yenisei", "Konza", "Wharton", "Sequoia", "KingsCanyon", "Escambia", "Lombard", "JosephJones")
cover <- filter(cover, !grepl(paste(exclude, collapse = "|"), Site))

df <- cover %>% calcs() %>% match_sitedata()
# responses_all <- as.data.frame(RR(df))
# responses_all <- match_sitedata(responses_all)
# figs(responses_all)


## Tree basal area ##

# ba <- read.csv("woodyBasalArea.csv")
# ba <- filter(ba, !grepl(paste(exclude, collapse = "|"), Site))
# ba <- ba %>% mutate(Site = replace(Site, Site == "DryPeachester", "Peachester North"),
#                     Site = replace(Site, Site == "WetPeachester", "Peachester South"),
#                     Site = replace(Site, Site == "Kisatchia", "Kisatchie")) %>% match_sitedata()
# 
# df_ba <- ba %>% calcs() %>% match_sitedata()

stems <- read.csv("woodyStemCounts.csv")
stems <- filter(stems, !grepl(paste(exclude, collapse = "|"), Site))
stems <- stems %>% mutate(Site = replace(Site, Site == "DryPeachester", "Peachester North"),
                    Site = replace(Site, Site == "WetPeachester", "Peachester South"),
                    Site = replace(Site, Site == "Kisatchia", "Kisatchie")) %>% match_sitedata()

df_stems <- stems %>% calcs() %>% match_sitedata()

# responses_ba <- as.data.frame(RR(df_ba))
# responses_ba <- match_sitedata(responses_ba)
# figs(responses_ba)


## Herbaceous ##

# separating out herbaceous by subtracting species included in the woody basal area surveys from the percent cover surveys
# woody_sp <- unique(ba$GenusSpecies)
woody_sp2 <- unique(stems$GenusSpecies)
herbaceous <- cover[!(cover$GenusSpecies %in% woody_sp2),]
df_herbaceous <- herbaceous %>% calcs() %>% match_sitedata()

# responses_herbaceous <- as.data.frame(RR(df_herbaceous))
# responses_herbaceous <- match_sitedata(responses_herbaceous)
# figs(responses_herbaceous)


## Grass cover ##

grass <- cover %>% filter(Family == "Poaceae") %>% filter(!is.na(GenusSpecies))
# grass <- grass %>% filter(Site != "ChimneySpring" & Site != "Limestone" & Site != "JosephJones" & Site != "Yenisei" & Site != "Wharton" & Site != "Konza" & Site != "Bauple")

df_grass <- calcs(grass)
df_grass <- match_sitedata(df_grass)

# responses_grass <- as.data.frame(RR(df_grass))
# responses_grass <- match_sitedata(responses_grass)
# figs(responses_grass)


## C3 vs C4 cover ##

# specieslist <- read.csv("grassSpeciesList.csv")
# grass$Pathway <- specieslist$Pathway[match(grass$GenusSpecies, specieslist$Species)]
# C3 <- grass %>% filter(Pathway == "C3")
# C4 <- grass %>% filter(Pathway == "C4")

# df_C3 <- calcs(C3)
# df_C3 <- match_sitedata(df_C3)
# responses_C3 <- as.data.frame(RR(df_C3))
# figs(responses_C3)

# df_C4 <- calcs(C4)
# df_C4 <- match_sitedata(df_C4)
# responses_C4 <- as.data.frame(RR(df_C4))
# figs(responses_C4)


## Additional sites herbaceous richness ##

herbaceous2 <- df_herbaceous %>% RR2()
responses_herbaceous2 <- as.data.frame(herbaceous2[1]) %>% mutate(logRR2 = coalesce(logRR2, logRR),
                                                                  logRR_H2 = coalesce(logRR_H2, logRR_H),
                                                                  logRR_J2 = coalesce(logRR_J2, logRR_J)) %>% 
  match_sitedata()

additional <- read.csv("additionalSites.csv") # import all additional sites
add_herbs <- additional %>% filter(DataType == "HerbaceousRichness") %>% rename(S = Richness) # subset to sites with herbaceous richness data
herbs_indiv <- add_herbs %>% filter(IndivObs  == "Y") %>% RR3() # function for data that are individual observations
responses_herbs_indiv <- as.data.frame(herbs_indiv[1]) %>% mutate(logRR2 = coalesce(logRR2, logRR))
herbs_agg <- add_herbs %>% filter(IndivObs == "N") %>% RR4() # function for data that are aggregate data (with mean and SD)
responses_herbs_agg <- as.data.frame(herbs_agg[1]) %>% mutate(logRR2 = coalesce(logRR2, logRR))

responses_herbs_all <- as.data.frame(rbind(responses_herbaceous2[,c(1,2,9,12,17,18)], responses_herbs_indiv[,c(1,2,5,6,9,10)], responses_herbs_agg[,c(1,3,8,11,15, 16)]) %>% match_sitedata())

# use the square of the weighted average CV for sites with missing SD (in this case, Solar Village)
tmp <- na.omit(responses_herbs_all)
avgCV_herbresponses <- (sum(tmp$n * tmp$CV) / sum(tmp$n))

control_herbs_all <- na.omit(rbind(as.data.frame(herbaceous2[2]) %>% select(Site, n, CV), as.data.frame(herbs_indiv[2]) %>% select(Site, n, CV), as.data.frame(herbs_agg[2]) %>% select(Site, n, CV)))
avgCV_herbcontrols <- (sum(control_herbs_all$n * control_herbs_all$CV) / sum(control_herbs_all$n))

responses_herbs_all[responses_herbs_all$Site == "Solar Village", 5] <- responses_herbs_all[responses_herbs_all$Site == "Solar Village", 5] + .5 * ( (avgCV_herbresponses^2)/5 -  (avgCV_herbcontrols^2)/5)
responses_herbs_all[responses_herbs_all$Site == "Solar Village", 6] <- (avgCV_herbresponses^2)/5 + (avgCV_herbcontrols^2)/5 + (avgCV_herbresponses^4)/(2*5^2) + (avgCV_herbcontrols^2)/(2*5^2)

supp_figure(responses_herbs_all)

## Additional sites grass richness ##

grass2 <- df_grass %>% RR2()
responses_grass2 <- as.data.frame(grass2[1]) %>% mutate(logRR2 = coalesce(logRR2, logRR),
                                                                                      logRR_H2 = coalesce(logRR_H, logRR_H2),
                                                                                      logRR_J2 = coalesce(logRR_J, logRR_J2)) %>% 
  match_sitedata()

add_grass <- additional %>% filter(DataType == "GrassRichness") %>% rename(S = Richness)
add_grass_indiv <- add_grass %>% filter(IndivObs  == "Y") %>% RR3() 
responses_grass_indiv <- as.data.frame(add_grass_indiv[1]) %>% mutate(logRR2 = coalesce(logRR2, logRR))
add_grass_agg <- add_grass %>% filter(IndivObs == "N") %>% RR4() 
responses_grass_agg <- as.data.frame(add_grass_agg[1]) %>% mutate(logRR2 = coalesce(logRR2, logRR))

responses_grass_all <- rbind(responses_grass2[,c(1,2,9,12,17, 18)], responses_grass_indiv[,c(1,2,5,6,9, 10)], responses_grass_agg[,c(1,3,8,11,15, 16)]) %>% match_sitedata()

# use the square of the weighted average CV for sites with missing SD (in this case, Solar Village)
tmp <- na.omit(responses_grass_all)
avgCV_grassresponses <- (sum(tmp$n * tmp$CV) / sum(tmp$n))

control_grass_all <- na.omit(rbind(as.data.frame(grass2[2]) %>% dplyr::select(Site, n, CV), as.data.frame(add_grass_indiv[2]) %>% dplyr::select(Site, n, CV), as.data.frame(add_grass_agg[2]) %>% dplyr::select(Site, n, CV)))
avgCV_grasscontrols <- (sum(control_grass_all$n * control_grass_all$CV) / sum(control_grass_all$n))

responses_grass_all[responses_grass_all$Site == "Solar Village", 5] <- responses_grass_all[responses_grass_all$Site == "Solar Village", 5] + .5 * ( (avgCV_grassresponses^2)/5 -  (avgCV_grasscontrols^2)/5)
responses_grass_all[responses_grass_all$Site == "Solar Village", 6] <- (avgCV_grassresponses^2)/5 + (avgCV_grasscontrols^2)/5 + (avgCV_grassresponses^4)/(2*5^2) + (avgCV_grasscontrols^2)/(2*5^2)

supp_figure(responses_grass_all)

## Non-grass ##
nongrass <- herbaceous %>% filter(Family != "Poaceae") %>% filter(!is.na(GenusSpecies))
df_nongrass <- nongrass %>% calcs() %>% match_sitedata()
nongrass2 <- RR2(df_nongrass)

responses_nongrass2 <- as.data.frame(nongrass2[1]) %>% mutate(logRR2 = coalesce(logRR2, logRR),
                                                        logRR_H2 = coalesce(logRR_H, logRR_H2),
                                                        logRR_J2 = coalesce(logRR_J, logRR_J2)) %>% 
  match_sitedata()

add_nongrass <- inner_join(add_herbs, add_grass, by = join_by(Site, Plot, FireFreq, IndivObs, Duration, n)) %>% 
  mutate(S = S.x - S.y, Variance = Variance.x + Variance.y, SD = sqrt(Variance), CV = SD/S) %>% 
  dplyr::select(-c(Notes.x, Notes.y)) 
  
add_nongrass_indiv <- add_nongrass %>% filter(IndivObs  == "Y") %>% RR3() 
responses_nongrass_indiv <- as.data.frame(add_nongrass_indiv[1]) %>% mutate(logRR2 = coalesce(logRR2, logRR))
add_nongrass_agg <- add_nongrass %>% filter(IndivObs == "N") %>% RR4() 
responses_nongrass_agg <- as.data.frame(add_nongrass_agg[1]) %>% mutate(logRR2 = coalesce(logRR2, logRR))

responses_nongrass_all <- rbind(responses_nongrass2[,c(1,2,9,12,17, 18)], responses_nongrass_indiv[,c(1,2,5,6,9, 10)], responses_nongrass_agg[,c(1,3,8,22, 24,25)]) %>% match_sitedata()
responses_nongrass_all[which(responses_nongrass_all$Site == "Nylsvley"),5] <- 0

tmp <- na.omit(responses_nongrass_all)
avgCV_nongrassresponses <- (sum(tmp$n * tmp$CV) / sum(tmp$n))

control_nongrass_all <- na.omit(rbind(as.data.frame(nongrass2[2]) %>% dplyr::select(Site, n, CV), as.data.frame(add_nongrass_indiv[2]) %>% dplyr::select(Site, n, CV), as.data.frame(add_nongrass_agg[2]) %>% dplyr::select(Site, n, CV)))
avgCV_nongrasscontrols <- (sum(control_nongrass_all$n * control_nongrass_all$CV) / sum(control_nongrass_all$n))

responses_nongrass_all[responses_nongrass_all$Site == "Solar Village", 5] <- responses_nongrass_all[responses_nongrass_all$Site == "Solar Village", 5] + .5 * ( (avgCV_nongrassresponses^2)/5 -  (avgCV_nongrasscontrols^2)/5)
responses_nongrass_all[responses_nongrass_all$Site == "Solar Village", 6] <- (avgCV_nongrassresponses^2)/5 + (avgCV_nongrasscontrols^2)/5 + (avgCV_nongrassresponses^4)/(2*5^2) + (avgCV_nongrasscontrols^2)/(2*5^2)

supp_figure(responses_nongrass_all)

## Additional sites woody richness ##

# woody2 <- df_ba %>% RR2()
woody2 <- df_stems %>% RR2()
responses_woody2 <- as.data.frame(woody2[1]) %>% mutate(logRR2 = coalesce(logRR2, logRR),
                                                        logRR_H2 = coalesce(logRR_H, logRR_H2),
                                                        logRR_J2 = coalesce(logRR_J, logRR_J2)) %>% 
  match_sitedata()

add_woody <- additional %>% filter(DataType == "WoodyRichness") %>% rename(S = Richness)
add_woody_indiv <- add_woody %>% filter(IndivObs  == "Y") %>% RR3()
responses_woody_indiv <- as.data.frame(add_woody_indiv[1]) %>% mutate(logRR2 = coalesce(logRR2, logRR))
add_woody_agg <- add_woody %>% filter(IndivObs == "N") %>% RR4() 
responses_woody_agg <- as.data.frame(add_woody_agg[1]) %>% mutate(logRR2 = coalesce(logRR2, logRR))

responses_woody_all <- rbind(responses_woody2[,c(1,2,9,12,17,18)], responses_woody_indiv[,c(1,2,5,6,9,10)], responses_woody_agg[,c(1,3,8,11,15,16)]) %>% match_sitedata()

# use the square of the weighted average CV for sites with missing SD (in this case, Solar Village)
tmp <- na.omit(responses_woody_all)
avgCV_woodyresponses <- (sum(tmp$n * tmp$CV) / sum(tmp$n))

control_woody_all <- na.omit(rbind(as.data.frame(woody2[2]) %>% select(Site, n, CV), as.data.frame(add_woody_indiv[2]) %>% select(Site, n, CV), as.data.frame(add_woody_agg[2]) %>% select(Site, n, CV)))
avgCV_woodycontrols <- (sum(control_woody_all$n * control_woody_all$CV) / sum(control_woody_all$n))

responses_woody_all[responses_woody_all$Site == "Solar Village", 5] <- responses_woody_all[responses_woody_all$Site == "Solar Village", 5] + .5 * ( (avgCV_woodyresponses^2)/5 -  (avgCV_woodycontrols^2)/5)
responses_woody_all[responses_woody_all$Site == "Solar Village", 6] <- (avgCV_woodyresponses^2)/5 + (avgCV_woodycontrols^2)/5 + (avgCV_woodyresponses^4)/(2*5^2) + (avgCV_woodycontrols^2)/(2*5^2)

supp_figure(responses_woody_all)



# All

# df_both <- inner_join(df_herbaceous, df_stems, by = join_by(Site, Plot, FireFreq, FireSeason)) %>% 
#   mutate(S = S.x + S.y) 
# 
# df_both <- inner_join(df_herbaceous, df_ba, by = join_by(Site, Plot, FireFreq, FireSeason)) %>% 
#   mutate(S = S.x + S.y) 
# 
# responses_both <- df_both %>% RR3() %>% match_sitedata()
# responses_both <- as.data.frame(responses_both[1]) %>% mutate(logRR2 = coalesce(logRR2, logRR))



### Calculating community (dis)similarity ###

# Note: The goal here is (within each site) to compare each plot to every unburned plot.
# This code is messy but I'm doing it this way because each site has a different number of unburned plots

dissimilarity <- function(cover) {
  cover$ID <- paste(cover$Plot, cover$FireFreq, sep = ".") # Create IDs that combine plot name and fire frequency (because some sites have plots with same name but different fire frequency)

  S_frame <- data.frame(matrix(ncol = 2)) # Create a dataframe to keep track of the S values (total number of specimens) for each unburned site
  colnames(S_frame) <- c('name', 'S')
  BC_frame <- data.frame(matrix(ncol = 3))
  colnames(BC_frame) <- c('Site', 'FireFreq', 'BC') # Create a dataframe for the Bray-Curtis Dissimilarity values

  for(site in unique(cover$Site)) { # Loop through all sites

    x <- 0 # Reset no burn plot counter
    tmp_df <- cover %>% filter(Site == site) # Create a temporary dataframe for the current site
    site_ids <- unique(cover$ID[cover$Site == site]) # List the unique plot IDs (plot # + FireFreq)

    for(id in site_ids) { # Loop through all plot IDs

      if(unique(cover$FireFreq[cover$ID == id] == 0)) { # Enter this section if the plot is unburned
        x <- x + 1 # Add to counter of no burn plots

        tmp_df <- tmp_df %>%
          group_by(GenusSpecies) %>%
          mutate("temp.{x}" := max(0, Abundance[ID == id])) %>% # Create a temporary variable for the species value for the unburned site
          group_by(ID, GenusSpecies) %>%
          mutate("min.{x}" := min(Abundance, get(paste("temp",x,sep = ".")))) # Compare the abundance of the species in each plot to the abundance of the same species in the unburned plot and take the minimum

        S_frame <- S_frame %>% add_row(name = paste("S", site, x, sep = "."), S = sum(tmp_df$Abundance[tmp_df$ID == id])) # Add to S value (total specimen) dataframe for unburned plots
      }
    }
    S_frame <- na.omit(S_frame)

    for (i in (1:x)) { # Loop through the number of unburned plots to calculate dissimilarity of each plot to each each unburned plot
      tmp_df2 <- tmp_df %>%
        group_by(Site, ID, FireFreq) %>%
        summarise(S = sum(Abundance), C = sum(!!rlang::sym(paste("min", i, sep = ".")))) %>% # Calculate the number of specimens at site i (S) and the sum of the lesser values for the species found at site i
        mutate(BC = 1 - ((2*C)/((S_frame[S_frame$name == paste("S", site, i, sep = "."),]$S)+S))) # Calculate Bray Curtis Dissimilarity
      BC_frame <- rbind(BC_frame, tmp_df2[,c(1,3,6)]) # Add to data frame
    }
  }

  BC_frame <- BC_frame %>% filter(BC != 0) # remove values that are exactly the same (i.e. the unburned plots when compared to themselves)
  BC_frame2 <- BC_frame %>% group_by(Site, FireFreq) %>% summarise(BC_dissimilarity = mean(BC), SD = sd(BC)) # Get averages for each fire frequency at each site

  BC_frame2 <- BC_frame2 %>% match_sitedata() %>% drop_na(Site)

  plot(ggplot(BC_frame2, aes(x = FireFreq, y = BC_dissimilarity, color = Vegetation)) +
    geom_point() + theme_classic() +
    scale_y_continuous(expand = c(0, 0), limits = c(0.01, 1.01)) +
    scale_x_continuous(expand = c(0, 0), limits = c(-0.01, 1.01), breaks = seq(0,1,.2)) +
    scale_color_manual(values=c("green4", "gold")))

  return(BC_frame2)
}

herbaceous_dissimilarity <- dissimilarity(herbaceous)
grass_dissimilarity <- dissimilarity(grass)
nongrass_dissimilarity <- dissimilarity(nongrass)
stem_dissimilarity <- dissimilarity(stems)


### Statistical Analysis ###

## All sites herbaceous species richness ##

# correlation matrix
rcorr(as.matrix(responses_herbs_all[,c(9:12)]),type="pearson")

# full model and model selection
herbs_mod_all <- lmer(logRR2 ~ (scale(FireFreq) + scale(FireFreq^2) + scale(Aridity) + scale(Seasonality) + scale(Duration))*Vegetation + (1|Site), data = responses_herbs_all, REML = F, na.action = "na.fail")
head(dredge(herbs_mod_all))

# best model and summary
best_herbs_all <- lme(logRR2 ~ scale(FireFreq) + scale(FireFreq^2), random = ~1|Site, data = responses_herbs_all)
summary(best_herbs_all)
r.squaredGLMM(best_herbs_all)

# model with unscaled predictors
best_herbs_all <- lme(logRR2 ~ FireFreq + I(FireFreq^2), random = ~1|Site, data = responses_herbs_all)
summary(best_herbs_all)

examples <- data.frame(FireFreq= c(0.25, 0.5, 1), 
                       Site = rep("Cedar Creek"), 3)

pct_change <- 100 * (exp(predict(best_herbs_all, examples))-1)

# visualization
gg_herbaceous <- visreg(lmer(logRR2 ~ scale(FireFreq) + scale(FireFreq^2) + (1|Site), data = responses_herbs_all, na.action = "na.fail"),
                        "FireFreq", gg = T, line = list(col = "black")) + 
  ggtitle("Herbaceous vegetation") +
  geom_point(aes(color=responses_herbs_all$Vegetation, size = I(3)))
figure(gg_herbaceous)

## All sites grass species richness ##

# correlation matrix
rcorr(as.matrix(responses_grass_all[,c(8:11)]),type="pearson")

# full model and model selection
grass_mod_all <- lmer(logRR2 ~ (scale(FireFreq) + scale(FireFreq^2) + scale(MAP) + scale(Seasonality) + scale(Duration))*Vegetation + (1|Site), data = responses_grass_all, REML = F, na.action = "na.fail")
head(dredge(grass_mod_all))

# best model and summary
best_grass_all <- lme(logRR2 ~ scale(FireFreq) + scale(FireFreq^2) * Vegetation, random = ~1|Site, data = responses_grass_all)
summary(best_grass_all)
r.squaredGLMM(best_grass_all)

# model with unscaled predictors
best_grass_all <- lme(logRR2 ~ FireFreq + I(FireFreq^2) * Vegetation, random = ~1|Site, data = responses_grass_all)
summary(best_grass_all)

examples <- data.frame(FireFreq= c(0.25, 0.5, 1, 0.25, 0.5, 1), 
                       Vegetation = as.factor(c("Savanna", "Savanna", "Savanna", "Forest", "Forest", "Forest")),
                       Site = rep("Cedar Creek"), 6)

pct_change <- 100 * (exp(predict(best_grass_all, examples))-1)

# visualization
gg_grass <- visreg(lmer(logRR2 ~ FireFreq + I(FireFreq^2) * Vegetation + (1|Site), data = responses_grass_all, na.action = "na.fail"),
                   "FireFreq", by = "Vegetation", overlay = T, gg = T) +  
    ggtitle("Grasses") +
    geom_point(aes(color=Vegetation, size = I(3)))
figure(gg_grass)

## All sites nongrass species richness

# correlation matrix
rcorr(as.matrix(responses_nongrass_all[,c(8:11)]),type="pearson")

# full model and model selection
nongrass_mod_all <- lmer(logRR2 ~ (scale(FireFreq) + scale(FireFreq^2) + scale(MAP) + scale(Seasonality) + scale(Duration))*Vegetation + (1|Site), data = responses_nongrass_all, REML = F, na.action = "na.fail")
head(dredge(nongrass_mod_all))

# best model and summary
best_nongrass_all <- lme(logRR2 ~ scale(FireFreq) + scale(FireFreq^2), random = ~1|Site, data = responses_nongrass_all)
summary(best_nongrass_all)
r.squaredGLMM(best_nongrass_all)

# model with unscaled predictors
best_nongrass_all <- lme(logRR2 ~ FireFreq + I(FireFreq^2), random = ~1|Site, data = responses_nongrass_all)
summary(best_nongrass_all)

gg_nongrass <- visreg(lmer(logRR2 ~ FireFreq + I(FireFreq^2) + (1|Site), data = responses_nongrass_all, na.action = "na.fail"),
                   "FireFreq", overlay = T, gg = T, line = list(col = "black")) +  
  ggtitle("Non-grass herbaceous vegetation") +
  geom_point(aes(color=responses_nongrass_all$Vegetation, size = I(3)))
figure(gg_nongrass)


## All sites woody species richness ##

# correlation matrix
rcorr(as.matrix(responses_woody_all[,c(8:11)]),type="pearson")

# full model and model selection

woody_mod_all <- lmer(logRR2 ~ (scale(FireFreq) + scale(FireFreq^2) + scale(MAP) + scale(Seasonality))*Vegetation + (1|Site), data = responses_woody_all, REML = F, na.action = "na.fail")
head(dredge(woody_mod_all))

# best model and summary
best_woody_all <- lme(logRR2 ~ scale(FireFreq^2) + (scale(FireFreq) + scale(Seasonality)) * Vegetation, random = ~1|Site, data = responses_woody_all)
summary(best_woody_all)
r.squaredGLMM(best_woody_all)

# model with unscaled predictors
best_woody_all <- lme(logRR2 ~ scale(FireFreq^2) + (scale(FireFreq) + scale(Seasonality)) * Vegetation, random = ~1|Site, data = responses_woody_all)
summary(best_woody_all)

examples <- data.frame(FireFreq= c(0.25, 0.5, 1, 0.25, 0.5, 1), 
                   Seasonality = rep(median(responses_woody_all$Seasonality), 6),
                   Vegetation = as.factor(c("Savanna", "Savanna", "Savanna", "Forest", "Forest", "Forest")),
                   Site = rep("Cedar Creek"), 6)

pct_change <- 100 * (exp(predict(best_woody_all, examples))-1)


gg_woody <- visreg(lmer(logRR2 ~ scale(FireFreq^2) + (scale(FireFreq) + scale(Seasonality)) * Vegetation + (1|Site), data = responses_woody_all, na.action = "na.fail"),
                   "FireFreq", by = "Vegetation", overlay = T, gg = T) +  
  ggtitle("Woody vegetation") +
  geom_point(aes(color=Vegetation, size = I(3)))
figure(gg_woody)

visreg(lmer(logRR2 ~ scale(FireFreq^2) + (scale(FireFreq) + scale(Seasonality)) * Vegetation + (1|Site), data = responses_woody_all, na.action = "na.fail"),
       "Seasonality", by = "Vegetation", gg = T)

## Plot-level herbaceous models ##

# correlation matrix 
rcorr(as.matrix(responses_herbaceous2[,c(24:27)]),type="pearson")

# Species richness # 
herbs_mod_richness <- lmer(logRR2 ~ scale(FireFreq) * Vegetation + scale(FireFreq^2) * Vegetation + scale(MAP) + scale(MAT) + scale(Duration) + (1|Site), data = responses_herbaceous2, REML = F, na.action = "na.fail")
head(dredge(herbs_mod_richness))

best_mod_richness <- lme(logRR2 ~ scale(FireFreq) + scale(FireFreq^2), random = ~1|Site, data = responses_herbaceous2)
summary(best_mod_richness)
r.squaredGLMM(best_mod_richness)

best_mod_richness <- lme(logRR2 ~ FireFreq + I(FireFreq^2), random = ~1|Site, data = responses_herbaceous2)
summary(best_mod_richness)

gg_richness <- visreg(lmer(logRR2 ~ scale(FireFreq) + scale(FireFreq^2) + (1|Site), data = responses_herbaceous2, REML = F, na.action = "na.fail"), 
                      "FireFreq", gg = T, line = list(col = "black")) +
  ggtitle("Species richness") +
  geom_point(aes(color=responses_herbaceous2$Vegetation), size = I(3)) 
figure(gg_richness)


# Shannon diversity # 

herbs_mod_shannon <- lmer(logRR_H2 ~ (scale(FireFreq) + scale(FireFreq^2) + scale(MAP) + scale(MAT) + scale(Duration))* Vegetation + (1|Site), data = responses_herbaceous2, REML = F, na.action = "na.fail")
head(dredge(herbs_mod_shannon))

best_mod_shannon <- lme(logRR_H2 ~ scale(FireFreq) + scale(FireFreq^2), random = ~1|Site, data = responses_herbaceous2)
summary(best_mod_shannon)
r.squaredGLMM(best_mod_shannon)

best_mod_shannon <- lme(logRR_H2 ~ FireFreq + I(FireFreq^2), random = ~1|Site, data = responses_herbaceous2)
summary(best_mod_shannon)

gg_shannon <- visreg(lmer(logRR_H2 ~ scale(FireFreq) + scale(FireFreq^2) + (1|Site), data = responses_herbaceous2, REML = F, na.action = "na.fail"),
                     "FireFreq", gg = T, line = list(col = "black")) +  
  ggtitle ("Shannon diversity") +
  geom_point(aes(color=responses_herbaceous2$Vegetation), size = I(3))
figure(gg_shannon)

# Evenness

herbs_mod_evenness <- lmer(logRR_J2 ~ scale(FireFreq) * Vegetation + scale(FireFreq^2) * Vegetation + scale(MAP) + scale(MAT) + scale(Duration) + (1|Site), data = responses_herbaceous2, REML = F, na.action = "na.fail")
head(dredge(herbs_mod_evenness))


best_mod_evenness <- lme(logRR_J2 ~ scale(FireFreq) + scale(FireFreq^2), random = ~1|Site, data = responses_herbaceous2)
summary(best_mod_evenness)
r.squaredGLMM(best_mod_evenness)

best_mod_evenness <- lme(logRR_J2 ~ FireFreq + I(FireFreq^2), random = ~1|Site, data = responses_herbaceous2)
summary(best_mod_evenness)

gg_evenness <- visreg(lmer(logRR_J2 ~ scale(FireFreq) + scale(FireFreq^2) + scale(MAP) + (1|Site), data = responses_herbaceous2, REML = F, na.action = "na.fail"), "FireFreq", gg = T, line = list(col = "black")) +  
  ggtitle("Evenness") +
  geom_point(aes(color=responses_herbaceous2$Vegetation), size = I(3))
figure(gg_evenness)

# Community dissimilarity herbaceous

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
# 
# 
# # Dissimilarity trees
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

herb_traits <- read.csv("Final_Herbaceous_Trait_List.csv")

seeds <- herbaceous %>% group_by(Site, Plot, FireFreq, FireStart, FireLast, FireSeason, SurveyYear) %>% 
  mutate(seed = herb_traits$DrySeedStdValue[match(GenusSpecies, herb_traits$Species)]) %>% 
  filter(!is.na(seed)) %>% 
  mutate(N = sum(Abundance), p = Abundance/N, seed_frac = p * seed) %>% 
  summarise(seed_cwm = mean(seed_frac))

heights <- herbaceous %>% group_by(Site, Plot, FireFreq, FireStart, FireLast, FireSeason, SurveyYear) %>% 
  mutate(height = herb_traits$HeightStdValue[match(GenusSpecies, herb_traits$Species)]) %>% 
  filter(!is.na(height)) %>% 
  mutate(N = sum(Abundance), p = Abundance/N, height_frac = p * height) %>% 
  summarise(height_cwm = mean(height_frac))

nitrogen <- herbaceous %>% group_by(Site, Plot, FireFreq, FireStart, FireLast, FireSeason, SurveyYear) %>% 
  mutate(nitrogen = herb_traits$NitrogenStdValue[match(GenusSpecies, herb_traits$Species)]) %>% 
  filter(!is.na(nitrogen)) %>% 
  mutate(N = sum(Abundance), p = Abundance/N, nitrogen_frac = p * nitrogen) %>% 
  summarise(nitrogen_cwm = mean(nitrogen_frac))

phosphorus <- herbaceous %>% group_by(Site, Plot, FireFreq, FireStart, FireLast, FireSeason, SurveyYear) %>% 
  mutate(phosphorus = herb_traits$PhosphorusStdValue[match(GenusSpecies, herb_traits$Species)]) %>% 
  filter(!is.na(phosphorus)) %>% 
  mutate(N = sum(Abundance), p = Abundance/N, phosphorus_frac = p * phosphorus) %>% 
  summarise(phosphorus_cwm = mean(phosphorus_frac))

stems_trait <- herbaceous %>% group_by(Site, Plot, FireFreq, FireStart, FireLast, FireSeason, SurveyYear) %>% 
  mutate(stem = herb_traits$StemStdValue[match(GenusSpecies, herb_traits$Species)]) %>% 
  filter(!is.na(stem)) %>% 
  mutate(N = sum(Abundance), p = Abundance/N, stem_frac = p * stem) %>% 
  summarise(stem_cwm = mean(stem_frac))

library(tidyverse)
cwm_list <- list(seeds, heights, nitrogen, phosphorus, stems_trait)
cwm_herbaceous <- cwm_list %>% reduce(left_join, by = c("Site" = "Site", "Plot" = "Plot", "FireFreq" = "FireFreq", "FireStart" = "FireStart", "FireLast" = "FireLast", "FireSeason" = "FireSeason", "SurveyYear" = "SurveyYear"))

cwm_herbaceous <- cwm_herbaceous  %>% arrange(Site, Plot, FireFreq, desc(SurveyYear)) %>% 
  distinct(Site, Plot, FireFreq, .keep_all = TRUE)

responses_herbaceoustraits <- cwm_herbaceous %>% filter(FireFreq != 0) %>% group_by(Site, FireFreq) %>% 
  summarise(meanSeed = mean(seed_cwm), meanHeight = mean(height_cwm), meanN = mean(nitrogen_cwm), meanP = mean(phosphorus_cwm), meanStem = mean(stem_cwm),
            CV_Seed = (sd(seed_cwm))/meanSeed , CV_Height = (sd(height_cwm)) / meanHeight, CV_N = (sd(nitrogen_cwm)) / meanN, CV_P = (sd(phosphorus_cwm))/meanP, CV_Stem = (sd(stem_cwm)) / meanStem,
            n = n())
controls_herbaceoustraits <- cwm_herbaceous %>% filter(FireFreq == 0) %>% group_by(Site) %>% 
  summarise(meanSeed = mean(seed_cwm), meanHeight = mean(height_cwm), meanN = mean(nitrogen_cwm), meanP = mean(phosphorus_cwm), meanStem = mean(stem_cwm),
            CV_Seed = (sd(seed_cwm))/meanSeed , CV_Height = (sd(height_cwm)) / meanHeight, CV_N = (sd(nitrogen_cwm)) / meanN, CV_P = (sd(phosphorus_cwm))/meanP, CV_Stem = (sd(stem_cwm)) / meanStem,
            n = n())


for (i in 1:nrow(responses_herbaceoustraits)) {
  for (j in 1:nrow(controls_herbaceoustraits)) {
    if (responses_herbaceoustraits$Site[i] == controls_herbaceoustraits$Site[j]) {
      
      responses_herbaceoustraits$logRR_Seed1[i] <- log(responses_herbaceoustraits$meanSeed[i]/controls_herbaceoustraits$meanSeed[j])
      responses_herbaceoustraits$logRR_Seed[i] <- log(responses_herbaceoustraits$meanSeed[i]/controls_herbaceoustraits$meanSeed[j]) + 0.5 * ( (((responses_herbaceoustraits$CV_Seed[i])^2)/controls_herbaceoustraits$n[j]) - (((controls_herbaceoustraits$CV_Seed[j])^2)/responses_herbaceoustraits$n[i]) )
      responses_herbaceoustraits$logRR_v_Seed[i] <- ((responses_herbaceoustraits$CV_Seed[i])^2)/responses_herbaceoustraits$n[i] + ((controls_herbaceoustraits$CV_Seed[j])^2)/controls_herbaceoustraits$n[j] + ((responses_herbaceoustraits$CV_Seed[i])^4)/(2*(responses_herbaceoustraits$n[i])^2) + ((controls_herbaceoustraits$CV_Seed[j])^4)/(2*(controls_herbaceoustraits$n[j])^2)
      
      responses_herbaceoustraits$logRR_Height1[i] <- log(responses_herbaceoustraits$meanHeight[i]/controls_herbaceoustraits$meanHeight[j])
      responses_herbaceoustraits$logRR_Height[i] <- log(responses_herbaceoustraits$meanHeight[i]/controls_herbaceoustraits$meanHeight[j]) + 0.5 * ( (((responses_herbaceoustraits$CV_Height[i])^2)/controls_herbaceoustraits$n[j]) - (((controls_herbaceoustraits$CV_Height[j])^2)/responses_herbaceoustraits$n[i]) )
      responses_herbaceoustraits$logRR_v_Height[i] <- ((responses_herbaceoustraits$CV_Height[i])^2)/responses_herbaceoustraits$n[i] + ((controls_herbaceoustraits$CV_Height[j])^2)/controls_herbaceoustraits$n[j] + ((responses_herbaceoustraits$CV_Height[i])^4)/(2*(responses_herbaceoustraits$n[i])^2) + ((controls_herbaceoustraits$CV_Height[j])^4)/(2*(controls_herbaceoustraits$n[j])^2)
      
      responses_herbaceoustraits$logRR_N1[i] <- log(responses_herbaceoustraits$meanN[i]/controls_herbaceoustraits$meanN[j])
      responses_herbaceoustraits$logRR_N[i] <- log(responses_herbaceoustraits$meanN[i]/controls_herbaceoustraits$meanN[j]) + 0.5 * ( (((responses_herbaceoustraits$CV_N[i])^2)/controls_herbaceoustraits$n[j]) - (((controls_herbaceoustraits$CV_N[j])^2)/responses_herbaceoustraits$n[i]) )
      responses_herbaceoustraits$logRR_v_N[i] <- ((responses_herbaceoustraits$CV_N[i])^2)/responses_herbaceoustraits$n[i] + ((controls_herbaceoustraits$CV_N[j])^2)/controls_herbaceoustraits$n[j] + ((responses_herbaceoustraits$CV_N[i])^4)/(2*(responses_herbaceoustraits$n[i])^2) + ((controls_herbaceoustraits$CV_N[j])^4)/(2*(controls_herbaceoustraits$n[j])^2)
      
      responses_herbaceoustraits$logRR_P1[i] <- log(responses_herbaceoustraits$meanP[i]/controls_herbaceoustraits$meanP[j])
      responses_herbaceoustraits$logRR_P[i] <- log(responses_herbaceoustraits$meanP[i]/controls_herbaceoustraits$meanP[j]) + 0.5 * ( (((responses_herbaceoustraits$CV_P[i])^2)/controls_herbaceoustraits$n[j]) - (((controls_herbaceoustraits$CV_P[j])^2)/responses_herbaceoustraits$n[i]) )
      responses_herbaceoustraits$logRR_v_P[i] <- ((responses_herbaceoustraits$CV_P[i])^2)/responses_herbaceoustraits$n[i] + ((controls_herbaceoustraits$CV_P[j])^2)/controls_herbaceoustraits$n[j] + ((responses_herbaceoustraits$CV_P[i])^4)/(2*(responses_herbaceoustraits$n[i])^2) + ((controls_herbaceoustraits$CV_P[j])^4)/(2*(controls_herbaceoustraits$n[j])^2)
      
      responses_herbaceoustraits$logRR_Stem1[i] <- log(responses_herbaceoustraits$meanStem[i]/controls_herbaceoustraits$meanStem[j])
      responses_herbaceoustraits$logRR_Stem[i] <- log(responses_herbaceoustraits$meanStem[i]/controls_herbaceoustraits$meanStem[j]) + 0.5 * ( (((responses_herbaceoustraits$CV_Stem[i])^2)/controls_herbaceoustraits$n[j]) - (((controls_herbaceoustraits$CV_Stem[j])^2)/responses_herbaceoustraits$n[i]) )
      responses_herbaceoustraits$logRR_v_Stem[i] <- ((responses_herbaceoustraits$CV_Stem[i])^2)/responses_herbaceoustraits$n[i] + ((controls_herbaceoustraits$CV_Stem[j])^2)/controls_herbaceoustraits$n[j] + ((responses_herbaceoustraits$CV_Stem[i])^4)/(2*(responses_herbaceoustraits$n[i])^2) + ((controls_herbaceoustraits$CV_Stem[j])^4)/(2*(controls_herbaceoustraits$n[j])^2)
      
    }
  }
}

responses_herbaceoustraits <- responses_herbaceoustraits %>% mutate(logRR_Seed = coalesce(logRR_Seed, logRR_Seed1),
                                                                    logRR_Height = coalesce(logRR_Height, logRR_Height1),
                                                                    logRR_N = coalesce(logRR_N, logRR_N1),
                                                                    logRR_P = coalesce(logRR_P, logRR_P1),
                                                                    logRR_Stem = coalesce(logRR_Stem, logRR_Stem1)) %>% 
  match_sitedata()

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

rcorr(as.matrix(responses_herbaceoustraits[,c(30:33)]),type="pearson")

# Seed stats
options(na.action = "na.fail")
herbs_mod_seed <- lmer(logRR_Seed ~ scale(FireFreq) * Vegetation + scale(FireFreq^2) * Vegetation + scale(MAP) + scale(MAT) + scale(Duration) + (1 | Site), REML = F, data = responses_herbaceoustraits)
head(dredge(herbs_mod_seed))
best_mod_seed <- lme(logRR_Seed ~ scale(FireFreq) + scale(MAP), random = ~1|Site, data = responses_herbaceoustraits)
summary(best_mod_seed)
r.squaredGLMM(best_mod_seed)
best_mod_seed <- lme(logRR_Seed ~ FireFreq + MAP, random = ~1|Site, data = responses_herbaceoustraits)
summary(best_mod_seed)

gg_seed <- visreg(lmer(logRR_Seed ~ scale(FireFreq) + scale(MAP) + (1 | Site), data = responses_herbaceoustraits), "FireFreq", gg = T, line = list(col = "black"), band = F) +
  ggtitle("Seed mass") +
  geom_point(aes(color=responses_herbaceoustraits$Vegetation), size = I(3))
trait_figure(gg_seed)

# Height stats

herbs_mod_height <- lmer(logRR_Height ~ scale(FireFreq) * Vegetation + scale(FireFreq^2) * Vegetation + scale(MAP) + scale(MAT) + scale(Duration) + (1 | Site), REML = F, data = responses_herbaceoustraits)
head(dredge(herbs_mod_height))

best_mod_herbheight <- lme(logRR_Height ~ scale(FireFreq) + scale(FireFreq^2), random = ~1|Site, data = responses_herbaceoustraits)
summary(best_mod_herbheight)
r.squaredGLMM(best_mod_herbheight)
best_mod_herbheight <- lme(logRR_Height ~ (FireFreq) + I(FireFreq^2), random = ~1|Site, data = responses_herbaceoustraits)
summary(best_mod_herbheight)

gg_heightherb <- visreg(lmer(logRR_Height ~ scale(FireFreq) + scale(FireFreq^2) + (1 | Site), data = responses_herbaceoustraits), "FireFreq", gg = T, line = list(col = "black"), band = F) +
  ggtitle("Theoretical maximum height") +
  geom_point(aes(color=responses_herbaceoustraits$Vegetation), size = I(3))
trait_figure(gg_heightherb)

# N stats
herbs_mod_N <- lmer(logRR_N ~ scale(FireFreq) * Vegetation + scale(FireFreq^2) * Vegetation + scale(MAP) + scale(MAT) + scale(Duration) + (1 | Site), REML = F, data = responses_herbaceoustraits)
head(dredge(herbs_mod_N))

best_mod_herbN <- lme(logRR_N ~ scale(MAP), random = ~1|Site, data = responses_herbaceoustraits)
summary(best_mod_herbN)
r.squaredGLMM(best_mod_herbN)
best_mod_herbN <- lme(logRR_N ~ MAP, random = ~1|Site, data = responses_herbaceoustraits)
summary(best_mod_herbN)

gg_Nherb <- visreg(lmer(logRR_N ~ scale(MAP) + (1 | Site), data = responses_herbaceoustraits), "MAP", gg = T, line = list(col = "black"), band = F) +
  ggtitle("Nitrogen per unit mass") +
  geom_point(aes(color=responses_herbaceoustraits$Vegetation), size = I(3))

gg_Nherb +
  ylab('log response ratio') +
  xlab('MAP (mm)') +
  theme_classic() + 
  theme(axis.ticks.length = unit(.3, 'cm'), 
        axis.title = element_text(size = 16), axis.text = element_text(size = 12, color = "black"), 
        axis.ticks = element_line(color = "black"),
        plot.title = element_text(hjust = 0.5),
        legend.position = "none") +
  scale_y_continuous(expand = c(0, 0), limits = c(-3,1), breaks = seq(-3,1,0.5)) +
  scale_x_continuous(expand = c(0, 0), limits = c(200,1880), breaks = seq(200, 1800, by = 400)) +
  geom_hline(yintercept=0, linetype='dashed') +
  scale_color_manual(values=c("green4", "gold")) +
  coord_cartesian(clip = "off")

# P stats
herbs_mod_P <- lmer(logRR_P ~ scale(FireFreq) * Vegetation + scale(FireFreq^2) * Vegetation + scale(MAP) + scale(MAT) + scale(Duration) + (1 | Site), REML = F, data = responses_herbaceoustraits)
head(dredge(herbs_mod_P))

best_mod_herbP <- lme(logRR_P ~ scale(MAP), random = ~1|Site, data = responses_herbaceoustraits)
summary(best_mod_herbP)
r.squaredGLMM(best_mod_herbP)
best_mod_herbP <- lme(logRR_P ~ MAP, random = ~1|Site, data = responses_herbaceoustraits)
summary(best_mod_herbP)

gg_Pherb <- visreg(lmer(logRR_P ~ scale(MAP) + (1 | Site), data = responses_herbaceoustraits), "MAP", gg = T, line = list(col = "black"), band = F) +
  ggtitle("Phosphorus per unit mass") +
  geom_point(aes(color=responses_herbaceoustraits$Vegetation), size = I(3))

gg_Pherb +
  ylab('log response ratio') +
  xlab('MAP (mm)') +
  theme_classic() + 
  theme(axis.ticks.length = unit(.3, 'cm'), 
        axis.title = element_text(size = 16), axis.text = element_text(size = 12, color = "black"), 
        axis.ticks = element_line(color = "black"),
        plot.title = element_text(hjust = 0.5),
        legend.position = "none") +
  scale_y_continuous(expand = c(0, 0), limits = c(-3,1), breaks = seq(-3,1,0.5)) +
  scale_x_continuous(expand = c(0, 0), limits = c(200,1880), breaks = seq(200, 1800, by = 400)) +
  geom_hline(yintercept=0, linetype='dashed') +
  scale_color_manual(values=c("green4", "gold")) +
  coord_cartesian(clip = "off")

# Stem stats
responses_herbaceoustraits2 <- responses_herbaceoustraits %>%   filter(!is.na(logRR_Stem))
herbs_mod_stem <- lmer(logRR_Stem ~ scale(FireFreq) * Vegetation + scale(FireFreq^2) * Vegetation + scale(MAP) + scale(MAT) + scale(Duration) + (1 | Site), REML = F, data = responses_herbaceoustraits2)
head(dredge(herbs_mod_stem))

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



## Woody traits ##

woody_traits <- read.csv("TraitListforWoodySpecies.csv")

seeds_woody <- stems %>% group_by(Site, Plot, FireFreq, FireStart, FireLast, FireSeason, SurveyYear) %>% 
  mutate(seed = woody_traits$DrySeedMass[match(GenusSpecies, woody_traits$oldspec)]) %>% 
  filter(!is.na(seed)) %>% 
  mutate(N = sum(Abundance), p = Abundance/N, seed_frac = p * seed) %>% 
  summarise(seed_cwm = mean(seed_frac))

height_woody <- stems %>% group_by(Site, Plot, FireFreq, FireStart, FireLast, FireSeason, SurveyYear) %>% 
  mutate(height = woody_traits$MaxHeight[match(GenusSpecies, woody_traits$oldspec)]) %>% 
  filter(!is.na(height)) %>% 
  mutate(N = sum(Abundance), p = Abundance/N, height_frac = p * height) %>% 
  summarise(height_cwm = mean(height_frac))

Nitrogen_woody <- stems %>% group_by(Site, Plot, FireFreq, FireStart, FireLast, FireSeason, SurveyYear) %>% 
  mutate(Nitrogen = woody_traits$LeafNitrogen[match(GenusSpecies, woody_traits$oldspec)]) %>% 
  filter(!is.na(Nitrogen)) %>% 
  mutate(N = sum(Abundance), p = Abundance/N, Nitrogen_frac = p * Nitrogen) %>% 
  summarise(Nitrogen_cwm = mean(Nitrogen_frac))

LMA_woody <- stems %>% group_by(Site, Plot, FireFreq, FireStart, FireLast, FireSeason, SurveyYear) %>% 
  mutate(LMA = woody_traits$LMA[match(GenusSpecies, woody_traits$oldspec)]) %>% 
  filter(!is.na(LMA)) %>% 
  mutate(N = sum(Abundance), p = Abundance/N, LMA_frac = p * LMA) %>% 
  summarise(LMA_cwm = mean(LMA_frac))

WoodDensity_woody <- stems %>% group_by(Site, Plot, FireFreq, FireStart, FireLast, FireSeason, SurveyYear) %>% 
  mutate(density = woody_traits$WoodDensity[match(GenusSpecies, woody_traits$oldspec)]) %>% 
  filter(!is.na(density)) %>% 
  mutate(N = sum(Abundance), p = Abundance/N, density_frac = p * density) %>% 
  summarise(density_cwm = mean(density_frac))

cwm_list_woody <- list(seeds_woody, height_woody, Nitrogen_woody, LMA_woody, WoodDensity_woody)
cwm_woody <- cwm_list_woody %>% reduce(left_join, by = c("Site" = "Site", "Plot" = "Plot", "FireFreq" = "FireFreq", "FireStart" = "FireStart", "FireLast" = "FireLast", "FireSeason" = "FireSeason", "SurveyYear" = "SurveyYear"))

cwm_woody <- cwm_woody  %>% arrange(Site, Plot, FireFreq, desc(SurveyYear)) %>% 
  distinct(Site, Plot, FireFreq, .keep_all = TRUE)

responses_woodytraits <- cwm_woody %>% filter(FireFreq != 0) %>% group_by(Site, FireFreq) %>% 
  summarise(meanSeed = mean(seed_cwm), meanHeight = mean(height_cwm), meanN = mean(Nitrogen_cwm), meanLMA = mean(LMA_cwm), meanDensity = mean(density_cwm),
            CV_Seed = (sd(seed_cwm))/meanSeed , CV_Height = (sd(height_cwm)) / meanHeight, CV_N = (sd(Nitrogen_cwm)) / meanN, CV_LMA = (sd(LMA_cwm))/meanLMA, CV_Density = (sd(density_cwm)) / meanDensity,
            n = n())
controls_woodytraits <- cwm_woody %>% filter(FireFreq == 0) %>% group_by(Site) %>% 
  summarise(meanSeed = mean(seed_cwm), meanHeight = mean(height_cwm), meanN = mean(Nitrogen_cwm), meanLMA = mean(LMA_cwm), meanDensity = mean(density_cwm),
            CV_Seed = (sd(seed_cwm))/meanSeed , CV_Height = (sd(height_cwm)) / meanHeight, CV_N = (sd(Nitrogen_cwm)) / meanN, CV_LMA = (sd(LMA_cwm))/meanLMA, CV_Density = (sd(density_cwm)) / meanDensity,
            n = n())

for (i in 1:nrow(responses_woodytraits)) {
  for (j in 1:nrow(controls_woodytraits)) {
    if (responses_woodytraits$Site[i] == controls_woodytraits$Site[j]) {
      
      responses_woodytraits$logRR_Seed1[i] <- log(responses_woodytraits$meanSeed[i]/controls_woodytraits$meanSeed[j])
      responses_woodytraits$logRR_Seed[i] <- log(responses_woodytraits$meanSeed[i]/controls_woodytraits$meanSeed[j]) + 0.5 * ( (((responses_woodytraits$CV_Seed[i])^2)/controls_woodytraits$n[j]) - (((controls_woodytraits$CV_Seed[j])^2)/responses_woodytraits$n[i]) )
      responses_woodytraits$logRR_v_Seed[i] <- ((responses_woodytraits$CV_Seed[i])^2)/responses_woodytraits$n[i] + ((controls_woodytraits$CV_Seed[j])^2)/controls_woodytraits$n[j] + ((responses_woodytraits$CV_Seed[i])^4)/(2*(responses_woodytraits$n[i])^2) + ((controls_woodytraits$CV_Seed[j])^4)/(2*(controls_woodytraits$n[j])^2)
      
      responses_woodytraits$logRR_Height1[i] <- log(responses_woodytraits$meanHeight[i]/controls_woodytraits$meanHeight[j])
      responses_woodytraits$logRR_Height[i] <- log(responses_woodytraits$meanHeight[i]/controls_woodytraits$meanHeight[j]) + 0.5 * ( (((responses_woodytraits$CV_Height[i])^2)/controls_woodytraits$n[j]) - (((controls_woodytraits$CV_Height[j])^2)/responses_woodytraits$n[i]) )
      responses_woodytraits$logRR_v_Height[i] <- ((responses_woodytraits$CV_Height[i])^2)/responses_woodytraits$n[i] + ((controls_woodytraits$CV_Height[j])^2)/controls_woodytraits$n[j] + ((responses_woodytraits$CV_Height[i])^4)/(2*(responses_woodytraits$n[i])^2) + ((controls_woodytraits$CV_Height[j])^4)/(2*(controls_woodytraits$n[j])^2)
      
      responses_woodytraits$logRR_N1[i] <- log(responses_woodytraits$meanN[i]/controls_woodytraits$meanN[j])
      responses_woodytraits$logRR_N[i] <- log(responses_woodytraits$meanN[i]/controls_woodytraits$meanN[j]) + 0.5 * ( (((responses_woodytraits$CV_N[i])^2)/controls_woodytraits$n[j]) - (((controls_woodytraits$CV_N[j])^2)/responses_woodytraits$n[i]) )
      responses_woodytraits$logRR_v_N[i] <- ((responses_woodytraits$CV_N[i])^2)/responses_woodytraits$n[i] + ((controls_woodytraits$CV_N[j])^2)/controls_woodytraits$n[j] + ((responses_woodytraits$CV_N[i])^4)/(2*(responses_woodytraits$n[i])^2) + ((controls_woodytraits$CV_N[j])^4)/(2*(controls_woodytraits$n[j])^2)
      
      responses_woodytraits$logRR_LMA1[i] <- log(responses_woodytraits$meanLMA[i]/controls_woodytraits$meanLMA[j])
      responses_woodytraits$logRR_LMA[i] <- log(responses_woodytraits$meanLMA[i]/controls_woodytraits$meanLMA[j]) + 0.5 * ( (((responses_woodytraits$CV_LMA[i])^2)/controls_woodytraits$n[j]) - (((controls_woodytraits$CV_LMA[j])^2)/responses_woodytraits$n[i]) )
      responses_woodytraits$logRR_v_LMA[i] <- ((responses_woodytraits$CV_LMA[i])^2)/responses_woodytraits$n[i] + ((controls_woodytraits$CV_LMA[j])^2)/controls_woodytraits$n[j] + ((responses_woodytraits$CV_LMA[i])^4)/(2*(responses_woodytraits$n[i])^2) + ((controls_woodytraits$CV_LMA[j])^4)/(2*(controls_woodytraits$n[j])^2)
      
      responses_woodytraits$logRR_Density1[i] <- log(responses_woodytraits$meanDensity[i]/controls_woodytraits$meanDensity[j])
      responses_woodytraits$logRR_Density[i] <- log(responses_woodytraits$meanDensity[i]/controls_woodytraits$meanDensity[j]) + 0.5 * ( (((responses_woodytraits$CV_Density[i])^2)/controls_woodytraits$n[j]) - (((controls_woodytraits$CV_Density[j])^2)/responses_woodytraits$n[i]) )
      responses_woodytraits$logRR_v_Density[i] <- ((responses_woodytraits$CV_Density[i])^2)/responses_woodytraits$n[i] + ((controls_woodytraits$CV_Density[j])^2)/controls_woodytraits$n[j] + ((responses_woodytraits$CV_Density[i])^4)/(2*(responses_woodytraits$n[i])^2) + ((controls_woodytraits$CV_Density[j])^4)/(2*(controls_woodytraits$n[j])^2)
      
    }
  }
}

responses_woodytraits <- responses_woodytraits %>% mutate(logRR_Seed = coalesce(logRR_Seed, logRR_Seed1),
                                                                    logRR_Height = coalesce(logRR_Height, logRR_Height1),
                                                                    logRR_N = coalesce(logRR_N, logRR_N1),
                                                                    logRR_LMA = coalesce(logRR_LMA, logRR_LMA1),
                                                                    logRR_Density = coalesce(logRR_Density, logRR_Density1)) %>% 
  match_sitedata()


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

### Statistics

options(na.action = "na.fail")
rcorr(as.matrix(responses_woodytraits[,c(30:33)]),type="pearson")

# Seed
woody_mod_seed <- lmer(logRR_Seed ~ scale(FireFreq) * Vegetation + scale(FireFreq^2) * Vegetation + scale(MAP) + scale(MAT) + (1 | Site), REML = F, data = responses_woodytraits)
head(dredge(woody_mod_seed))

best_mod_woodyseed <- lme(logRR_Seed ~ scale(FireFreq) + scale(FireFreq^2), random = ~1|Site, data = responses_woodytraits)
summary(best_mod_woodyseed)
r.squaredGLMM(best_mod_woodyseed)
best_mod_woodyseed <- lme(logRR_Seed ~ FireFreq + I(FireFreq^2), random = ~1|Site, data = responses_woodytraits)
summary(best_mod_woodyseed)

gg_woodyseed <- visreg(lmer(logRR_Seed ~ scale(FireFreq) + scale(FireFreq^2) + (1 | Site), data = responses_woodytraits), "FireFreq", gg = T, line = list(col = "black"), band = F) +
  ggtitle("Seed mass") +
  geom_point(aes(color=responses_woodytraits$Vegetation), size = I(3))
trait_figure_woody(gg_woodyseed)

# Height
woody_mod_height <- lmer(logRR_Height ~ scale(FireFreq) * Vegetation + scale(FireFreq^2) * Vegetation + scale(MAP) + scale(MAT) + (1 | Site), REML = F, data = responses_woodytraits)
head(dredge(woody_mod_height))
best_mod_woodyheight <- lme(logRR_Height ~ scale(FireFreq) + scale(FireFreq^2) + scale(MAT), random = ~1|Site, data = responses_woodytraits)
summary(best_mod_woodyheight)
r.squaredGLMM(best_mod_woodyheight)
best_mod_woodyheight <- lme(logRR_Height ~ FireFreq + I(FireFreq^2) + MAT, random = ~1|Site, data = responses_woodytraits)
summary(best_mod_woodyheight)

gg_woodyheight <- visreg(lmer(logRR_Height ~ scale(FireFreq) + scale(FireFreq^2) + scale(MAT) + (1 | Site), data = responses_woodytraits), "FireFreq", gg = T, line = list(col = "black"), band = F) +
  ggtitle("Maximum theoretical height") +
  geom_point(aes(color=responses_woodytraits$Vegetation), size = I(3))
trait_figure_woody(gg_woodyheight)

# N

responses_woodytraits2 <- responses_woodytraits %>% filter(!is.na(logRR_N))
woody_mod_N <- lmer(logRR_N ~ scale(FireFreq) * Vegetation + scale(FireFreq^2) * Vegetation + scale(MAP) + scale(MAT) + (1 | Site), REML = F, data = responses_woodytraits2)
head(dredge(woody_mod_N))

best_mod_woodyN <- lme(logRR_N ~ scale(FireFreq) + scale(FireFreq^2) + scale(MAT) + Vegetation, random = ~1|Site, data = responses_woodytraits2)
summary(best_mod_woodyN)
r.squaredGLMM(best_mod_woodyN)
best_mod_woodyN <- lme(logRR_N ~ FireFreq + I(FireFreq^2) + MAT + Vegetation, random = ~1|Site, data = responses_woodytraits2)
summary(best_mod_woodyN)

gg_woodyN <- visreg(lmer(logRR_N ~ scale(FireFreq) + scale(FireFreq^2) + scale(MAT) + (1 | Site), data = responses_woodytraits2), "FireFreq", gg = T, line = list(col = "black"), band = F) +
  ggtitle("Leaf Nitrogen") +
  geom_point(aes(color=responses_woodytraits2$Vegetation), size = I(3))
trait_figure_woody(gg_woodyN)

# LMA
responses_woodytraits3 <- responses_woodytraits %>% filter(!is.na(logRR_LMA))
# removed MAT bc singular
woody_mod_LMA <- lmer(logRR_LMA ~ scale(FireFreq) * Vegetation + scale(FireFreq^2) * Vegetation + scale(MAP) + (1 | Site), REML = F, data = responses_woodytraits3)
head(dredge(woody_mod_LMA))

best_mod_LMAwoody <- lme(logRR_LMA ~ scale(FireFreq) * Vegetation + scale(FireFreq^2), random = ~1|Site, data = responses_woodytraits3)
summary(best_mod_LMAwoody)
r.squaredGLMM(best_mod_LMAwoody)

best_mod_LMAwoody <- lme(logRR_LMA ~ FireFreq * Vegetation + I(FireFreq^2), random = ~1|Site, data = responses_woodytraits3)
summary(best_mod_LMAwoody)

gg_woodyLMA <- visreg(lmer(logRR_LMA ~ scale(FireFreq) * Vegetation + scale(FireFreq^2) + scale(MAT) + (1 | Site), data = responses_woodytraits3), "FireFreq", by = "Vegetation", overlay = T, gg = T, band = F) +
  ggtitle("LMA") +
  geom_point(aes(color=Vegetation), size = I(3))
trait_figure_woody(gg_woodyLMA)


# Wood Density

responses_woodytraits4 <- responses_woodytraits %>% filter(!is.na(logRR_Density))
woody_mod_density <- lmer(logRR_Density ~ scale(FireFreq) * Vegetation + scale(FireFreq^2) * Vegetation + scale(MAP) + scale(MAT) + (1 | Site), REML = F, data = responses_woodytraits4)
head(dredge(woody_mod_density))

best_mod_woodydensity <- lme(logRR_Density ~ scale(FireFreq) * Vegetation + scale(FireFreq^2) + scale(MAT), random = ~1|Site, data = responses_woodytraits4)
summary(best_mod_woodydensity)
r.squaredGLMM(best_mod_woodydensity)
best_mod_woodydensity <- lme(logRR_Density ~ FireFreq * Vegetation + I(FireFreq^2) + MAT, random = ~1|Site, data = responses_woodytraits4)
summary(best_mod_woodydensity)

gg_woodydensity <- visreg(lmer(logRR_Density ~ scale(FireFreq) * Vegetation + scale(FireFreq^2) + scale(MAT) + (1 | Site), data = responses_woodytraits4), "FireFreq", by = "Vegetation", overlay = T, gg = T, band = F) +
  ggtitle("Wood density") +
  geom_point(aes(color=Vegetation), size = I(3))
trait_figure_woody(gg_woodydensity)





### Rarefaction

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
  coord_cartesian(clip = "off")

ggplot(responses_herbaceous2, aes(x = cover_adj, y = logRR2, color = Vegetation)) + geom_point() +
  ylab ('log response ratio') + xlab("")
rarefaction <- lm(logRR2 ~ cover_adj*Vegetation, data = responses_herbaceous2)


  