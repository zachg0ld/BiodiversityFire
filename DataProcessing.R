setwd("~/Desktop/Research/FireBiodiversity/Data")

library(dplyr)
library(tidyr)

### PART I: Data Processing for biodiversity metrics ###


# Calculate plot-level richness (S), Shannon diversity (H), evenness (J)
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

# Calculate response ratios for each fire frequency at each site
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

# Match site and climate data
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
  
  return(BC_frame2)
}

herbaceous_dissimilarity <- dissimilarity(herbaceous)
grass_dissimilarity <- dissimilarity(grass)
nongrass_dissimilarity <- dissimilarity(nongrass)
woody_dissimilarity <- dissimilarity(stems)

# Write outputs to .csv files
write.csv(responses_herbs_all, file = "ProcessedData/herbaceous_responses_all.csv", row.names = F)
write.csv(responses_woody_all, file = "ProcessedData/woody_responses_all.csv", row.names = F)
write.csv(responses_grass_all, file = "ProcessedData/grass_responses_all.csv", row.names = F)
write.csv(responses_nongrass_all, file = "ProcessedData/nongrass_responses_all.csv", row.names = F)

write.csv(responses_herbaceous2, file = "ProcessedData/herbaceous_responses_plotlevel.csv", row.names = F)

write.csv(herbaceous_dissimilarity, file = "ProcessedData/herbaceous_dissimilarity.csv", row.names = F)
write.csv(grass_dissimilarity, file = "ProcessedData/grass_dissimilarity.csv", row.names = F)
write.csv(nongrass_dissimilarity, file = "ProcessedData/nongreass_dissimilarity.csv", row.names = F)
write.csv(woody_dissimilarity, file = "ProcessedData/woody_dissimilarity.csv", row.names = F)




### PART II: Traits ###

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


## Woody ##
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


write.csv(responses_herbaceoustraits, file = "ProcessedData/herbaceous_traits.csv", row.names = F)
write.csv(responses_woodytraits, file = "ProcessedData/woody_traits.csv", row.names = F)
