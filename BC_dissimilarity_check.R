## Hitchiti check (i'm just doing the same thing here but outside of the loop to make sure I get the same values)
cover_hitchiti <- cover %>% filter(Site == "Hitchiti")
df_hitchiti1 <- cover_hitchiti %>% 
  group_by(GenusSpecies) %>% 
  mutate(tmp_15 = max(0,Abundance[Plot == "15"]), tmp_21 = max(0,Abundance[Plot == "21"]), tmp_6 = max(0,Abundance[Plot == "6"])) %>% 
  group_by(Plot, GenusSpecies) %>% 
  mutate(min_15 = min(Abundance, tmp_15), min_21 = min(Abundance, tmp_21), min_6 = min(Abundance, tmp_6)) 

S_15 <- sum(df_hitchiti1$Abundance[df_hitchiti1$Plot == "15"])
S_21 <- sum(df_hitchiti1$Abundance[df_hitchiti1$Plot == "21"])
S_6 <- sum(df_hitchiti1$Abundance[df_hitchiti1$Plot == "6"])

df_hitchiti <- df_hitchiti1 %>% 
  group_by(Site, Plot, FireFreq) %>% 
  summarise(S = sum(Abundance), C_15 = sum(min_15), C_21 = sum(min_21), C_6 = sum(min_6)) %>% 
  mutate(BC_15 = 1 - ((2*C_15)/(S_15+S)), BC_21 = 1 - ((2*C_21)/(S_21+S)), BC_6 = 1 - ((2*C_6)/(S_6+S))) 

df_hitchiti <- df_hitchiti %>% select(c(Site, Plot, FireFreq, BC_2_0, BC_7_0)) %>% pivot_longer(cols = c(BC_2_0, BC_7_0))