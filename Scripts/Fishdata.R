# Clean and average the Fish Data $$

# read libraries
library(tidyverse)
library(here)

# read in the fish data
fish<-read_csv(here("Data","raw_data","MCR_LTER_Annual_Fish_Survey_20260304.csv"))

# This calculates the total fish data
fish_clean<-fish %>%
  filter(Site == "LTER_1",
         Habitat == "Backreef") %>%
  filter(Biomass < 8000) %>% # there are 3 big sharks in the entire dataset that are biassing the biomass data.  Dropping them
  mutate(Biomass = ifelse(Biomass<0, NA, Biomass)) %>% # negative values are used for missing data
  #mutate(fish_g_m2 = (Biomass/50*Swath)) %>%
  group_by(Year)%>%
  summarise(total_fish_g_m2 = sum(Biomass, na.rm  =TRUE)/1200) # 300 m2 by 4 transects is 1200 m2 per year

fish_clean %>%  
ggplot(aes(x = Year, y = total_fish_g_m2))+
  geom_point()+
  geom_smooth(method = "lm", color= "black")+
  labs(y = "Mean total fish biomass (g/ m2)",
       x = "Year")+
  theme_minimal()

modfish<-lm(total_fish_g_m2~Year, fish_clean)
  #0.85 g per year
## group by trophic levels
fish_summary<-fish %>%
  filter(Site == "LTER_1",
         Habitat == "Backreef") %>%
  filter(Biomass < 8000) %>% # there are 3 big sharks in the entire dataset that are biassing the biomass data.  Dropping them
  mutate(Biomass = ifelse(Biomass<0, NA, Biomass)) %>% # negative values are used for missing data
  #mutate(fish_g_m2 = (Biomass/50*Swath)) %>%
  mutate(trophic_new = case_when(Fine_Trophic %in% c("Brusher", "Browser","Excavator","Concealed Cropper","Cropper", 
                                                     "Scraper","Herbivore/Detritivore" )~"Herbivores",
                                 Fine_Trophic  == "Corallivore"~"Corallivore",
                                 .default = "Other"
                                )
  )%>%
 # mutate(fish_kg_m2 = (Biomass/50*Swath)/1000)%>%
  group_by(Year, trophic_new) %>%
  summarise(fish_g_m2 = sum(Biomass, na.rm = TRUE)/1200)

fish_summary %>%
  left_join(fish_clean) %>%
  mutate(percent = fish_g_m2/total_fish_g_m2*100) %>%
  group_by(trophic_new)%>%
  summarise(mean_per = mean(percent, na.rm = TRUE))

fish_summary %>%
   ggplot(aes(x = Year, y = fish_g_m2))+
  geom_point()+
  geom_smooth(method = "lm")+
  # add_fishape(#family = "Labridae",
  #   option = "Chaetodon_plebeius",
  #   xmin = 2010, xmax = 2014, ymin = 0.15, ymax = 0.16,
  #   #fill = fish(option = "Chaetodon_plebeius"),
  #   alpha = 0.8) +
  facet_wrap(~trophic_new, scales = "free")

#write_csv(fish_clean, here("Data","fish_clean.csv"))

write_csv(fish_summary, here("Data","fish_summary.csv"))


fish2<-fish %>%
  filter(Site == "LTER_1",
         Habitat == "Backreef") %>%
  filter(Biomass < 8000) %>% # there are 3 big sharks in the entire dataset that are biassing the biomass data.  Dropping them
  mutate(Biomass = ifelse(Biomass<0, NA, Biomass)) %>% # negative values are used for missing data
  #mutate(fish_g_m2 = (Biomass/50*Swath)) %>%
  mutate(trophic_new = case_when(Fine_Trophic %in% c("Brusher", "Browser","Excavator","Concealed Cropper","Cropper", 
                                                     "Scraper","Herbivore/Detritivore" )~"Herbivores",
                                 Fine_Trophic  == "Corallivore"~"Corallivore",
                                 .default = "Other"
  )
  )

fish2 %>%
  separate(Taxonomy, " ", into = c("Genus", "Species")) %>%
  group_by(trophic_new) %>%
  summarise(
    families = paste(sort(unique(Fine_Trophic)), collapse = ", "),
    .groups = "drop"
  ) %>%
write_csv("fish_in_group2.csv")


unique(fish2$Taxonomy[fish2$trophic_new == "Herbivore"])
