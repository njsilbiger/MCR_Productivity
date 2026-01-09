# Clean and average the Fish Data $$

# read libraries
library(tidyverse)
library(here)

# read in the fish data
fish<-read_csv(here("Data","raw_data","MCR_LTER_Annual_Fish_Survey_20250324.csv"))

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
  geom_smooth(method = "lm")+
  labs(y = "Mean total fish biomass (g/ m2)",
       x = "Year")+
  theme_bw()

## group by trophic levels
fish_summary<-fish %>%
  filter(Site == "LTER_1",
         Habitat == "Backreef") %>%
  filter(Biomass < 8000) %>% # there are 3 big sharks in the entire dataset that are biassing the biomass data.  Dropping them
  mutate(Biomass = ifelse(Biomass<0, NA, Biomass)) %>% # negative values are used for missing data
  #mutate(fish_g_m2 = (Biomass/50*Swath)) %>%
  mutate(trophic_new = case_when(Fine_Trophic %in% c("Brusher", "Browser","Excavator","Scraper")~"Herbivores",
                                 Fine_Trophic  == "Corallivore"~"Corallivore",
                                 .default = "Other"
                                )
  )%>%
 # mutate(fish_kg_m2 = (Biomass/50*Swath)/1000)%>%
  group_by(Year, trophic_new) %>%
  summarise(fish_g_m2 = sum(Biomass, na.rm = TRUE)/1200)


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
