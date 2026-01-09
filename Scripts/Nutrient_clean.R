## Clean and Process the Turbinaria and water column N data data
library(here)
library(tidyverse)

# Turb %N
Turb<-read_csv(here("Data","raw_data","MCR_LTER_Macroalgal_CHN_2005_to_2024_20250616.csv")) %>%
  filter(Habitat == "Backreef",
         Genus == "Turbinaria",
         Site == "LTER_1") %>%
  drop_na(N)%>%
  group_by(Year) %>%
  summarise(N_percent = mean(N, na.rm = TRUE),
            C_percent = mean(C, na.rm = TRUE),
            CN = mean(CN_ratio, na.rm = TRUE)) %>%
  ungroup()

write_csv(Turb, here("Data","Turb_mean.csv"))

## Water nutrients - TAKE the mean from the bimonthly surveys (2005:2018)
water<-read_csv(here("Data","raw_data","WaterColumnN.csv"))%>% 
  mutate(Date = mdy(Date), Year = year(Date))%>%
  group_by(Year)%>%
  summarise_at(vars(Phosphate:Nitrite_and_Nitrate), "mean") %>%
  mutate(Site = "LTER 1")

write_csv(water, here("Data","Water_N_mean.csv"))

# Inpute the remaining water based on the turb data
AllNutrients<-Turb %>%
  full_join(water) %>%
  arrange(Year)

modN<-lm(Nitrite_and_Nitrate~N_percent, data = AllNutrients)
anova(modN)
summary(modN)

## predicted N+N from the turbinaria data which has a significant p value
AllNutrients$NN_impute<-predict(modN, newdata = AllNutrients)
