### Extract Yearly Average and max daily temp values from the timeseries data ###

### load libraries ####
library(here)
library(tidyverse)

## Temperature Backreef LTER 1/2 
Temp_LTER1<-read_csv(here("Data","raw_data","MCR_LTER02_BTM_Backreef_Forereef_20251114.csv")) %>%
  filter(reef_type_code == "Backreef") 

# calculate the average daily temperature
Mean_daily_Temp<-Temp_LTER1 %>%
  mutate(date = as.Date(time_local)) %>%
  group_by(date)%>%
  summarise(daily_temp =mean(temperature_c, na.rm = TRUE))

# Get the max and mean daily temperature per year
Mean_year_Temp<-Mean_daily_Temp %>%
  mutate(Year = year(date)) %>%
  group_by(Year)%>%
  summarise(Mean_temp =mean(daily_temp  , na.rm = TRUE),
            Max_temp = max(daily_temp  , na.rm = TRUE))

# write this to use in analysis
write_csv(x = Mean_year_Temp, file = here("Data","InSituTemp.csv"))

