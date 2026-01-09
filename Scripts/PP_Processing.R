# Clean the raw production data and calculate yearly averages
## We first run a Bayesian analysis to extract Pmax and ER for each year##

# load libraries
library(tidyverse)
library(here)
library(brms)

### Read in the data ######

filedir<-here("Data","raw_data","QC_PP") # file path for the QC PP files
files<-dir(path = filedir, pattern = ".csv", full.names = TRUE)

# bring in all the PP Data
All_PP_data<-files %>%
  set_names()%>% # set's the id of each list to the file name
  map_df(read_csv,.id = "filename") %>%
  mutate(DateTime = mdy_hm(DateTime))


# clean the data
# From 2007-2014 the PP units are in g O2/m2/h. 
#From June 2014 on, the units are mmol O2/m2/h. So those early rates will need to be converted
# 32 g/mol of O2

All_PP_data<-All_PP_data %>%
  mutate(UP_Oxy = ifelse(DateTime< ymd_hms("2014-04-01 00:00:00"), (UP_Oxy/32)*1000, UP_Oxy),
         DN_Oxy = ifelse(DateTime< ymd_hms("2014-04-01 00:00:00"), (DN_Oxy/32)*1000, DN_Oxy),
         PP = ifelse(DateTime< ymd_hms("2014-04-01 00:00:00"), (PP/32)*1000, PP))%>%
  mutate(Date = as_date(DateTime),
         DielDateTime = DateTime+hours(12), # all sampling started at noon such that midnight was the middle of a "day". Use this to extract the "daily" R to calculating GP
         DielDate = as_date(DielDateTime),
         Year = year(DateTime)) %>%
  filter(!Date %in% mdy("5/27/2011","5/28/2011","1/21/2014","05/25/2024" ))   # the respiration rate is incorrect these days from instrument failure

# find all the incomplete datasets (less than 24 hours) and only keep the complete ones
complete_dates<-All_PP_data %>% 
  group_by(DielDate)%>% 
  count() %>%
  filter(n == 24)

All_PP_data<-complete_dates %>%
  left_join(All_PP_data)


# calculate hourly GP and R data 
Daily_R <-All_PP_data %>%
  #filter(PAR==0)%>% # pull out all the night data
  group_by(Year, Season, DielDate) %>% # get the average nighttime respiration by day to add to NEP to calcualte GP
  summarise(R_average = mean(PP[PAR==0], na.rm = TRUE))

#Calculate GP and hourly in situ temperature and flow
All_PP_data<-All_PP_data %>%
  left_join(Daily_R) %>%
  mutate(GP = PP - R_average) %>%
  mutate(#GP = ifelse(PAR == 0, NA, GP),# remove GP from any of the night data 
    GP = ifelse(PAR == 0, NA, GP),
    GP = ifelse(GP<0, NA, GP),
    Temperature_mean = (UP_Temp+ DN_Temp)/2, # average temperature for the site
    Flow_mean = (UP_Velocity_mps+DN_Velocity_mps)/2) # average flow for the site

# Calculate yearly averages of Raw PP data
Year_Averages_PP<- All_PP_data %>%
  mutate(NP = PP,# remove nighttime respiration for average NP
         R = ifelse(PP<0, PP, NA) # only include night data for R
  ) %>% 
  group_by(Year) %>%
  summarise(NP_mean = mean(NP, na.rm = TRUE),
            NP_SE = sd(NP, na.rm = TRUE)/sqrt(n()),
            NP_max = max(NP, na.rm = TRUE),
            GP_mean = mean(GP, na.rm = TRUE),
            GP_SE = sd(GP, na.rm = TRUE)/sqrt(n()),
            R_mean = mean(R, na.rm = TRUE),
            R_SE = sd(R, na.rm = TRUE)/sqrt(n()),
            Temperature_mean = mean(Temperature_mean, na.rm = TRUE),
            Flow_mean = mean(Flow_mean, na.rm = TRUE),
            PAR_mean = mean(PAR[PAR>0], na.rm = TRUE)
  ) 

#### Bayesian analysis ###########
### Run PI curves allowing Rd and Pmax to change by year

# make a dataframe with names the brms likes
LTER1_Pnet <-All_PP_data %>%
  mutate(tempc = Temperature_mean,
         Flowmean = Flow_mean,
         flow_log = log(Flowmean), # log scale the flow data since it fits in a power function
         UPDN = as.factor(UPDN), # make factors for the analysis
         Year = as.factor(Year))

# Run a bayesian PI curve that varies by year -  We are doing this so that our "GP" values are light agnostic 
# This allows us to get good estimates without worrying about cloudy days
fit1_f<-brm(
  bf(PP ~ ((alpha*Pmax*PAR)/(alpha*PAR+Pmax))+Rd, 
     nl = TRUE, alpha~0+Year,Pmax~0+Year,
     Rd~0+Year)+ student(),
  data = LTER1_Pnet,
  set_rescor(FALSE),
  prior = c(
    prior(normal(0.1,10), nlpar = "alpha", lb = 0),
    prior(normal(-100,10), nlpar = "Rd", ub = 0), 
    prior(normal(100,10), nlpar = "Pmax", lb = 0) 
  ), 
  control = list(adapt_delta = 0.95, max_treedepth = 20), 
  cores = 3, 
  chains = 3, seed = 223, iter = 8000, warmup = 2000
  #silent = TRUE
) 

# extract the fixed effects
names<-rownames(fixef(fit1_f))

#extract the coefficients
coefficients<-as_tibble(fixef(fit1_f)) %>%
  mutate(params = names)

PRValues<-coefficients %>%
  separate(params, sep = "_",into = c("param_type","yearname"), remove = FALSE) %>%
  filter(param_type %in% c("Pmax","Rd")) %>% # extract Pmax and Rd
  separate(yearname, sep = "r", into = c("y","year")) %>% # clean up the years
  mutate(year = as.numeric(year))  %>% # make year a number again
  select(Year = year, Estimate, param_type) %>% # pivot wider so pmax and r in own column
  pivot_wider(names_from = param_type,
              values_from = Estimate)

# Bring it together with the raw averages 
Year_Averages_PP<-Year_Averages_PP %>%
  full_join(PRValues) %>%
  mutate(Rd = - Rd)  # make Rd positive

## write the yearly averages for PP data 
write_csv(Year_Averages_PP,here("Data","Year_Averages_PP.csv"))
