### Updated and cleaned analysis for primary production paper###

## Then we look at drivers of Pmax and ER with a focus on coral cover, temperature,
## and flow. Then we end with a Bayesian SEM that also includes NEC ###

## By Nyssa Silbiger ###
## Created on 2025-10-22 ####

##### load libraries #############
library(tidyverse)
library(here)
library(brms)
library(tidybayes)
library(lme4)
library(lmerTest)
library(scales)
library(patchwork)
library(ggtext)
library(viridis)
library(mgcv)
library(blavaan)
#library(semPlot)
library(bayesplot)
library(ggridges)
library(ggsci)
library(psych)
library(ggcorrplot)
library(ggcorrplot2)
library(ggeffects)


### Read in the data ######

## PP data 
Year_Averages_PP<-read_csv(here("Data","Year_Averages_PP.csv"))

# Turb %N
Turb<-read_csv(here("Data","Turb_mean.csv")) 

## Water nutrients 
water<-read_csv(here("Data","Water_N_mean.csv"))

#InSituData
insitu_temp<-read_csv(here("Data","InSituTemp.csv"))

# Total Fish Biomass - Trophic groups
#fish<-read_csv(here("Data","fish_clean.csv"))
fish_trophic<-read_csv(here("Data","fish_summary.csv")) %>%
  pivot_wider(names_from = trophic_new, 
              values_from = fish_g_m2) %>%
  select(!Other)

## read in the Benthic data
Benthic_summary_Algae<-read_csv(here("Data","Benthic_summary_algae.csv"))

############## ANALYSIS #######################
## Make a plot of the Benthic Data

myPal <- c(Coral = "#c38370", `Fleshy Macroalgae` = "#9da",
           `Turf/Cyanobacteria` =  "#9da993",
           #   millepora = "#bdc3cb", 
             `Crustose Corallines` = "#523a28", 
           Other = "#bdc3cb", 
           Sand = "#d6ad60")

# Make a plot
Benthic_summary_Algae %>%
  mutate(name = factor(name, levels = c("Coral","Fleshy Macroalgae",
                                        "Crustose Corallines","Turf/Cyanobacteria",
                                        "Sand","Other")))%>%
  filter(Site == "LTER 1")%>%
  ggplot(aes(x = Year, y = mean_cover, fill = name))+
  geom_bar(stat = "identity") +
  scale_fill_manual(values = myPal)+
  theme_classic() +
  labs(fill = "",
       y = "% Cover") +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.position = "bottom") 
  

ggsave(filename = here("Output","BenthicBob.pdf"), width = 8, height = 6)

# what is the coral cover in year 1 and year 2. Multiple that by 75.625 for N flux
Benthic_summary_Algae %>%
  filter(Year == 2006|Year == 2025) %>%
  filter(name == "Coral", Site == "LTER 1")

## Note: Linda's paper says Porites produce ~ 0.25 umol L-1 N+N per 8 hour incubation
# x 3 would give per 24 hours = ~0.75 umol L-1 - coral surface area was 124 cm2
# 0.75/124 gives umol L-1 Cm-2 x 1000 = m2 6.05 umol L-1 m-2 - multiply this by the
# percent coral cover and the surface area of the transect (12.5m2)
# 75.625 umol L-1 tranect of 100% Porites-1
#31 % coral cover to 6.8%
0.31*75.625 #23.44 umol L-1 transect-1
0.068*75.625 #5.14 umol L-1 transect-1 - 78% reduction in N production

## Calculate the total percent of calcifiers
Total_Calc<-Benthic_summary_Algae %>%
  filter(name %in% c("Coral","Crustose Corallines"))%>%
  group_by(Year, Site)%>%
  reframe(total_Calc = mean_cover[name == "Coral"]+
            mean_cover[name == "Crustose Corallines"])

# Total macroproducers
TotalLiving<-Benthic_summary_Algae %>%
  filter(name %in% c("Coral","Crustose Corallines","Fleshy Macroalgae"))%>%
  group_by(Year, Site)%>%
  summarise(mean_alive = sum(mean_cover),
            mean_fleshy = sum(mean_cover[name == "Fleshy Macroalgae"]),
            mean_coral = sum(mean_cover[name == "Coral"]))


#### yearly averages all together ####
Year_Averages <-  Year_Averages_PP %>% # production
  full_join(TotalLiving %>%
              filter(Site == "LTER 1")) %>% # Benthic Data
  full_join(water)%>% # water column N
  full_join(Turb)  %>% # add in the turbinaria data
  full_join(insitu_temp) %>% # read in annual temp data
  full_join(insitu_temp) %>% # bring in the in situ temperature data
  mutate(log_coral = log(mean_coral)) %>% # log transform the coral data
  left_join(fish_trophic) %>%
  filter(Year <2026) %>%
  arrange(Year)

### How is everything changing over time ####

# create a dataframe of standardized data
std_data<- Year_Averages %>%
  select(mean_coral, mean_fleshy, Pmax, Rd, 
          N_percent,
        Herbivores, Corallivore, NP_mean, Max_temp, GP_mean, Nitrite_and_Nitrate,
         Phosphate) %>%
  mutate(across(everything(), 
                ~as.numeric(scale(.x)))) %>%
  bind_cols(Year_Averages %>% select(Year))

## Everything below is a linear model showing standardized change per year

# Coral
coral_year<-brm(mean_coral~Year, data = std_data)
posterior_coral <- as_tibble(as.matrix(coral_year)) %>%
  select(Year = b_Year)%>%
  mutate(Parameter = "% Coral Cover")

# Macroalgae
fleshy_year<-brm(mean_fleshy~Year, data = std_data)
posterior_algae <- as_tibble(as.matrix(fleshy_year)) %>%
  select(Year = b_Year)%>%
  mutate(Parameter = "% Macroalgae Cover")

# Corallivore biomass
Corallivore_year<-brm(Corallivore~Year, data = std_data)
posterior_Corallivore <- as_tibble(as.matrix(Corallivore_year)) %>%
  select(Year = b_Year)%>%
  mutate(Parameter = "Corallivore Biomass")

# Herbivore Fish Biomass
dead_coral_fish_year<-brm(Herbivores~Year, data = std_data)
posterior_dead_coral_fish <- as_tibble(as.matrix(dead_coral_fish_year)) %>%
  select(Year = b_Year)%>%
  mutate(Parameter = "Herbivore/bioeroding fish biomass")

# In situ average daily max temperature
temp_year<-brm(Max_temp~Year, data = std_data)
posterior_temp <- as_tibble(as.matrix(temp_year)) %>%
  select(Year = b_Year)%>%
  mutate(Parameter = "Max Temperature")

# Percent N tubrinaria
N_year<-brm(N_percent~Year, data = std_data)
posterior_N <- as_tibble(as.matrix(N_year)) %>%
  select(Year = b_Year)%>%
  mutate(Parameter = "%N Content")

# water N + N
Nwater_year<-brm(Nitrite_and_Nitrate~Year, data = std_data)
posterior_Nwater <- as_tibble(as.matrix(Nwater_year)) %>%
  select(Year = b_Year)%>%
  mutate(Parameter = "Nitrate + Nitrite")

# Rd
Rd_year<-brm(Rd~Year, data = std_data)
posterior_Rd <- as_tibble(as.matrix(Rd_year)) %>%
  select(Year = b_Year)%>%
  mutate(Parameter = "Ecosystem Respiration")

# Pmax
Pmax_year<-brm(Pmax~Year, data = std_data)
posterior_Pmax <- as_tibble(as.matrix(Pmax_year)) %>%
  select(Year = b_Year)%>%
  mutate(Parameter = "Maximum Photosynthetic Capacity")

# NEP
NEP_year<-brm(NP_mean~Year, data = std_data)
posterior_NEP <- as_tibble(as.matrix(NEP_year)) %>%
  select(Year = b_Year)%>%
  mutate(Parameter = "Net Ecosystem Production")


# Bring together all the posterior data
All_posterior<-bind_rows(posterior_coral,
                         posterior_algae,
                         posterior_Corallivore,
                         posterior_dead_coral_fish,
                         posterior_N,
                         posterior_Rd,
                         posterior_Pmax,
                         posterior_NEP,
                         posterior_temp,
                         posterior_Nwater
                        )

# make a plot showing the change in each parameter over time
# Get the default NPG palette colors
npg_colors <- pal_npg("nrc")(10) # Default NPG palette has 10 colors

# Create a color ramp function to extend the palette to 11 colors
# This will interpolate between the existing NPG colors to generate an 11th color
npg_extended_palette_function <- colorRampPalette(npg_colors)

# Generate 10 colors from the extended palette
npg_11_colors <- npg_extended_palette_function(10)

All_posterior %>%
  ggplot(aes(x = Year, y = fct_reorder(Parameter, Year, mean), 
             fill = Parameter))+
  geom_vline(xintercept = 0)+
  geom_density_ridges(alpha = 0.5, color = NA)+
  scale_fill_manual(values = npg_11_colors)+
  annotate("text",x = 0.1, y = 0.75, label = "Increasing over time")+
  annotate("text",x = -0.1, y = 0.75, label = "Decreasing over time")+
  labs(x = "Standardized change per year",
       y = "")+
  lims(x = c(-0.3,0.3))+
  theme_ridges()+
  theme(legend.position = "none",
        axis.text.y = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        axis.text.x = element_text(size = 12))

ggsave(here("Output","ChangePosterior.pdf"), height = 6, width = 10)

### make a plot showing the bulk % change from start to end of timeseries

## tricky because the timeseries is different for each variable

get_first_non_na_per_column <- function(tbl) {
  # Apply a function to each column of the tibble
  # The function finds the index of the first non-NA value
  # and then extracts that value from the column.
  # If a column contains only NAs, NA is returned for that column.
  sapply(tbl, function(col) {
    first_non_na_index <- which(!is.na(col))[1]
    if (is.na(first_non_na_index)) {
      NA  # Return NA if no non-NA value is found
    } else {
      col[first_non_na_index]
    }
  })
}

get_last_non_na_per_column <- function(data_tibble) {
  # Ensure the input is a tibble
  if (!inherits(data_tibble, "tbl_df")) {
    stop("Input must be a tibble.")
  }
  
  # Apply a function to each column to find the last non-NA value
  last_values <- purrr::map(data_tibble, function(col) {
    non_na_values <- col[!is.na(col)]
    if (length(non_na_values) > 0) {
      return(tail(non_na_values, 1))
    } else {
      return(NA) # Return NA if all values in the column are NA
    }
  })
  
  # Convert the list to a named vector
  unlist(last_values)
}

# get the first value
first<-get_first_non_na_per_column(Year_Averages %>%
                              select(NP_mean, #GP_mean, 
                                     mean_fleshy,
                                     mean_coral, #NEC_mean_Day, 
                                     Pmax, Rd,
                                     N_percent, #C_percent, #mean_biomass,
                                     Herbivores, Corallivore, Max_temp,
                                     Nitrite_and_Nitrate #, Phosphate
                                     ))
first<-as_tibble(first) %>%
  mutate(Params = names(first)) %>%
  rename(first = value)


last<-get_last_non_na_per_column(Year_Averages %>%
                              select(NP_mean, #GP_mean, 
                                     mean_fleshy,
                                     mean_coral, #NEC_mean_Day, 
                                     Pmax, Rd,
                                     N_percent, #C_percent, #mean_biomass, 
                                     Herbivores, Corallivore, Max_temp,
                                     Nitrite_and_Nitrate, #Phosphate
                                     ))

last<-as_tibble(last) %>%
  mutate(Params = names(last))%>%
  rename(last = value)

# calculate the percent change in the variables
Per_change_var <- first %>%
  left_join(last) %>%
  mutate(percent_change = (last-first)/first*100) %>%
  mutate(percent_change = ifelse(Params == "NP_mean", -percent_change, percent_change)) %>%
  arrange(desc(percent_change))%>%
  mutate(nicenames = case_when(Params== "mean_fleshy" ~ "% Macroalgae Cover",
                               Params== "NP_mean" ~ "Net Ecosystem Production",
                               Params== "Max_temp" ~ "Max Temperature",
                               Params== "mean_biomass"~ "Fish Biomass",
                               Params== "Herbivores"~ "Hebivore/bioeroding fish",
                               Params == "Corallivore"~"Corallivore",
                               Params== "C_percent" ~ "% C Content",
                               Params== "Pmax" ~ "Maximum Photosynthetic Capacity",
                               Params== "GP_mean" ~ "Gross Ecosystem Production",
                               Params== "N_percent"~ "%N Content",
                               Params== "Rd" ~ "Ecosystem Respirataion",
                               Params== "NEC_mean_Day" ~ "Net Ecosystem Calcification",
                               Params== "mean_coral" ~ "% Coral Cover",
                               Params == "Nitrite_and_Nitrate" ~"Nitrate + Nitrite",
                               Params == "Phosphate" ~ "Phosphate")
  )


## Make a lillipop plot

## What are the number of years with data from each group
Num_years<-Year_Averages %>%
  select(NP_mean, GP_mean, mean_fleshy,
         mean_coral, #NEC_mean_Day,
         Pmax, Rd,
         N_percent, C_percent, #mean_biomass,
         Corallivore, Herbivores, Max_temp,
         Nitrite_and_Nitrate, Phosphate) %>%
  summarise_all(.funs = function(x){sum(!is.na(x))}) %>%
  pivot_longer(NP_mean:Phosphate) %>%
  rename(Params = name,
         N = value)

Per_change_var %>%
  left_join(Num_years) %>%
  mutate(nicenames = factor(nicenames, levels = c("Hebivore/bioeroding fish",
                                                  "% Macroalgae Cover",
                                                  "Fish Biomass",
                                                  "Net Ecosystem Production",
                                                  "Max Temperature",
                                                  "Gross Ecosystem Production",
                                                  "Phosphate",
                                                  "Maximum Photosynthetic Capacity",
                                                  "Net Ecosystem Calcification",
                                                  "Ecosystem Respirataion",
                                                  "Corallivore",
                                                  "% Coral Cover",
                                                  "%N Content",
                                                  "% C Content",
                                                  "Nitrate + Nitrite"
    
  ))) %>%
  mutate(color = ifelse(percent_change>0, "pos", "neg"))%>%
  ggplot(aes(x = percent_change, y = fct_rev(nicenames)))+
  geom_segment(aes(x = 0, xend = percent_change, yend = fct_rev(nicenames)), size = 1.1)+
  geom_point(aes(x = percent_change, color = color), size = 4)+
  geom_vline(xintercept = 0)+
  geom_text(aes(x = 180, y = nicenames, label = N), size = 6)+
  scale_color_manual(values = c("#B22222","#22B2B2"))+
  labs(x = "Percent change over collected time series (%)",
       y = "")+
  theme_minimal()+
  lims(x = c(-100,185))+
  theme(legend.position = "none",
        axis.text = element_text(size = 12),
        axis.title.x = element_text(size = 14))

ggsave(here("Output","lollipop.pdf"), height = 6, width = 8)  


ct <- corr.test(Year_Averages %>% 
                  mutate(log_fleshy = log(mean_fleshy)
                                             )%>%
                  select( NP = NP_mean, #GP = GP_mean,
                          Rd, Pmax,
                                          #NEC=NEC_mean_Day,
                                          `%N`=N_percent, 
                                          #`%C`=C_percent, 
                                          `N+N`=Nitrite_and_Nitrate, 
                                         # PO = Phosphate, 
                                          `% Coral`=mean_coral,
                                         `% Algae`= mean_fleshy, 
                                        # `Fish` = log_fish, 
                                          `Herbs/eroder` = Herbivores,
                                          `Corallivore` = Corallivore,
                                         `Max Temp`= Max_temp,
                                          #`Current Speed` = Flow_mean
                                        ), adjust = "none")
corr <- ct$r
p.mat <- ct$p

ggcorrplot.mixed(corr, 
                 upper = "ellipse", 
                 lower = "number", 
                 p.mat = p.mat, 
                insig = "label_sig", 
                sig.lvl = c(0.05, 0.01, 0.001))

ggsave(here("Output","correlations.pdf"), width = 10, height = 10)

## Run a Bayesian SEM to see how the different parameters are related
############### A BRMS example #############
# get standardized SEM data and rename for brns

# add in a variable to account for yearly disturbance that are above and beyond the effects of continuous temp change


yearresid<-resid(lm(Year ~ Max_temp, data = Year_Averages, na.action = na.pass))

semdata2<-Year_Averages %>%
  mutate( logfleshy = log(mean_fleshy), # log transform the data
          #logfish = log(mean_biomass),
          logrd = log(Rd),
          logpmax = log(Pmax),
          logproducers = log(mean_alive), # total producers
          N_percent = log(N_percent),
          yearresid = yearresid)%>%
  select(Year,meancoral = mean_coral, 
         logproducers,
         meanfleshy = mean_fleshy,
         temperature = Max_temp,
         flow = Flow_mean,
         Npercent = N_percent,
         #fish = mean_biomass,
         #GP = GP_mean,
         logrd,logpmax,
         logcoral = log_coral,
         logfleshy, #logfish, 
         Corallivore, 
         herbs = Herbivores,
         yearresid
        ) %>%
  drop_na(logcoral) %>%
  mutate(across(everything(), # scale the data
                ~as.numeric(scale(.x)))) 

# ----------------------------
# SEM with smooth terms + mi()
# ----------------------------

# imputation for missing data - 
# we hypothesize that temperature is changing by year and this allows us to impute the missing data by year
bf_temp <- bf(
  temperature | mi() ~ 1 + Year
)
# temp drives coral cover
bf_coraltemp <- bf( ### need to play with this because s makes a U and makes future predictions terrible
  logcoral | trunc(ub = 3.39)  ~ 1 + temperature + yearresid#+t2(Year)
)
## need to set the upper bound to be 100%
#(log(100)-modelmean$logcoral)/modelsd$logcoral
# 3.39

# LEVEL 2: Metabolic Rate Models
# GPP/Pmax is driven by the producers (corals, algae) and temperature
bf_logpmax <- bf(
  #logpmax | mi() ~ 1 +logproducers+ mi(temperature)
  logpmax | mi() ~ 1 +logcoral+ mi(temperature)
)

# ER is driven by corals as mixotrophs and abiotic factors 

bf_logrd <- bf(
  logrd | mi() ~ 1 +logcoral+ mi(temperature) #+mi(logpmax)
 # mvbind(logrd, logpmax)| mi() ~ 1 +logcoral+ mi(temperature) # rd and pmax are correlated errors
)
# Community models: Bottom up macroalgae drives herbivore biomass
bf_herbs <- bf(
  herbs | mi() ~ 1 +logfleshy
)

# Coral cover drives corallivore abundance
bf_corallivore <- bf(
  Corallivore | mi() ~ 1 +logcoral
)

# coral cover drives fleshy algae
bf_fleshy  <- bf(
  logfleshy | trunc(ub = 2.26) ~ logcoral)
# (log(100)-modelmean$logfleshy)/modelsd$logfleshy

# LEVEL 3: Ecosystem Function Models (Original Hypotheses)
# Nutrient recycling is driven by mean coral cover to Rd
bf_nutrients<- bf(
  Npercent | mi() ~ 1 +mi(logrd)
)


# COVARIANCES: Allow unexplained parts of external drivers to correlate
#Temperature_mean ~~ Flow_mean
#log_pmax ~~ log_rd
#log_coral~~log_fleshy

# algae and coral are correlated


# ----------------------------
# Auxiliary imputation models
# (intercept-only is fine to start; you can also add covariates/smooths here)
# ----------------------------
#bf_temp    <- bf(temperature | mi() ~ 1)
#bf_flow    <- bf(flow        | mi() ~ 1)
#bf_coral   <- bf(meancoral       | mi() ~ 1)
#bf_fleshy  <- bf(logfleshy | mi() ~ 1) # to allow for correlated errors between coral and macroalgae since algae is not a response variable
#bf_biomass <- bf(fish     | mi() ~ 1)

# ----------------------------
# Fit
# ----------------------------
brms_sem_full <- brm(
  bf_temp+bf_coraltemp+ bf_logrd + bf_herbs +bf_corallivore+bf_nutrients + bf_fleshy + bf_logpmax +
  set_rescor(FALSE),             # residual correlations among responses
  data = semdata2,
  chains = 3,
  iter = 10000,
  warmup = 5000,
  seed = 11,
  sample_prior = "no",
  prior = c(set_prior("normal(0, 2)", class = "b"),
            #set_prior("normal(0, 2)", class = "b", resp = "logcoral", ub = 3.39), # set upper bound to 100%
            set_prior("normal(0, 0.25)", resp = "temperature", class = "sigma", lb = 0),
            set_prior("normal(0, 2)", class = "Intercept")),
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  family = gaussian()
)

summary(brms_sem_full) # Rhat and ESS good for everything
# look at the trace plots
#plot(brms_sem_full) # looks good all chains well mixed
# look at PP checks (looks good)
pp_check(brms_sem_full, resp = "temperature", ndraws = 100)
pp_check(brms_sem_full, resp = "herbs", ndraws = 100)
pp_check(brms_sem_full, resp = "logcoral", ndraws = 100)
pp_check(brms_sem_full, resp = "Npercent", ndraws = 100)
pp_check(brms_sem_full, resp = "logpmax", ndraws = 100)
pp_check(brms_sem_full, resp = "logrd", ndraws = 100)
pp_check(brms_sem_full, resp = "logfleshy", ndraws = 100)

# plot all the posteriors
brms_sem_full %>%
  spread_draws(`b.*`, regex = TRUE) %>%
  select(b_temperature_Year, #bsp_logrd_milogpmax,#b_logpmax_logproducers, 
         b_logpmax_logcoral,
         b_logrd_logcoral,b_herbs_logfleshy,
         b_logfleshy_logcoral,
         b_Corallivore_logcoral, b_logcoral_temperature, bsp_logpmax_mitemperature,
         bsp_logrd_mitemperature, b_logcoral_yearresid,bsp_Npercent_milogrd) %>%
  pivot_longer(cols = b_temperature_Year:bsp_Npercent_milogrd)%>%
  ggplot(aes(x = value, y = name))+
  stat_halfeye(point_interval=median_hdi, .width=c(.95, .75), 
               fatten_point = 2, slab_alpha = 0.6, fill = scales::alpha("#009E73",0.6)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
#  facet_wrap(~name, ncol = 1)+
  theme_bw() + 
  xlab("")+
  ylab("")+
  theme(panel.grid.major = element_blank(), # remove gridlines
        panel.grid.minor = element_blank(), #remove gridlines
        strip.background = element_blank(), 
        legend.position = "none",  
        rect = element_rect(fill = "transparent"),  
        plot.background = element_rect(fill = "transparent", color = NA),
      #  axis.text.y = element_blank(),
       # axis.ticks.y = element_blank()
        )

# get the nonscaled values for the model
Yearmodelvalues<-Year_Averages %>%
  mutate( logfleshy = log(mean_fleshy),
          logproducers = log(mean_alive), # total producers
          #logfish = log(mean_biomass),
          logrd = log(Rd),
          logpmax = log(Pmax),
         # Rd =-R_mean,
         # NEC = NEC_mean_Day+NEC_mean_Night
         N_percent = log(N_percent),
         yearresid = yearresid
         )%>%
  select(Year,meancoral = mean_coral, 
         meanfleshy = mean_fleshy,
         temperature = Max_temp,
         flow = Flow_mean,
         Npercent = N_percent,
         #fish = mean_biomass,
         #GP = GP_mean,
         logrd,logpmax,
         logcoral = log_coral,logproducers,
         logfleshy, #logfish, 
         Corallivore, herbs = Herbivores,
         yearresid
  ) %>%
  drop_na(logcoral)

# get the means and SD to unscale the plots
modelmean<-Yearmodelvalues%>%
  mutate(across(everything(), # scale the data
                ~as.numeric(mean(.x, na.rm = TRUE)))) %>%
  distinct()

modelsd<-Yearmodelvalues%>%
  mutate(across(everything(), # scale the data
                ~as.numeric(sd(.x, na.rm = TRUE)))) %>%
  distinct()

# set the theme for everything
theme_regression<-
  theme_bw()+
  theme(text = element_text(size = 14, family = "Arial"),
        rect = element_rect(fill = "transparent"),  
        plot.background = element_rect(fill = "transparent", color = NA)
  )

theme_set(theme_regression) 

coral_temp<-conditional_effects(brms_sem_full, resp = "logcoral", effects = "temperature") 
ct_plot<-as_tibble(coral_temp$logcoral.logcoral_temperature) %>%
  mutate(temperature = temperature*modelsd$temperature+modelmean$temperature,
         estimate__ = exp(estimate__*modelsd$logcoral+modelmean$logcoral),
         lower__ = exp(lower__*modelsd$logcoral+modelmean$logcoral),
         upper__ = exp(upper__*modelsd$logcoral+modelmean$logcoral))  %>% # unscale the data
  ggplot(aes(x = temperature, y = estimate__))+
  geom_line(linewidth = 1)+
  geom_ribbon(aes(ymin =lower__, ymax =  upper__), alpha = 0.25, fill = "#67bed9")+
  geom_point(data = Year_Averages, aes(x = Max_temp, y = exp(log_coral)), alpha = 0.5)+
  coord_transform(y = "log")+
  scale_y_continuous(breaks = c(1,5,10,20,40))+
  labs(x = expression("Temperature ("*degree*"C)"),
       y = expression(atop("Coral Cover", "(%)")))
  
  ### Temperature Year
temp_year<-conditional_effects(brms_sem_full, resp = "temperature", effects = "Year") 
ty_plot<-as_tibble(temp_year$temperature.temperature_Year) %>%
  mutate(Year = Year*modelsd$Year+modelmean$Year,
         estimate__ = estimate__*modelsd$temperature+modelmean$temperature,
         lower__ = lower__*modelsd$temperature+modelmean$temperature,
         upper__ = upper__*modelsd$temperature+modelmean$temperature)  %>% # unscale the data
  ggplot(aes(x = Year, y = estimate__))+
  geom_line(linewidth = 1)+
  geom_ribbon(aes(ymin =lower__, ymax =  upper__), alpha = 0.25, fill = "#67bed9")+
  geom_point(data = Year_Averages, aes(x = Year, y = Max_temp), alpha = 0.5)+
  #coord_trans(y = "log")+
 # scale_y_continuous(breaks = c(1,5,10,20,50,100))+
  labs(y = expression(atop("Temperature", "("*degree*"C)")),
       x = "Year")

# log Pmax ~ total producers and temperature
pmax_coral<-conditional_effects(brms_sem_full, resp = "logpmax", effects = "logcoral") 
pc_plot<-as_tibble(pmax_coral$logpmax.logpmax_logcoral) %>%
  mutate(logcoral = exp(logcoral*modelsd$logcoral+modelmean$logcoral),
         estimate__ = exp(estimate__*modelsd$logpmax+modelmean$logpmax),
         lower__ = exp(lower__*modelsd$logpmax+modelmean$logpmax),
         upper__ = exp(upper__*modelsd$logpmax+modelmean$logpmax) ) %>% # unscale the data
  ggplot(aes(x = logcoral, y = estimate__))+
  geom_line(linewidth = 1)+
  geom_ribbon(aes(ymin =lower__, ymax =  upper__), alpha = 0.25, fill = "#67bed9")+
  geom_point(data = Year_Averages, aes(x = mean_coral, y = Pmax), alpha = 0.5)+
  coord_trans(x = "log", y = "log")+
  # scale_y_continuous(breaks = c(1,5,10,20,50,100))+
  labs(y = expression(atop("Pmax", "(mmol O "[2]*" m"^2*" hr"^-1*")")),
       x = "Total Producer Cover (%)")

# log Rd ~ coral and temperature
rd_coral<-conditional_effects(brms_sem_full, resp = "logrd", effects = "logcoral") 
rc_plot<-as_tibble(rd_coral$logrd.logrd_logcoral) %>%
  mutate(logcoral = exp(logcoral*modelsd$logcoral+modelmean$logcoral),
         estimate__ = exp(estimate__*modelsd$logrd+modelmean$logrd),
         lower__ = exp(lower__*modelsd$logrd+modelmean$logrd),
         upper__ = exp(upper__*modelsd$logrd+modelmean$logrd) ) %>% # unscale the data
  ggplot(aes(x = logcoral, y = estimate__))+
  geom_line(linewidth = 1)+
  geom_ribbon(aes(ymin =lower__, ymax =  upper__), alpha = 0.25, fill = "#67bed9")+
  geom_point(data = Year_Averages, aes(x = exp(log_coral), y = Rd), alpha = 0.5)+
  coord_trans(x = "log", y = "log")+
  # scale_y_continuous(breaks = c(1,5,10,20,50,100))+
  labs(y = expression(atop("Ecosystem Respiration", "(mmol O "[2]*" m"^2*" hr"^-1*")")),
       x = "Coral (%)")

# herbs ~ log fleshy
herbs_fleshy<-conditional_effects(brms_sem_full, resp = "herbs", effects = "logfleshy") 
hf_plot<-as_tibble(herbs_fleshy$herbs.herbs_logfleshy) %>%
  mutate(logfleshy = exp(logfleshy*modelsd$logfleshy+modelmean$logfleshy),
         estimate__ = estimate__*modelsd$herbs+modelmean$herbs,
         lower__ = lower__*modelsd$herbs+modelmean$herbs,
         upper__ = upper__*modelsd$herbs+modelmean$herbs ) %>% # unscale the data
  ggplot(aes(x = logfleshy, y = estimate__))+
  geom_line(linewidth = 1)+
  geom_ribbon(aes(ymin =lower__, ymax =  upper__), alpha = 0.25, fill = "#67bed9")+
  geom_point(data = Year_Averages, aes(x = mean_fleshy, y = Herbivores), alpha = 0.5)+
  coord_trans(x = "log")+
  # scale_y_continuous(breaks = c(1,5,10,20,50,100))+
  labs(y = expression(atop("Herbivore Fish Biomass", "(g m "^-1*")")),
       x = "Fleshy Macroalgae (%)")

# corallivore ~ corals
corallivore_coral<-conditional_effects(brms_sem_full, resp = "Corallivore", effects = "logcoral") 
cc_plot<-as_tibble(corallivore_coral$Corallivore.Corallivore_logcoral) %>%
  mutate(logcoral = exp(logcoral *modelsd$logcoral +modelmean$logcoral ),
         estimate__ = estimate__*modelsd$Corallivore+modelmean$Corallivore,
         lower__ = lower__*modelsd$Corallivore+modelmean$Corallivore,
         upper__ = upper__*modelsd$Corallivore+modelmean$Corallivore ) %>% # unscale the data
  ggplot(aes(x = logcoral , y = estimate__))+
  geom_line(linewidth = 1)+
  geom_ribbon(aes(ymin =lower__, ymax =  upper__), alpha = 0.25, fill = "#67bed9")+
  geom_point(data = Year_Averages, aes(x = exp(log_coral), y = Corallivore), alpha = 0.5)+
  coord_trans(x = "log")+
  # scale_y_continuous(breaks = c(1,5,10,20,50,100))+
  labs(y = expression(atop("Corallivore Fish Biomass", "(g m "^-1*")")),
       x = "Coral Cover (%)")

# Fleshy algae ~ corals
fleshy_coral<-conditional_effects(brms_sem_full, resp = "logfleshy", effects = "logcoral") 
fc_plot<-as_tibble(fleshy_coral$logfleshy.logfleshy_logcoral) %>%
  mutate(logcoral = exp(logcoral *modelsd$logcoral +modelmean$logcoral ),
         estimate__ = exp(estimate__*modelsd$logfleshy+modelmean$logfleshy),
         lower__ = exp(lower__*modelsd$logfleshy+modelmean$logfleshy),
         upper__ = exp(upper__*modelsd$logfleshy+modelmean$logfleshy )) %>% # unscale the data
  ggplot(aes(x = logcoral , y = estimate__))+
  geom_line(linewidth = 1)+
  geom_ribbon(aes(ymin =lower__, ymax =  upper__), alpha = 0.25, fill = "#67bed9")+
  geom_point(data = Year_Averages, aes(x = exp(log_coral), y = mean_fleshy), alpha = 0.5)+
  coord_trans(x = "log",
              y = "log")+
  # scale_y_continuous(breaks = c(1,5,10,20,50,100))+
  labs(y = expression(atop("Fleshy Macroalgae Cover", "(%)")),
       x = "Coral Cover (%)")

# N% ~ respiration
N_Rd<-conditional_effects(brms_sem_full, resp = "Npercent", effects = "logrd") 
nr_plot<-as_tibble(N_Rd$Npercent.Npercent_logrd) %>%
  mutate(logrd = exp(logrd*modelsd$logrd+modelmean$logrd),
         estimate__ = exp(estimate__*modelsd$Npercent+modelmean$Npercent),
         lower__ = exp(lower__*modelsd$Npercent+modelmean$Npercent),
         upper__ = exp(upper__*modelsd$Npercent+modelmean$Npercent )) %>% # unscale the data
  ggplot(aes(x = logrd, y = estimate__))+
  geom_line(linewidth = 1)+
  geom_ribbon(aes(ymin =lower__, ymax =  upper__), alpha = 0.25, fill = "#67bed9")+
  geom_point(data = Year_Averages, aes(x = Rd, y = N_percent), alpha = 0.5)+
  coord_trans(x = "log", y = "log")+
  # scale_y_continuous(breaks = c(1,5,10,20,50,100))+
  labs(x = expression("Ecosystem Respiration (mmol O"[2]*" m"^2*" hr"^-1*")"),
       y = expression(atop("Nitrogen Content", "(%)")))


(((plot_spacer()|nr_plot|plot_spacer())+plot_layout(widths = c(1,2,1)))/((plot_spacer()|rc_plot|plot_spacer())+plot_layout(widths = c(1,2,1)))/(cc_plot|hf_plot)/(ct_plot|fc_plot)/(((plot_spacer()|ty_plot|plot_spacer()))+
    plot_layout(widths = c(1,2,1))))+plot_annotation(tag_levels = "a")

ggsave(here("Output","RegressionSEMplots.pdf"), height = 12, width = 8, device = cairo_pdf)

######### Standardarized effect size plot 
gp_edges<-as_tibble(fixef(brms_sem_full)) %>%
  mutate(coef =rownames(fixef(brms_sem_full))) %>%
  separate(col = coef, into = c("resp","pred"))%>%
  filter(pred != "Intercept") %>%
  rename(est = Estimate, lo = `Q2.5`, hi = `Q97.5`)%>%
  mutate(sig = ifelse(sign(lo)==sign(hi),1,0.5))  %>% # determine significance
  mutate(resp = case_when( 
    resp == "temperature" ~"Temperature",
    resp == "logfleshy"~ "Fleshy Macroalgae",
    resp == "logpmax"~"Ecosystem Pmax",
    resp == "logrd"~"Ecosystem Rd",
    resp == "herbs" ~ "Herbivorous Fish" ,
    resp == "Corallivore" ~ "Corallivore Fish" ,
    resp == "logcoral"~ "Coral",
    resp == "Npercent" ~"N stock")) %>% 
  mutate(resp = factor(resp, c("Temperature","Coral","Fleshy Macroalgae","Herbivorous Fish",
                               "Corallivore Fish", "Ecosystem Rd",
                               "Ecosystem Pmax","N stock"))) %>%
  mutate(pred = case_when( 
    pred == "temperature" ~"Temperature",
    pred == "mitemperature" ~"Temperature",
    pred == "logproducers"~"Total Producer Cover",
    pred == "milogrd"~"Ecosystem Rd Rate",
    pred == "logcoral"~ "Coral Cover",
    pred == "logfleshy"~ "Fleshy Macroalgae Cover",
    pred == "Year" ~"Year",
    pred == "yearresid"~"Non-Temperature Yearly Variability")) %>%
  mutate(pred = factor(pred, c("Ecosystem Rd Rate","Total Producer Cover",
                               "Fleshy Macroalgae Cover","Coral Cover","Non-Temperature Yearly Variability",
                               "Temperature","Year")))

gp_edges %>%
  ggplot(aes(x = est, y = pred))+
  geom_point(size = 3, color = "firebrick", aes(alpha = sig))+
  geom_errorbarh(aes(xmin = lo, xmax = hi,alpha = sig), height = 0, color = "firebrick")+
  geom_vline(xintercept = 0)+
  scale_alpha(range = c(0.25,1))+
  labs(x = "Standardized Effect Size",
       y = "Predictor Variables")+
  theme_minimal()+
  facet_wrap(~resp, nrow = 1)+
  theme(text = element_text(size = 14),
        axis.text = element_text(size = 14),
        strip.text = element_text(size = 14),
        legend.position = "none")

ggsave(here("Output","SEMCoeffs.pdf"), height = 6, width = 18, device = cairo_pdf)

## get the correlated errors
#a<-summary(brms_sem_full) 

#cor_res<-as_tibble(a$rescor_pars) %>%
#  mutate(coef = rownames(a$rescor_pars)) %>%
#  mutate(coef = str_remove_all(coef, "[rescor()]")) %>%
#  separate(coef, into = c("x","y"), sep = ",")

#cor_res %>%
#  ggplot(aes(x=x, y=y, fill = Estimate))+
#  geom_tile()+
#  scale_fill_gradient(low = "white", high = "pink")+
#  theme_minimal()


#lgpmax,lgd
#lgal,lgflhy

# Make Predictions for different years
test_years<-c(2010,2020,2030,2040,2050)
test_years_scale<-(test_years-modelmean$Year)/modelsd$Year

# hindcast values 
hind_values<-Year_Averages %>%
  mutate(yearresid = yearresid)%>%
  filter(Year %in% c(2010,2020))

# Predict the temperature in past and future years/hindcast and forcast
new_data <- data.frame(Year = test_years_scale)
pred_new <- posterior_predict(brms_sem_full, newdata = new_data, resp = "temperature")

#convert back to temperature
Temp_pred_plot<-as_tibble(pred_new*modelsd$temperature+modelmean$temperature) %>%
  pivot_longer(cols = V1:V5) %>%
  mutate(year = case_when(name == "V1" ~"2010",
                          name == "V2" ~"2020",
                          name == "V3" ~"2030",
                          name == "V4" ~"2040",
                          name == "V5" ~"2050")) %>%
  ggplot(aes(y = value, x = year))+
  stat_halfeye(point_interval=mean_hdi, .width=c(.95, .75), 
               fatten_point = 2, slab_alpha = 0.6, fill = scales::alpha("#009E73",0.6))+
  geom_point(data = hind_values, aes(x = as.character(Year), y = Max_temp), 
             size = 4, color = "firebrick", position = position_dodge2(width = 1))+
  labs(x = "Year",
       y = "Temperature "~degree~"C")+
  lims(y = c(28,33))

## Take the predicted temperature from each year and calculate predicted coral cover
#newdata_temp<-data.frame(temperature = colMeans(pred_new), Year = test_years_scale) 
# get the expected "yearresidual" from future years
modyear<-lm(yearresid~Year, data = semdata2)
test_years_resid<-predict(modyear, newdata = data.frame(Year = test_years_scale))

newdata_temp<-data.frame(temperature = colMeans(pred_new), yearresid = test_years_resid) 
pred_new_coral <- posterior_predict(brms_sem_full, newdata = newdata_temp, resp = "logcoral")

Coral_pred_plot<-exp(as_tibble(pred_new_coral*modelsd$logcoral+modelmean$logcoral)) %>%
  pivot_longer(cols = V1:V5) %>%
  mutate(year = case_when(name == "V1" ~"2010",
                          name == "V2" ~"2020",
                          name == "V3" ~"2030",
                          name == "V4" ~"2040",
                          name == "V5" ~"2050")) %>%
  ggplot(aes(y = value, x = year))+
  stat_halfeye(point_interval=mean_hdi, .width=c(.95, .75), 
               fatten_point = 2, slab_alpha = 0.6, fill = scales::alpha("#009E73",0.6))+
  geom_point(data = hind_values, aes(x = as.character(Year), y = mean_coral), 
             size = 4, color = "firebrick", position = position_dodge2(width = 1))+
  labs(x = "Year",
       y = "Coral Cover (%)")+
  lims(y = c(0,25))

# predict Rd
newdata_logcoral<-data.frame(logcoral = colMeans(pred_new_coral),
                             temperature = newdata_temp$temperature) 

pred_new_Rd<- posterior_predict(brms_sem_full, newdata = newdata_logcoral, resp = "logrd")
Rd_pred_plot<-exp(as_tibble(pred_new_Rd*modelsd$logrd+modelmean$logrd)) %>%
  pivot_longer(cols = V1:V5) %>%
  mutate(year = case_when(name == "V1" ~"2010",
                          name == "V2" ~"2020",
                          name == "V3" ~"2030",
                          name == "V4" ~"2040",
                          name == "V5" ~"2050")) %>%
  ggplot(aes(y = value, x = year))+
  stat_halfeye(point_interval=mean_hdi, .width=c(.95, .75), 
               fatten_point = 2, slab_alpha = 0.6, fill = scales::alpha("#009E73",0.6))+
  geom_point(data = hind_values, aes(x = as.character(Year), y = Rd), 
             size = 4, color = "firebrick", position = position_dodge2(width = 1))+
  labs(x = "Year",
       y = expression(atop("Ecosystem Respiration", "(mmol O "[2]*" m"^2*" hr"^-1*")")))+
  lims(y = c(0,80))

# Predict N -
newdata_N<-data.frame(logrd = colMeans(pred_new_Rd)) 
pred_new_N<- posterior_epred(brms_sem_full, newdata = newdata_N, resp = "Npercent")
N_pred_plot<-exp(as_tibble(pred_new_N*modelsd$Npercent+modelmean$Npercent)) %>%
  pivot_longer(cols = V1:V5) %>%
  mutate(year = case_when(name == "V1" ~"2010",
                          name == "V2" ~"2020",
                          name == "V3" ~"2030",
                          name == "V4" ~"2040",
                          name == "V5" ~"2050")) %>%
  ggplot(aes(y = value, x = year))+
  stat_halfeye(point_interval=mean_hdi, .width=c(.95, .75), 
               fatten_point = 2, slab_alpha = 0.6, fill = scales::alpha("#009E73",0.6))+
  geom_point(data = hind_values, aes(x = as.character(Year), y = N_percent), 
             size = 4, color = "firebrick", position = position_dodge2(width = 1))+
  labs(x = "Year",
       y = "N (%)")

  N_pred_plot/Rd_pred_plot/Coral_pred_plot/Temp_pred_plot
  ggsave(here("Output","Predictions.pdf"), height = 10, width = 5,device = cairo_pdf)
  