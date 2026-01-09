# Calculate percent cover for Back Reef from Carpenter Data set ##

##### load libraries #############
library(tidyverse)
library(here)

# Bob's transect percent cover data
## Note these are 0.25 m2 each and there are 50 quads per site (12.5 m2)
BenthicCover_Algae<-read_csv(here("Data","raw_data","MCR_LTER_Annual_Survey_Benthic_Cover_20251009.csv")) %>%
  filter(Habitat == "Backreef")

Benthic_summary_Algae<-BenthicCover_Algae %>%
  rename(name = Taxonomy_Substrate_Functional_Group)%>%
  mutate(name = ifelse(name %in%c("Amansia rhodantha",         
                                  "Turbinaria ornata" ,        
                                  "Dictyota sp.",              
                                  "Halimeda sp.",              
                                  "Galaxaura sp.",             
                                  "Liagora ceranoides",        
                                  "Cyanophyta",                
                                  "Halimeda minima",           
                                  "Amphiroa fragilissima",     
                                  "Caulerpa serrulata",        
                                  "Corallimorpharia",          
                                  "Dictyota friabilis",        
                                  "Galaxaura rugosa",          
                                  "Cladophoropsis membranacea",
                                  "Galaxaura filamentosa",     
                                  "Halimeda discoidea",        
                                  "Peyssonnelia inamoena",     
                                  "Caulerpa racemosa",         
                                  "Valonia ventricosa",        
                                  "Actinotrichia fragilis",    
                                  "Dictyota bartayresiana",    
                                  "Microdictyon umbilicatum",  
                                  "Halimeda distorta",         
                                  "Halimeda incrassata",       
                                  "Halimeda macroloba",        
                                  "Dictyota implexa",          
                                  "Gelidiella acerosa",        
                                  "Dictyosphaeria cavernosa",  
                                  "Valonia aegagropila",       
                                  "Microdictyon okamurae",     
                                  "Halimeda opuntia",          
                                  "Dichotomaria obtusata",     
                                  "Chlorodesmis fastigiata",   
                                  "Phyllodictyon anastomosans",
                                  "Phormidium sp.",            
                                  "Cladophoropsis luxurians",  
                                  "Sargassum pacificum",       
                                  "Chnoospora implexa",        
                                  "Halimeda taenicola",        
                                  "Boodlea kaeneana",          
                                  "Padina boryana",            
                                  "Coelothrix irregularis",    
                                  "Gelidiella sp.",            
                                  "Hydroclathrus clathratus",  
                                  "Dictyota divaricata",       
                                  "Hypnea spinella",           
                                  "Dichotomaria marginata",    
                                  "Sporolithon sp.",           
                                  "Chaetomorpha antennina",    
                                  "Asparagopsis taxiformis"
  ), "Fleshy Macroalgae",name))%>%
  mutate(name = ifelse(name %in%c("Algal Turf", 
                                  "Damselfish Turf", 
                                  "Coral Rubble",
                                  "Lobophora variegata",
                                  "Shell Debris",  
                                  "Bare Space"
  ),"Turf/Cyanobacteria", name))%>%
  mutate(name = ifelse(name %in%c("Peyssonnelia bornetii", 
                                  "Peyssonnelia sp.", 
                                  "Sponge",
                                  "Tridacna sp.",
                                  "No data"
  ),"Sand", name))%>% # there is barely any of this in the dataset so group it with sand for the visual
  group_by(Year, Site, name)%>%
  summarise(total_cover = sum(Percent_Cover, na.rm = TRUE),
            mean_cover = 100*total_cover/5000) %>%
  mutate(Site = str_replace_all(Site, "_", " ")) # replace all the _ with spaces

write_csv(x = Benthic_summary_Algae, file = here("Data","Benthic_summary_algae.csv"))
