#' ---
#' title: "Estimating drought tolerance"
#' author: "Julie Larson"
#' date: "17 August 2025"
#' output: github_document
#' ---


#' **Welcome!**
#'
#' This R script contains code to conduct analyses and create figures included 
#' in the following manuscript in preparation:
#' 
#' Larson, J.E., B.J. Butterfield, S.M. Munson, D. Neuhaus, S.M. Copeland. (In Prep) 
#' Functional traits can explain coordinated seedling responses to drought and freezing stress


#' **This code: Estimaing drought tolerance as drydown survival duration**
#' 
#' In this code, we use survival data from the drydown experiment to estimate 
#' 'drought tolerance' as the mean duration of seedling survival for each species.
#' We also characterize soil water dynamics during the experiment.



########################
#'
#' ** Load packages, options & functions **
#'
#'#####################

library(tidyverse)
library(GGally)


########################
#'
#' **Load or Source Data**
#'
#'#####################

#' Specify working directory
path <- "C:/Users/larsonju/OneDrive - UW/Documents/Seed_Team/Data/Data_for_R/R_Scripts_and_Input/EcoLetters_Archive_Zenodo/"

#' Data from seedling drydown experiment
drought <- read.csv(paste0(path,'Drought_Mortality.csv'))
drought_water <- read.csv(paste0(path,'Drought_WaterPots.csv'))
moisture_at_death <- read.csv(paste0(path,'Drought_Soil_Moisture_At_Death.csv'))


########################
#' 
#'   **Prepare Data**
#'
#'######################


#' *Seedling Mortality*
#' 

#'  View structure 
str(drought)

#'  Change species to a factor 
drought$species <- as.factor(drought$species)
unique(drought$species)

#'  Estimate means, standard deviations, and n
drought_means <- drought %>% 
  select(species, days_to_dead) %>%
  group_by(species ) %>% 
  summarise(days_to_dead_mean = mean(days_to_dead, na.rm=T), 
            days_to_dead_sd = sd(days_to_dead, na.rm=T),
            days_to_dead_n = length(days_to_dead))



#'  *Soil water loss*
#'  

#' View structure
str(drought_water)

#'   Change species to a factor 
drought_water$species <- as.factor(drought_water$species)
unique(drought_water$species)

#'  Estimate means
drought_water_means <- drought_water %>%
  select(-rep) %>%
  group_by(species) %>%
  summarise_all(mean, na.rm=T)

#' Separate join variables
drought_water_join <- drought_water_means %>% 
  select(species, water_loss_wk1_g_per_day)



#' *Estimated soil water content on day of death*
#' 

#' This dataframe contains an estimate of the mean soil moisture at death
#' for each species. To do this, we first estimated the mean drought survival
#' duration for each species ('drought_means$days_to_dead'). We then summarized 
#' an average pot drydown curve for each species, based on the mean gravimetric
#' soil moisture of three replicate pots on each day of the drydown experiment
#' (stored in 'drought_water_means'). In an externa workbook, we used these data
#' to estimate soil moisture at the fractional date of death, using linear 
#' interpolation between daily soil moisture estimates.
#' 
#' A few species had average drought survival times that were 1-2 days longer than the 
#' available soil moisture record (since we only recorded soil moisture in a subset
#' of 3 pots per species). In these cases, we modeled soil moisture on the average
#' day of death using a linear model of soil moisture loss based on the last two
#' soil moisture measurement days (again, typically 1-2 days prior to average day of death).

#' View structure
str(moisture_at_death)

#'   Change species to a factor 
moisture_at_death$species <- as.factor(moisture_at_death$species)
unique(drought_water$species)


#'  *Joining of dataframes*
#' 
#' Join summary dataframes of species' mean values for other analyses in R
list_df = list(drought_means, drought_water_join, moisture_at_death)
drought_df <- list_df %>% reduce(inner_join, by='species')


#'  *Preview drought data distributions*
#' 
drought_df %>% select(-species) %>% ggpairs()
#' 
#' Note that rate of water loss (water_loss_wk1_g_per_day) is not related to anything,
#'  nor does it vary much among species.
#'  
#' It appears that days_to_dead_mean is slightly left-skewed, so we will try ln_transformation
#' 
drought_df$ln_days_to_dead <- log(drought_df$days_to_dead_mean)
#' 
#' Check data - Looks better
drought_df %>% select(-species) %>% ggpairs()



#' *Export dataframe for use in full analysis*
#' 
write.csv(drought_df, 'drought_metrics.csv')



#'  *Create figure of soil water loss through the experiment, by species*
#' 
drought_fig_dat <- drought_water_means %>% select(-water_loss_wk1_g_per_day) %>%
  gather(X0:X18, key='day', value='g_water_per_g_soil')
drought_fig_dat$day_numeric <- gsub('X','', drought_fig_dat$day)
drought_fig_dat$day_numeric <- as.numeric(drought_fig_dat$day_numeric)

#' Specify the estimated wilting point (estimated water content when water potential = -1.5)
#'    
#'    This was estimated for our specific soil mix, using established protocols
#'    provided by METER Group for use with their WP4-C Water Potential Meter.
#'    This was the average estimate of gravimetric soil moisture at a water potential
#'    of -1.5 MPa, averaged across 3 replicates.
#'    
wilting_point <- 0.0289

#' *FIGURE S1-A*
#' 
#' Create figure with points showing soil moisture drydown curves for each species (in gray),
#'  highlighting soil moisture at the average day of death for each species (in black)
#' 
water_fig <- ggplot() + 
  geom_point(aes(x=day_numeric, y=g_water_per_g_soil), color='gray75', alpha=0.5, data= drought_fig_dat) + 
  geom_line(aes(x=day_numeric, y=g_water_per_g_soil, group=species), color='gray75', alpha=0.5, data= drought_fig_dat) +
  geom_point(data=drought_df, aes(x=days_to_dead_mean, y=mean_soil_moisture_at_death))+
  scale_shape_manual(values=c(15,1,16,17,0))+
  geom_text(aes(x=1.5,y=0.026,label="wilting point"), color='gray60') +
  geom_hline(aes(yintercept=wilting_point), lty=2, color='gray10') +
  labs(x='\nDay of Drought',y="Gravimetric Water Content\n(g water g-1 dry soil)\n",
       shape="Species", lty="Species", color="Species") +
  theme_classic()
water_fig

#' *FIGURE S1-A*
#' 
#' Similar figure, faceted by species
#' 
water_fig_b <- ggplot() + 
  geom_point(aes(x=day_numeric, y=g_water_per_g_soil), color='gray75', alpha=0.5, data= drought_fig_dat) + 
  geom_line(aes(x=day_numeric, y=g_water_per_g_soil, group=species), color='gray75', alpha=0.5, data= drought_fig_dat) +
  geom_point(data=drought_df, aes(x=days_to_dead_mean, y=mean_soil_moisture_at_death), cex=3.5)+
  geom_point(data=drought_df, aes(x=days_to_dead_mean, y=mean_soil_moisture_at_death, color=days_to_dead_mean), cex=3)+
  scale_color_gradientn(colors = c("ivory","#67001F"),
                        values = scales::rescale(c(9, max(drought_df[4]))))+
  labs(x='\nDay of Drought',y="Gravimetric Water Content\n(g water g-1 dry soil)\n",
       shape="Species", lty="Species", color="Drought\nSurvival\n(days)") +
  ylim(0,0.24) +
  geom_hline(yintercept=wilting_point, lty=2, color='gray10') +
  facet_wrap(~species, ncol=7) +
  theme_classic()
water_fig_b




