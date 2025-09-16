#' ---
#' title: "Estimating freezing tolerance"
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


#' **This code: Estimating freezing tolerance as lethal freezing temperature**
#' 
#' In this code, we use survival data from the freezing experiment to estimate 
#' 'freezing tolerance' as the lethal freezing temperature for each species (LT50).
#'  We also characterize freezing temperature dynamics during the experiment.



########################
#'
#' **Load or Source data**
#'
#'#####################

#' Specify working directory
path <- "C:/Users/larsonju/OneDrive - UW/Documents/Seed_Team/Data/Data_for_R/R_Scripts_and_Input/EcoLetters_Archive_Zenodo/"

#' Data from seedling freezing experiment
freeze_raw <- read.csv('Freezing_Mortality.csv')
freeze_temp <- read.csv('Freeze_Temps.csv')


########################
#'
#' **Load packages, options & functions**
#'
#'#####################

library(tidyverse)
library(GGally)
library(performance)
library(ggeffects)

options(contrasts = c("contr.sum", "contr.poly"))
select <- dplyr::select

########################
#' 
#'   **Freezing temperature data**
#'
#'######################


#' First, let's take a look at temperature data throughout the freezing experiment.
#' 

#' View data structure
#'    Note that data are separated by June A, June B, and July trials; we needed to
#'    utilize three different trials in order to assess all 49 species 
str(freeze_temp)

#' Gather data from different trials for a figure
freeze_temp_fig_dat <- freeze_temp %>% gather(Air_June_A:Soil_Average, key="Trial", value="Temp_C")

#' Subset only air temperature data for a figure
air_freeze_temp_fig_dat <- freeze_temp %>% 
  select(Phase:Air_Average) %>%
  gather(Air_June_A:Air_Average, key="Trial", value="Temp_C")


#' Subset air and soil data from day 7 to 8 only (during the freezing portion of experiment)
freezing_portion_fig_dat <- freeze_temp %>% gather(Air_June_A:Soil_Average, key="Trial", value="Temp_C") %>%
  filter( between(Day_Time,7,8))

#' *Figure S2-A*
#' 
#' Air temperatures throughout the whole experiment
#' 
air_temps <- ggplot( data=air_freeze_temp_fig_dat ) + 
  geom_point( aes(x=Day_Time, y=Temp_C, color=Trial), cex=1, alpha=0.4) +
  geom_line( aes(x=Day_Time, y=Temp_C, color=Trial), cex=1, alpha=0.4) +
  scale_color_manual(values=c("navyblue", "skyblue","dodgerblue","dodgerblue4"))+
  ylim(-12,30) + 
  labs(x="\nDay of Trial", y="Temperature (°C)\n") +
  annotate("text", x=3.5, y=-11,label='Growth',size=3) +
  annotate("text", x=6, y=-11,label='Hardening',size=3) +
  annotate("text", x=7.5, y=-11,label='Freeze',size=3) +
  annotate("text", x=9, y=-11,label='Reacclimation',size=3) +
  annotate("text", x=12.5, y=-11,label='Recovery',size=3) +
  theme_bw()
air_temps


#' *Figure S2-B*
#' 
#' Air and soil temps during the freezing experiment
#' 
freeze_period <- ggplot( data=freezing_portion_fig_dat) + 
  geom_point( aes(x=Day_Time, y=Temp_C, color=Trial), cex=1, alpha=0.4) +
  geom_line( aes(x=Day_Time, y=Temp_C, color=Trial), cex=1, alpha=0.4) +
  scale_color_manual(values=c("navyblue", "skyblue","dodgerblue","dodgerblue4",
                              'darkgoldenrod4','goldenrod3','goldenrod','gold'))+
  labs(x="\nDay of Trial", y="Temperature (°C)\n") +
  theme_bw()
freeze_period




########################
#' 
#'   **Estimating Freezing Tolerance**
#'
#'######################


#' *View dataframe*
#'  
str(freeze_raw)


#' *Check response data*
#'  
#' Check data to see whether certain categories of seedlings should be excluded, 
#'    specifically, seedings with:
#'    (1) abnormal development (e.g., a partially deformed or missing cotyledon), 
#'    (2) with leaves (grasses) or cotyledons (forbs) not exposed, meaning they had
#'        not yet emerged from the seed coat or coleoptile at the time of freezing
#'              -- For forbs, this only includes seedlings with seed coat >50% on
#'    (3) with small size (noted as much smaller than typical for a species) 
#'    (4) with true leaves already developing (forbs only)


#' Figure - Violin plot of survival according to freeze temperature and development category 
ggplot(data=freeze_raw) + 
  geom_violin(aes(x=development_category,y=survival, fill=development_category))+ 
  facet_wrap(~freeze_temp)+
  scale_fill_manual(values=c('darkred','darkorange3','gray85','blue','gray40')) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(x='\nDevelopment Category', fill='\nDevelopment Category', y='Survival')

# View n of each group
freeze_raw %>% select(freeze_temp, development_category, survival) %>%
  group_by(freeze_temp, development_category) %>%
  summarise(survival_n = length(survival)) 

#' Notes for each of the development categories:
#'  (1) Abnormal development: DO NOT EXCLUDE, Survival seems like it might be a little
#'                                            lower at -8C, but otherwise survival trends
#'                                            are similar to normal
#'  (2) Leaf not exposed:   EXCLUDE, Survival is really pretty different at most temperatures
#'  (3) Small size:  DO NOT EXCLUDE, Proportions may seem off by a lot but 
#'                                   the sample size of small seedlings is just much smaller
#'  (4) True leaves: DO NOT EXCLUDE, Very similar to 'normal'

#' Final dataset with seedlings removed where the leaf was not exposed
freeze <- freeze_raw %>% filter (development_category != "leaf not exposed")
str(freeze)



#' *Estimate LT50 for each species*
#' 


#' *Logistic regression: Single species*
#' 
#' Specify a logistic regression model to estimate survival as a function of temperature
#'     (following Meyer & Badaruddin 2001, Crop Sci.)


#' Here's some code to manually work through one species at a time
#' 
#' Select the species
#' Can use this code to set a different species quickly without changing the rest of the code below
test_species <- 'ACHMIL'

#' Subset the dataframe
species_freeze_test <- freeze %>% filter(species==test_species) %>%
  select(species,survival, freeze_temp, pot_id)
str(species_freeze_test)

#' Run model
freeze_mod <- glm(survival ~ freeze_temp, data=species_freeze_test, family=binomial)
summary(freeze_mod)
performance(freeze_mod)

#' Estimate LT50 and R2 values
mod_sum <- summary(freeze_mod)
a <- mod_sum$coefficients[1,1]
b <- mod_sum$coefficients[2,1]
freeze_LT50 <- -a/b
freeze_LT50
mod_perf <- performance(freeze_mod)
r2 <- mod_perf$R2_Tjur
r2

#' Test figure
freeze_predict <- ggpredict(freeze_mod, terms = c("freeze_temp"))
ggplot(freeze_predict, aes(x, predicted)) + 
  geom_point(color="dodgerblue3",cex=3) +
  geom_smooth(aes(x=freeze_temp, y=survival), method='glm',  method.args = list(family = "binomial"), 
              se = FALSE, color='dodgerblue3', data=species_freeze_test)+
  labs(x='\nFreezing temperature', y="Survival Probability\n") + 
  geom_point(aes(x=freeze_LT50,y=0.5),cex=4, color='darkblue') +
  geom_text(aes(x=freeze_LT50+0.15,y=0.45,label=round(freeze_LT50,3)),cex=4, color='darkblue') +
  geom_point(aes(x=freeze_temp, y=survival), position=position_jitter(height=0.05,width=0.05), 
             alpha=0.3, cex=2.5,data=species_freeze_test)+
  theme_bw()




#' *Logistic regression: All species*
#' 
#' After exploring models for several species, it is clear that we have too many singularity issues 
#' to include random effects, so we will utilize basic logistic regression.
#' 
#' Build a loop to cycle through all species

#' Setup
species_list <- unique(freeze$species)
lt50_logistic <- data.frame(matrix(NA,0,6), stringsAsFactors = F)
names(lt50_logistic) <- c("species","convergence","intercept","slope","R2","LT50")

#' For loop
for(i in species_list){
  model_df <- subset(freeze, species==i)
  logistic_model <- glm( survival ~ freeze_temp, family=binomial, data=model_df)
  mod_sum <- summary(logistic_model)
  mod_perf <- performance(logistic_model)
  lt50_logistic <- rbind(lt50_logistic, data.frame(species=i,
                                                   convergence= as.character(logistic_model$converged[1]),
                                                   intercept=mod_sum$coefficients[1,1], 
                                                   slope=mod_sum$coefficients[2,1], 
                                                   R2=mod_perf$R2_Tjur,
                                                   LT50= -1*mod_sum$coefficients[1,1]/mod_sum$coefficients[2,1])
  )
}

#' Save a rounded R2 value
lt50_logistic$R2_round <- round(lt50_logistic$R2, 2)



#' View results, including model convergence status
#' 
freeze_fig <-ggplot(data=freeze, aes(x=freeze_temp, y=survival)) + 
  geom_point(alpha=0.5, position=position_jitter(height=0.1, width=0.2),cex=1, alpha=0.3) +
  geom_smooth(method = 'glm',method.args = list(family = "binomial"), se = FALSE, color='gray20') + 
  geom_hline(aes(yintercept=0.5), color='black', lty=2, size=0.5) +
  geom_point(aes(x=LT50, y=0.5, color=convergence),cex=3, data=lt50_logistic) +
  scale_color_manual(values=c('red3','dodgerblue3')) + 
  labs(x="Temperature (C)", y="Survival Probability\n", color="Model Covergence") +
  geom_text(aes(x=-2.5, y=0.3,label=paste("R2: \n",R2_round)),color='black',size=2.75, data=lt50_logistic)+
  theme_bw() + 
  facet_wrap(~species)
freeze_fig


#' *Note on non-convergence:*
#'  Six of models did not converge, but these were all species with *Complete Separation* 
#'  -- at every temperature, ALL replicates either survived or died (1 or 0).
#'  
#'  LT50 is still reasonably predicted in these cases, so we will not worry about this warning.
#'
#' View non-convergence
lt50_logistic %>% filter(convergence==FALSE)  # View list of species



#' *Note on excluded species*
#' 
#' Two species stand out with poor and random survivorship across temperatures:
#'   CHADOU, SPOCRY
#'   
#'  We will not use LT50 estimates for these species (remove from dataset)
#'  
freeze2 <- freeze
lt50_logistic2 <- lt50_logistic
freeze2 <- freeze2 %>% filter (species != "CHADOU") %>% filter (species != "SPOCRY")
lt50_logistic2  <- lt50_logistic2 %>% filter (species != "CHADOU") %>% filter (species != "SPOCRY")


#' *Final -- Figure S3* 
#' 
#' Create final figure with excluded species removed.
#' Color points by lethal freezing temperature.
#' 
freeze_fig_update <-ggplot(data=freeze2, aes(x=freeze_temp, y=survival)) + 
  geom_point(alpha=0.5, position=position_jitter(height=0.1, width=0.2),cex=1, alpha=0.3) +
  geom_smooth(method = 'glm',method.args = list(family = "binomial"), se = FALSE, color='gray20') + 
  geom_hline(aes(yintercept=0.5), color='black', lty=2, size=0.5) +
  geom_point(aes(x=LT50, y=0.5),color='black',cex=3.5, data=lt50_logistic2) +
  geom_point(aes(x=LT50, y=0.5, color=LT50),cex=3, data=lt50_logistic2) +
  scale_color_gradientn(colors = c("azure","#053061"),
                        values = scales::rescale(c(-4, min(lt50_logistic2[6],na.rm=T) )),
                        na.value='white')+
  labs(x="Temperature (C)", y="Survival Probability\n", color="Freezing\nLethal\nTemp.\n(LT50,°C)") +
  geom_text(aes(x=-2.5, y=0.3,label=paste("R2: \n",R2_round)),color='black',size=2.75, data=lt50_logistic2)+
  theme_bw() + 
  facet_wrap(~species)
freeze_fig_update


#' *Save output*
#' 
#' View final structure of data for 47 species
str(lt50_logistic2)

#' Select and save LT50 variable (lethal freezing temp)
freeze_metric <- lt50_logistic2 %>% select(species, LT50)


