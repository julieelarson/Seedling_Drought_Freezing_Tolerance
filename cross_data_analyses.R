#' ---
#' title: "Cross-data analyses: Trait, climate and phylogenetic drivers of drought & freezing stress"
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


#' **This code: Conducting cross-data analyses**
#' 
#' In this code, we unite several species-level metrics generated in 
#' other R scripts to assess how traits, phylogenetic relatendness, and climate niche
#' explain variation in seedling drought and freezing tolerance. These synthetic analyses 
#' comprise most of the key figures and analyses included in the main text.
#' 
#' Most input data for analyses were generated in other R scripts available in 
#' this or other data repositories.
#' 
#' 
#' *Source of explanatory variables:*
#' 1) Seed and seedling traits: Linked from other repository 
#' 2) Phylogenetic PCoA scores: phylogenetic_metrics&analysis.R
#' 3) Climate niche metrics: climate_niche.R
#' 
#' *Source of response variables:* 
#' 1) Lethal freezing temperatures: freezing_tolerance.R
#' 2) Drought tolerance: drought.tolerance.R



#######################
#'
#' **Packages, options & functions**
#'
#'#####################

#' Packages
#' 
library(tidyverse)
library(GGally)
library(vegan)
library(patchwork)
library(rdacca.hp)
library(ggvenn)
library(ggrepel)
library(scales)
library(car)
library(PVR)
library(DescTools)
library(here)
library(hier.part)
library(rdryad)



#' Options 
#' 
options(contrasts = c("contr.sum", "contr.poly"))
select <- dplyr::select

#' Functions
#' 
#  Calculate standard errors
se <- function(x, na.rm=FALSE) {
    if (na.rm) x <- na.omit(x)
    sqrt(var(x)/length(x))
} 

# Change color of correlation plots in ggpairs()
# Reference:  https://stackoverflow.com/questions/45873483/ggpairs-plot-with-heatmap-of-correlation-values
ggpairs_fn <- function(data, mapping, method="p", use="pairwise", ...){
    # grab data
    x <- eval_data_col(data, mapping$x)
    y <- eval_data_col(data, mapping$y)
    # calculate correlation
    corr <- cor(x, y, method=method, use=use)
    # calculate colour based on correlation value
    # Here I have set a correlation of minus one to blue, 
    # zero to white, and one to red 
    # Change this to suit: possibly extend to add as an argument of `my_fn`
    colFn <- colorRampPalette(c("blue", "white", "red"), interpolate ='spline')
    fill <- colFn(100)[findInterval(corr, seq(-1, 1, length=100))]
    ggally_cor(data = data, mapping = mapping, ...) + 
        theme_void() +
        theme(panel.background = element_rect(fill=fill))
}





########################
#'
#' **Source and prepare data**
#'
#'#####################


#' *Specify working directory*
#' 
wd <- "C:/Users/larsonju/OneDrive - UW/Documents/Seed_Team/Data/Data_for_R/R_Scripts_and_Input/EcoLetters_Archive_Zenodo/"


#' *Phylogenetic metrics*
#' 
phylogenetic_raw <- read.csv(paste0(wd,'phylo_metrics.csv'))
#' Remove number column
phylogenetic <- phylogenetic_raw  %>% select (-X)



#' *Drought tolerance*
#' 
drought_raw <- read.csv(paste0(wd,"drought_metrics.csv"))
#' Subset drought tolerance only
drought <- drought_raw %>% select (species, days_to_dead_mean)



#' *Freeze Tolerance*
#' 
#' Read in file
freeze_raw <- read.csv(paste0(wd,"freeze_metric.csv"))
#' Subset freeze tolerance only
freeze <- freeze_raw %>% select(species, LT50)



#' *Species & Trait Data*
#' 
#' We will pull in a list of species and traits from archived data associated
#' with prior publication from this experiment:
#' 
#'  Larson, J., S. Copeland, D. Neuhuas. 2025. Seed and seedling traits suggest ontogenetic 
#'  coordination in the functional recruitment niche for dryland restoration species.
#'  Journal of Ecology.
#' 
#' Download data from Dryad into a local folder based on the dataset DOI
#'   Note that this downloads ALL files associated with the DOI,
#'   while we will only need one
data <- dryad_download(dois="10.5061/dryad.b5mkkwhpt")

#' View all of the downloaded files, then save the path to the last (8th) file, which 
#' contains species names, families, and other traits
data_df <- data.frame(data)
data_df
path <-data_df[8,1]
path

#' Read in the dataframe of species names, families and other traits
traits_raw <- read.csv(path)

#' View species data structure
str(traits_raw)

#' Save a species names dataframe
species_names <- traits_raw %>% select(species, species_name)


#' Subset to traits used in analysis
#' 
#' The dataset contains species names and families (columns 1:3) and many 
#' seed and seedling traits (columns 4:28). We will only use a subset of
#' traits selected a priori for these analyses based on expected links to
#' stress tolerance (seed mass, root diameter, SLA, RMR, RER Minimum Temperature)
#' and avoidance (seed dormancy index, germination minimum temperature), so we will 
#' subset those here along with species information.
#' 
#' Some traits have already been ln-transformed to improve normality.
#' 
traits <- traits_raw %>% select(species, species_name, family,
                                ln_seed_mass, 
                                ln_Rdiam, 
                                SLA, 
                                RMR, 
                                ln_RER, 
                                RER_min_temp,
                                seed_dorm_index,
                                germ_min_temp)

#' Missing trait values
#' 
#' Two species have missing germination minimum temperature values that were
#' imputed from a full trait dataset in Larson et al. 2025. We will use these
#' imputed values here:
#' 
#' LUPARG: -0.3282943 degC
#' SPOCRY: 8.4424510 deg C
#' 
traits[37,11] <- -0.3282943
traits[49,11] <- 8.4424510




#' *Climate niche metrics*
#' 
climate_raw <- read.csv(paste0(wd,'climate_niche_metrics.csv'))
#'  Adjust MAT values to deg C
climate_raw$MAT05 <- climate_raw$MAT05/10
#' Add a species code column
climate_join <- merge(climate_raw, species_names, by="species_name")
#' 
#'  Reduce to only MAT05 and MAP05 metrics which capture the coldest (MAT05, 5th percentile 
#'  mean annual temp) and driest (MAP 05, 5th percentile mean annual precip) conditions in
#'  each species' distribution, then rename the species column.
climate <- climate_join %>% select(species, MAT05, MAP05) 
#'
#' View data
ggplot(data=climate, aes(x=MAT05, y=MAP05)) + geom_point()
#' 
#' There is one clear outlying species, Asclepias tuberosa, which has been 
#' documented in the Colorado Plateau as well as most states east of the 
#' Mississippi River, leading to a very high mean annual precip (MAP05) estimate.
#' 
#' Rather than transform data, we've opted to remove this species from climate
#' analyses give the extreme deviation.
#' 
#' Remove MAP value
climate$MAP05[climate$MAP05>600] <- NA
#' Remove MAT value
climate$MAT05[climate$MAT05==6.7] <- NA




#' *Join dataframes*
#' 

#' Store dataframes in a list
df_list <- list(traits, climate, phylogenetic, drought, freeze)

#' Join all dataframes by 'ID' using left_join
species_dat_f <- df_list %>%
  reduce(left_join, by = "species")


#' *Add plant family variable*
#' 
#' Create a 'reduced' family variable that replaces all families outside of
#' Poaceae, Fabaceae, or Asteraceae with 'other'
#' 
unique(species_dat_f$family)  # View all unique families
`%notin%` <- Negate(`%in%`)  # Create function to refer to all values NOT in a list
family_list <- c('Poaceae', 'Fabaceae','Asteraceae')  # Create list of families to keep
#' 
#' Create reduced family column replacing anything that is not in family_list with 'other'
species_dat_f$family_reduced <- species_dat_f$family
species_dat_f$family_reduced[species_dat_f$family_reduced %notin% family_list] <- "Other"
species_dat_f$family_reduced <- as.factor(species_dat_f$family_reduced)






#'###
#'
#' **Trait PCA**
#'
#'###

#' We will first use PCA to reduce trait dimensionality to a limited number of
#' axes that can be used in analyses alongside climate and phylogenetic metrics.

#' Set rownames to store species codes
rownames(traits) <- traits$species

#' View trait distributions
species_traits_pca <- traits %>% select(-species, -species_name, -family) 
ggpairs(species_traits_pca)                            

#' Scale traits in trait dataframe
species_traits_pca_s <- scale(species_traits_pca)

#' Run PCA on trait dataframe
trait_pca <- rda(species_traits_pca_s, scale=FALSE)                         

#' View results
plot(trait_pca)
summary(trait_pca)
summary(trait_pca)$species

#' Save trait loading vectors and species scores
loadings <- data.frame(summary(trait_pca)$species)
species_scores <- data.frame(summary(trait_pca)$sites)
species_scores$species <- row.names(species_scores)

#' Save list of traits and PC axes
trait_list <- species_traits_pca  %>% colnames()
pc_list <- species_scores %>% select(-species) %>% colnames()
pc_trait_list <- c(trait_list, pc_list)

#' Bind species scores, traits, and species info
species_traits_pca$species <- rownames(species_traits_pca)
pca_bind <- merge(species_traits_pca, species_scores, by='species')

#' Estimate correlations between traits and PCs
pca_corr_dat <- pca_bind %>% select (all_of(pc_trait_list))
pca_corrs <- data.frame(cor(pca_corr_dat))[1:length(trait_list),] %>% select(all_of(pc_list))
pca_corrs$label <- rownames(pca_corrs)
pca_corrs #View correlations

#' Save nicer labels
pca_corrs$label <- dplyr::recode(pca_corrs$label, 
                                 "ln_seed_mass"="Seed mass",
                                 "seed_dormancy_index" ="Seed dormancy",
                                 "germ_min_temp_10_imputed"="Germination minimum temp.",
                                 "ln_RootDiam_18" ="Root diameter",
                                 "SLA_18"  ="SLA",
                                 "RMR_18"="RMR",
                                 "ln_RER_18" ="RER",
                                 "RER_min_temp"="RER Minimum Temp.")


#'*Figure: Trait PCA*
#'
#'
#' PCA - PC1 vs PC2
pca_fig_1_2 <- ggplot () +
  geom_segment(mapping=aes(x=0, y=0, xend=PC1, yend=PC2), linewidth=0.5, data=pca_corrs)+
  geom_text(aes(x=PC1, y=PC2,label=species), cex=2, color="gray35", data=pca_bind)+
  labs(x=paste("Trait PC1\n41.5% of variation\n"), y=paste("Trait PC2\n21.9% of variation\n")) +
  geom_text_repel(mapping=aes(x=PC1, y=PC2, label = label), size=2.5, fontface="bold", color="black", data=pca_corrs) +
  theme_classic() + 
  theme(text=element_text(size=10))
pca_fig_1_2


#' PCA - PC1 vs PC3
pca_fig_1_3 <- ggplot () +
  geom_segment(mapping=aes(x=0, y=0, xend=PC1, yend=PC3), linewidth=0.5, data=pca_corrs)+
  geom_text(aes(x=PC1, y=PC3,label=species), cex=2, color="gray35", data=pca_bind)+
  labs(x=paste("Trait PC1\n41.5% of variation\n"), y=paste("Trait PC3\n14.7% of variation\n")) +
  geom_text_repel(mapping=aes(x=PC1, y=PC3, label = label), size=2.5,fontface="bold", color="black", data=pca_corrs) +
  theme_classic() + 
  theme(text=element_text(size=10))
pca_fig_1_3

#' Joined figure
pca_fig_1_2 / pca_fig_1_3  + plot_annotation(tag_levels = 'A') &  theme(plot.tag = element_text(size=12))


#' *Merge Trait PC1, PC2, and PC3 with the full dataframe*
#' 
trait_pc_scores <- species_scores %>% select(species, PC1:PC3)
species_dat_full <- merge(species_dat_f, trait_pc_scores, by="species")





#' **Figure 2 in the MS**
#' 
#' **Combined Explanatory Variables**
#' 



#' Phylo
#' 
#' 
phylo_pcoa_1 <- ggplot (data=species_dat_full) + 
    geom_point(aes(x=phylo_pcoa1, y=phylo_pcoa2, shape=family_reduced, color=family_reduced), cex=2, position=position_jitter(width=0.002,height=0.002)) +
    labs(x='\nPhylo PCoA1 - 33.6%\n', y='\nPhylo PCoA2 - 14.5%\n', shape="Plant\nFamily", color="Plant\nFamily")+
    scale_shape_manual(values=c(15,16,17,18)) +
    scale_color_manual(values=c("gray20", "gray40", "gray60", "gray80")) +
    ggtitle('\nPhylogeny\n')+
    theme_classic() +
    theme(text=element_text(size=10), plot.title = element_text( face = "bold"))
phylo_pcoa_1



#' Traits
#'  
#'  
#' Save labels, but note that the sign of vectors can flip
#' when repeating the PCA
#' 
pca_corrs$PC1_label <- c(-1.11,-0.93,.95,.45,-.45,.19,-.38,.75)
pca_corrs$PC2_label <- c(-.04,.14,-.08,-.7,-.8,.52,.88,.37)
pca_corrs$PC3_label <- c(.23,.03,-.1,0.5,.4,.92,.16,.38)

#' PCA - PC1 vs PC2
pca_fig_1_2_basic <- ggplot () +
    geom_segment(mapping=aes(x=0, y=0, xend=PC1, yend=PC2), linewidth=0.5, data=pca_corrs)+
    geom_point(aes(x=PC1, y=PC2, shape=family_reduced, color=family_reduced), cex=2, data=species_dat_full)+
    labs(x=paste("Trait PC1\n41.5% of variation\n"), y=paste("Trait PC2\n21.9% of variation\n"), shape="Plant\nFamily", color="Plant\nFamily") +
    geom_text(mapping=aes(x=-PC1_label, y=-PC2_label, label = label), size=2.5, fontface="bold", color="black", data=pca_corrs) +
    scale_shape_manual(values=c(15,16,17,18)) +
    scale_color_manual(values=c("gray20", "gray40", "gray60", "gray80")) +
    theme_classic() + 
    ggtitle('\nTraits\n')+
    theme(text=element_text(size=10), legend.position = 'none', plot.title = element_text( face = "bold"))
pca_fig_1_2_basic


#' PCA - PC1 vs PC3
pca_fig_1_3_basic <- ggplot () +
    geom_segment(mapping=aes(x=0, y=0, xend=PC1, yend=PC3), linewidth=0.5, data=pca_corrs)+
    geom_point(aes(x=PC1, y=PC3, shape=family_reduced, color=family_reduced), cex=2,  data=species_dat_full)+
    labs(x=paste("Trait PC1\n41.5% of variation\n"), y=paste("Trait PC3\n14.7% of variation\n"), shape="Plant\nFamily", color="Plant\nFamily") +
    geom_text(mapping=aes(x=PC1_label, y=PC3_label, label = label), size=2.5,fontface="bold", color="black", data=pca_corrs) +
    scale_shape_manual(values=c(15,16,17,18)) +    #ylim(-1.75,1.75) +
    scale_color_manual(values=c("gray20", "gray40", "gray60", "gray80")) +
    #xlim(-1.75,1.75)+
    theme_classic() + 
    theme(text=element_text(size=10), legend.position = 'none')
pca_fig_1_3_basic



#' Climate
MAT_MAP_fig <- ggplot () +
    geom_point(aes(x=MAT05, y=MAP05, shape=family_reduced, color=family_reduced), cex=2,  data=species_dat_full)+
    labs(x=paste("\nMAT 05 (°C)"), y=paste("MAP 05 (mm)\n"), shape="Plant\nFamily", color="Plant\nFamily") +
    scale_shape_manual(values=c(15,16,17,18)) +   
    scale_color_manual(values=c("gray20", "gray40", "gray60", "gray80")) +
    ggtitle('\nClimate Niche\n')+
    theme_classic() + 
    theme(text=element_text(size=10), legend.position = 'none', plot.title = element_text( face = "bold"))
MAT_MAP_fig


pca_fig_1_2_basic+ pca_fig_1_3_basic   + guide_area()+ phylo_pcoa_1 + MAT_MAP_fig  + 
    plot_layout(guides = 'collect', widths=c(1,1,0.25))+ plot_annotation(tag_levels = 'A') &  theme(plot.tag = element_text(size=10))






#' **Analysis**
#' 
#' 


#'  First, let's view correlations among variables used in analysis.
#'

#' *Fig. S4*
#' 
#' *Trait correlations*
#' 
#' Prepare data with nice labels
species_dat_traits <- species_dat_full %>%
  select(ln_seed_mass, ln_Rdiam, SLA, , germ_min_temp,
         seed_dorm_index, ln_RER, RMR,
         RER_min_temp) %>%
  rename("Seed Mass (ln)" = "ln_seed_mass", 
         "Root Diam (ln)" = "ln_Rdiam",
         "SLA" = "SLA",
         "RMR" = "RMR", 
         "RER (ln)" = "ln_RER", 
         "RER Min. Temp." = "RER_min_temp",
         "Seed Dormancy" = "seed_dorm_index", 
         "Germ. Min. Temp." = "germ_min_temp")
#' Figure
species_dat_traits_corr_fig <- species_dat_traits  %>% 
  ggpairs(upper = list(continuous = ggpairs_fn), 
          lower = list(continuous = wrap("smooth", alpha=0.4, size=1.5)))
species_dat_traits_corr_fig + theme(strip.text.x = element_text(size = 6),
                                    strip.text.y = element_text(size = 6))




#' *Fig. S5*
#' 
#' *Explanatory and Response variables*
#' 

var_dat <- species_dat_full %>% select(days_to_dead_mean, LT50, PC1, PC2, PC3, 
                                       phylo_pcoa1, phylo_pcoa2, MAT05, MAP05) %>%
                         rename("Drought Survival" = "days_to_dead_mean",
                                "Freezing Lethal Temp." = "LT50",
                                "Trait PC1" = "PC1", 
                                "Trait PC2" = "PC2",
                                "Trait PC3" = "PC3",
                                "Phylo PCoA1" = "phylo_pcoa1",
                                "Phylo PCoA2" = "phylo_pcoa2",
                                "MAT 05" = "MAT05",
                                "MAP 05" = "MAP05")
# Figure
species_dat_fig <- var_dat %>% ggpairs(upper = list(continuous = ggpairs_fn), 
                                       lower = list(continuous = wrap("smooth", alpha=0.4, size=1.5)))
species_dat_fig + theme(strip.text.x = element_text(size = 6),
                          strip.text.y = element_text(size = 6))






#' 
#' 
#' **Drought ~ Freezing**
#' 
#' Are drought and freezing tolerance correlated? 
#' 
#' 


#' *Correlations*
#' 
#' View Pearson r correlation between drought and freezing
pearson_r <- cor.test(species_dat_full$LT50, species_dat_full$days_to_dead_mean)$estimate
pearson_r <- round(pearson_r,3)
pearson_r 
#' 
#'   Recall phylogenetic r, which was estimated via phylogenetic least squares analysis
#'   in the R script: phylogenetic_metrics&analysis.R
#'   
pgls_r <- -0.407


#' *Figure 1C in the MS*
#' 
#' 
drought_freeze <- ggplot (data = species_dat_full, aes(x=days_to_dead_mean, y=LT50)) + 
  geom_point (aes(color= family_reduced),cex=3, alpha=0.9) + 
  geom_smooth(method='lm', size=1, color='black', se=F) +
  geom_smooth(aes(color=family_reduced),method='lm', linetype=2, size=0.75, se=F) +
  scale_color_manual(values=c("gray20", "gray40", "gray60", "gray80")) +
  labs(x='\nDrought Survival (days)\n', y='\nFreezing Lethal Temp (°C)\n', color = 'Family') + 
  xlim(9,18) + 
  ylim(-9.5,-3.5)+
  geom_text(aes(x=16, y=-4, label=paste("Pearson r = ", pearson_r)), size=4.5)  +
  geom_text(aes(x=16, y=-4.25, label=paste("PGLS r = ", pgls_r)), size=4.5) +
  theme_bw() +
  theme(text=element_text(size=14), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
drought_freeze





#'
#' *Linear models, hierarchical and variance partitioning*
#' 
#' Constructing trait-only and full models to predict drought tolerance and freezing tolerance
#' 
#' 


#' *Scale explanatory variables*
#' 
#' Prior to analysis, put trait, phylogenetic and climatic variables on a similar scale
#' 
species_dat_full$MAT05_s <- scale(species_dat_full$MAT05)
species_dat_full$MAP05_s <- scale(species_dat_full$MAP05)
species_dat_full$phylo_pcoa1_s <- scale(species_dat_full$phylo_pcoa1)
species_dat_full$phylo_pcoa2_s <- scale(species_dat_full$phylo_pcoa2)
species_dat_full$PC1_s <- scale(species_dat_full$PC1)
species_dat_full$PC2_s <- scale(species_dat_full$PC2)
species_dat_full$PC3_s <- scale(species_dat_full$PC3)


  
#'  
#'  **Drought ~ Traits**
#'  
#'  How much variation in drought tolerance is explained by traits?
#'  
#'  


#' *Traits linear model*
#' 
#' 
d_mod_pc <- lm( days_to_dead_mean ~ PC1_s + PC2_s + PC3_s, data=species_dat_full)
summary(d_mod_pc)


#' *Traits hierarchical partitioning*
#' 
#' 
traits <- species_dat_full %>% select(PC1_s, PC2_s, PC3_s)  # Save trait dataframe
drought <-  species_dat_full$days_to_dead_mean
#' Set model
dt_rda <- rda(drought~.,traits)  # Same model as linear model above
RsquareAdj(dt_rda)               # R2 values are same as linear model
vif.cca(dt_rda)                  # Variance inflation factors (all 1 if orthogonal PCs)
#' Partitioning
dt_hp <- rdacca.hp(drought,traits,method="RDA", type="adjR2",var.part=T)
dt_hp



#' **Drought - All variables**
#'       
#'  How much additional variation do climate and phylogenetic variables explain?
#'  


#'  *Linear model*
#'   
d_mod_all <- lm( days_to_dead_mean ~ PC1_s + PC2_s + PC3_s + MAT05_s + MAP05_s + phylo_pcoa1_s + phylo_pcoa2_s, data=species_dat_full)
summary(d_mod_all)

#' 
#' *Hierarchical partitioning*
#' 
#' Remove missing data (row 7 - Asclepias tuberosa, ASCTUB) from analysis dataframe
species_dat_full_complete <- species_dat_full[-c(7), ] 
drought_complete <-  species_dat_full_complete$days_to_dead_mean

#' Capture explanatory and response variables
all_variables <- species_dat_full_complete %>% select(PC1_s,PC2_s,PC3_s,MAT05_s,MAP05_s, phylo_pcoa1_s, phylo_pcoa2_s)
drought <-  species_dat_full_complete$days_to_dead_mean

dall_rda <- rda(drought~.,all_variables)  # Same as linear model above
RsquareAdj(dall_rda)               # R2 values
vif.cca(dall_rda)                  # Variance inflation factors (all 1 if orthogonal PCs)
dall_hp <- rdacca.hp(drought,all_variables,method="RDA", type="adjR2",var.part=T)
dall_hp


#' **Create venn diagramm**
#' 
#' **Hierarchical partitioning by category**
#' 
#' 
traits_pc <- species_dat_full_complete %>% select(PC1_s, PC2_s, PC3_s)
climate <- species_dat_full_complete %>% select(MAT05_s, MAP05_s)
phylo <- species_dat_full_complete %>% select(phylo_pcoa1_s, phylo_pcoa2_s)

dall_hp_group <- rdacca.hp(drought_complete,list(traits_pc,climate,phylo),method="RDA", type="adjR2",var.part=T)
dall_hp_group

#' 
#' 
#' *Figure 4A*
#' 
#' Save variance partitioning list
d_var_part <- 100*dall_hp_group$Var.part[,1]
#' Save residual
d_residual <- 100-as.numeric(d_var_part[8])

d_var_part
d_residual

d_venn_list <- list( Traits = c(d_var_part[1],d_var_part[4], d_var_part[5],d_var_part[7]),
                     `Climate Niche` = c(d_var_part[2],d_var_part[4], d_var_part[6], d_var_part[7]),
                     Phylogeny = c(d_var_part[3], d_var_part[5], d_var_part[6], d_var_part[7]))

drought_venn <- ggvenn(d_venn_list, show_elements=TRUE, text_size=4.5) + 
    geom_text(aes(x=-1.25,y=-1.45,label=paste("Residual = ",d_residual)),cex=4) +
    geom_text(aes(x=.7,y=.65,label=d_var_part[2]),cex=4.5) +
    scale_fill_manual(values=c("#67001F",'tomato3','bisque2'))+ 
    ggtitle('Drought Survival\n')+
    theme(plot.title = element_text(hjust = 0.5, size=16))
drought_venn 





#'  
#'  
#'  **Freezing ~ Traits**
#'  
#'  How much variation in freezing tolerance is explained by traits?
#'  
#'  

#' *Traits linear model*
#'
f_mod_pc <- lm(  LT50 ~ PC1_s + PC2_s + PC3_s, data=species_dat_full)
summary(f_mod_pc)



#' **Traits hierarchical partitioning**
#' 
#' Filter out species with missing data
species_dat_full_freeze <- species_dat_full %>% filter(species != "CHADOU") %>% 
                                                filter(species != "SPOCRY")
traits_f <- species_dat_full %>% filter(species != "CHADOU") %>% 
                                filter(species != "SPOCRY") %>% 
                                select(PC1_s, PC2_s, PC3_s)  # Save trait dataframe
freezing <-  species_dat_full %>% filter(species != "CHADOU") %>% 
                                 filter(species != "SPOCRY") %>% 
                                 select(LT50)
#' Set model
ft_rda <- rda(freezing~.,traits_f)  # Same model as linear model above
RsquareAdj(ft_rda)               # R2 values are same as linear model
vif.cca(ft_rda)                  # Variance inflation factors (all 1 if orthogonal PCs)
#' Partitioning
ft_hp <- rdacca.hp(freezing,traits_f,method="RDA", type="adjR2",var.part=T)
ft_hp




#' **Freezing - All variables**
#' 
#'  How much additional variation do climate and phylogenetic variables explain?
#'  

#'       
#'   *Linear model*
#'   
species_dat_full_freeze <- species_dat_full %>% na.exclude()
f_mod_all <- lm(LT50 ~ PC1_s + PC2_s + PC3_s + MAT05_s + MAP05_s + phylo_pcoa1_s + phylo_pcoa2_s, data=species_dat_full_freeze)
summary(f_mod_all)


#' 
#' *Hierarchical partitioning*
#' 
#' Remove missing data (row 7 - Asclepias tuberosa, ASCTUB) from analysis dataframe
species_dat_freeze <- species_dat_full_freeze[-c(7), ] 

#' 
#' Model
all_variables_f <- species_dat_freeze %>% select(PC1_s,PC2_s,PC3_s, MAT05_s,MAP05_s, 
                                                           phylo_pcoa1_s, phylo_pcoa2_s)
freezing_f <- species_dat_freeze %>% select(LT50)
fall_rda <- rda(freezing_f~.,all_variables_f)  # Same as linear model above
RsquareAdj(fall_rda)               # R2 values
vif.cca(fall_rda)                  # Variance inflation factors (all 1 if orthogonal PCs)
fall_hp <- rdacca.hp(freezing_f,all_variables_f,method="RDA", type="adjR2",var.part=T)
fall_hp


#' **Create venn diagram**
#' 
#' 
#' *Hierarchical partitioning by category*
#' 
traits_pc <- species_dat_full_freeze %>% na.exclude %>% select(PC1_s, PC2_s, PC3_s)
climate <- species_dat_full_freeze %>% na.exclude %>% select(MAT05_s,MAP05_s)
phylo <- species_dat_full_freeze %>% na.exclude  %>% select(phylo_pcoa1_s, phylo_pcoa2_s)
fall_hp_group <- rdacca.hp(freezing_f,list(traits_pc,climate,phylo),method="RDA", type="adjR2",var.part=T)
fall_hp_group

#' 
#' 
#' *Figure 4B in the MS*
#' 
#' Save variance partitioning list
f_var_part <- 100*fall_hp_group$Var.part[,1]
#' Save residual
f_residual <- 100-as.numeric(f_var_part[8])

f_var_part
f_residual 

f_venn_list <- list( Traits = c(f_var_part[1],f_var_part[4], f_var_part[5],f_var_part[7]),
                     `Climate Niche`= c(f_var_part[2],f_var_part[4], f_var_part[6], f_var_part[7]),
                     Phylogeny = c(f_var_part[3], f_var_part[5], f_var_part[6], f_var_part[7]))

freeze_venn <- ggvenn(f_venn_list, show_elements=TRUE, text_size=4.5) + 
    geom_text(aes(x=-1.25,y=-1.45,label=paste("Residual = ",f_residual)),cex=4) +
    scale_fill_manual(values=c('navyblue','royalblue','azure2')) + 
    ggtitle('Freezing Lethal Temp.\n')+
    theme(plot.title = element_text(hjust = 0.5, size=16))
freeze_venn 








#' **Figure 3 -- Trait bivariate plots**
#' 
#' 


#' *Drought* 

#' Rename variables
#' 
#' 
species_dat_all <- species_dat_full %>%
  rename("Seed Mass (ln)" = "ln_seed_mass", 
         "Root Diam (ln)" = "ln_Rdiam",
         "SLA" = "SLA",
         "RMR" = "RMR", 
         "RER (ln)" = "ln_RER", 
         "RER Min. Temp." = "RER_min_temp",
         "Seed Dormancy" = "seed_dorm_index", 
         "Germ. Min. Temp." = "germ_min_temp",
         "Drought Survival" = "days_to_dead_mean",
         "Freezing Lethal Temp." = "LT50")


#' PC1 traits
#' 
#' Seed Mass
d_sm_fig <- ggplot(species_dat_all, aes(x=`Seed Mass (ln)`, y=`Drought Survival`))+ 
  geom_point(aes(color=`Drought Survival`), cex=2) +
  geom_smooth(method='lm',se=F,color='black') + 
  ylim(9,18)+
  scale_color_gradientn(colors = c("bisque2","#67001F")) +
  #geom_text(aes(y=17.75,x=-8.5, label=paste("r = ", round(d_sm,3))),size=3.5) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
d_sm_fig                            

#' SLA
d_sla_fig <- ggplot(species_dat_all, aes(x=SLA, y=`Drought Survival`))+ 
  geom_point(aes(color=`Drought Survival`),cex=2) +
  geom_smooth(method='lm',se=F,color='black') + 
  #geom_text(aes(y=17.75,x=900, label=paste("r = ", round(d_sla,3))),size=3.5) +
  scale_color_gradientn(colors = c("bisque2","#67001F")) +
  theme_bw() +
  ylim(9,18)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title.y = element_blank())
d_sla_fig  

#' Germ min temp
d_germ_fig <- ggplot(species_dat_all, aes(x=`Germ. Min. Temp.`, y=`Drought Survival`))+ 
  geom_point(aes(color=`Drought Survival`),cex=2) +
  geom_smooth(method='lm',se=F,color='black') + 
  #geom_text(aes(y=17.75,x=9.5, label=paste("r = ", round(d_germ,3))),size=3.5) +
  theme_bw() +
  scale_color_gradientn(colors = c("bisque2","#67001F")) +
  ylim(9,18)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title.y = element_blank())
d_germ_fig  

#' Root Diam
d_rdiam_fig <- ggplot(species_dat_all, aes(x=`Root Diam (ln)`, y=`Drought Survival`))+ 
  geom_point(aes(color=`Drought Survival`),cex=2) +
  geom_smooth(method='lm',se=F,color='black') + 
  #geom_text(aes(y=17.75,x=-1.5, label=paste("r = ", round(d_rdiam,3))),size=3.5) +
  theme_bw() +
  scale_color_gradientn(colors = c("bisque2","#67001F")) +
  ylim(9,18)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title.y = element_blank())
d_rdiam_fig  


#' PC2 traits
#' 

#' RER
d_rer_fig <- ggplot(species_dat_all, aes(x=`RER (ln)`, y=`Drought Survival`))+ 
  geom_point(aes(color=`Drought Survival`),cex=2) +
  geom_smooth(method='lm',se=F,color='black') + 
  #geom_text(aes(y=17.75,x=-2.6, label=paste("r = ", round(d_rer,3))),size=3.5) +
  theme_bw() +
  scale_color_gradientn(colors = c("bisque2","#67001F")) +
  ylim(9,18)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
d_rer_fig  


#' Seed Dormancy
d_dorm_fig <- ggplot(species_dat_all, aes(x=`Seed Dormancy`, y=`Drought Survival`))+ 
  geom_point(aes(color=`Drought Survival`),cex=2) +
  geom_smooth(method='lm',se=F,color='black') + 
  #geom_text(aes(y=17.75,x=0.25, label=paste("r = ", round(d_dorm,3))),size=3.5) +
  theme_bw() +
  scale_color_gradientn(colors = c("bisque2","#67001F")) +
  ylim(9,18)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title.y = element_blank())
d_dorm_fig  

#' RMR
d_rmr_fig <- ggplot(species_dat_all, aes(x=RMR, y=`Drought Survival`))+ 
  geom_point(aes(color=`Drought Survival`),cex=2) +
  geom_smooth(method='lm',se=F,color='black') + 
  #geom_text(aes(y=17.75,x=0.15, label=paste("r = ", round(d_rmr,3))),size=3.5) +
  ylim(9,18)+
  scale_color_gradientn(colors = c("bisque2","#67001F")) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title.y = element_blank())
d_rmr_fig  


#' PC3
#' 
#' RER min temp
d_rermin_fig <- ggplot(species_dat_all, aes(x=`RER Min. Temp.`, y=`Drought Survival`))+ 
  geom_point(aes(color=`Drought Survival`),cex=2) +
  geom_smooth(method='lm',se=F,color='black') + 
  #geom_text(aes(y=17.75,x=6.5, label=paste("r = ", round(d_rermin,3))),size=3.5) +
  ylim(9,18)+
  scale_color_gradientn(colors = c("bisque2","#67001F")) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
d_rermin_fig  



#' *Drought patchwork figure*
#' 
d_sm_fig + d_rdiam_fig + d_germ_fig+ d_sla_fig    + 
  d_rer_fig+ d_rmr_fig + d_dorm_fig  + plot_spacer() + 
  d_rermin_fig + guide_area() + plot_spacer() + plot_spacer() + 
  plot_layout(nrow = 3, widths=c(1.03,1,1,1), guides = 'collect')





#' **Freezing -- trait figure**
#' 
#' PC1 traits
#' 
#' Germ min temp
f_germ_fig <- ggplot(species_dat_all, aes(x=`Germ. Min. Temp.`, y=`Freezing Lethal Temp.`))+ 
    geom_point(aes(color=`Freezing Lethal Temp.`),cex=2) +
    geom_smooth(method='lm',se=F,color='black') + 
    #geom_text(aes(y=-3.5,x=1.5, label=paste("r = ", round(f_germ,3))),size=3.5) +
    theme_bw() +
    scale_color_gradientn(colors = c("#053061","azure2"))+
    ylim(-10,-3)+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title.y = element_blank())
f_germ_fig 


#' Seed Mass
f_sm_fig <- ggplot(species_dat_all, aes(x=`Seed Mass (ln)`, y=`Freezing Lethal Temp.`))+ 
    geom_point(aes(color=`Freezing Lethal Temp.`),cex=2) +
    geom_smooth(method='lm',se=F,color='black') + 
    ylim(-10,-3)+
    scale_color_gradientn(colors = c("#053061","azure2"))+
    #geom_text(aes(y=-3.5,x=-4.5, label=paste("r = ", round(f_sm,3))),size=3.5) +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
f_sm_fig                            

#' SLA
f_sla_fig <- ggplot(species_dat_all, aes(x=SLA, y=`Freezing Lethal Temp.`))+ 
    geom_point(aes(color=`Freezing Lethal Temp.`),cex=2) +
    geom_smooth(method='lm',se=F,color='black') + 
   # geom_text(aes(y=-3.5,x=300, label=paste("r = ", round(f_sla,3))),size=3.5) +
    theme_bw() +
    ylim(-10,-3)+
    scale_color_gradientn(colors = c("#053061","azure2"))+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title.y = element_blank())
f_sla_fig  


#' Root Diam
f_rdiam_fig <- ggplot(species_dat_all, aes(x=`Root Diam (ln)`, y=`Freezing Lethal Temp.`))+ 
    geom_point(aes(color=`Freezing Lethal Temp.`),cex=2) +
    geom_smooth(method='lm',se=F,color='black') + 
    #geom_text(aes(y=-3.5,x=-.14, label=paste("r = ", round(f_rdiam,3))),size=3.5) +
    theme_bw() +
    ylim(-10,-3)+
    scale_color_gradientn(colors = c("#053061","azure2"))+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title.y = element_blank())
f_rdiam_fig  


#' PC2 traits
#' 
#' RER
f_rer_fig <- ggplot(species_dat_all, aes(x=`RER (ln)`, y=`Freezing Lethal Temp.`))+ 
   geom_point(aes(color=`Freezing Lethal Temp.`),cex=2) +
    geom_smooth(method='lm',se=F,color='black') + 
    #geom_text(aes(y=-3.5,x=-0.3, label=paste("r = ", round(f_rer,3))),size=3.5) +
    theme_bw() +
    ylim(-10,-3)+
    scale_color_gradientn(colors = c("#053061","azure2"))+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
f_rer_fig  


#' Seed Dormancy
f_dorm_fig <- ggplot(species_dat_all, aes(x=`Seed Dormancy`, y=`Freezing Lethal Temp.`))+ 
     geom_point(aes(color=`Freezing Lethal Temp.`),cex=2) +
    geom_smooth(method='lm',se=F,color='black') + 
    #geom_text(aes(y=-3.45,x=0.2, label=paste("r = ", round(f_dorm,3))),size=3.5) +
    theme_bw() +
    ylim(-10,-3)+
    scale_color_gradientn(colors = c("#053061","azure2"))+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title.y = element_blank())
f_dorm_fig  

#' RMR
f_rmr_fig <- ggplot(species_dat_all, aes(x=RMR, y=`Freezing Lethal Temp.`))+ 
    geom_point(aes(color=`Freezing Lethal Temp.`),cex=2) +
    geom_smooth(method='lm',se=F,color='black') + 
   # geom_text(aes(y=-3.5,x=0.15, label=paste("r = ", round(f_rmr,3))),size=3.5) +
    ylim(-10,-3)+
    scale_color_gradientn(colors = c("#053061","azure2"))+
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title.y = element_blank())
f_rmr_fig  


#' PC3 traits
#' 
#' RER Min temp
f_rermin_fig <- ggplot(species_dat_all, aes(x=`RER Min. Temp.`, y=`Freezing Lethal Temp.`))+ 
    geom_point(aes(color=`Freezing Lethal Temp.`),cex=2) +
    geom_smooth(method='lm',se=F,color='black') + 
    #geom_text(aes(y=-3.5,x=1.5, label=paste("r = ", round(f_rermin,3))),size=3.5) +
    ylim(-10,-3)+
    scale_color_gradientn(colors = c("#053061","azure2"))+
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
f_rermin_fig  


                           
#' *Freezing patchwork figure*
#' 
#'  All three rows
f_sm_fig  +f_rdiam_fig + f_germ_fig + f_sla_fig + 
    f_rer_fig + f_rmr_fig + f_dorm_fig + plot_spacer() + 
    f_rermin_fig + guide_area() + plot_spacer() + plot_spacer() + 
    plot_layout(nrow = 3, widths=c(1.03,1,1,1), guides = 'collect')



