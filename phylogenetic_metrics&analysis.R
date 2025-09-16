#' ---
#' title: "Generating phylogenetic metrics for analysis of drought & freezing stress"
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


#' **This code: Conducting phylogenetic analyses**
#' 
#' In this code, we construct a phylogenetic tree and use this to estimate 
#' phylogenetic signal for lethal freezing temperature and drought survival duration.
#' 
#' We also estimate phylogenetic axes which capture major dimensions of phylogenetic
#' separation among study species (Phylo PCoS1 and PCoS2)




########################
#'
#' ** Load packages, options & functions **
#'
#'#####################
library(rdryad)
library(tidyverse)
library(phytools)
library(V.PhyloMaker2)  # Note -- must be installed via devtools
library(ggtree)
library(caper)
library(PVR)

select <- dplyr::select

########################
#'
#' **Source data**
#'
#'#####################


#' *Specify working directory*
#' 
wd <- "C:/Users/larsonju/OneDrive - UW/Documents/Seed_Team/Data/Data_for_R/R_Scripts_and_Input/EcoLetters_Archive_Zenodo/"



#' *Drought Tolerance*
#' 
#' Read in file
drought_raw <- read.csv("drought_metrics.csv")

#' View data
str(drought_raw)

#' Subset drought tolerance only
drought <- drought_raw %>% select (species, days_to_dead_mean)



#' *Freeze Tolerance*
#' 
#' Read in file
freeze_raw <- read.csv("freeze_metric.csv")

#' View data
str(freeze)

#' Subset freeze tolerance only
freeze <- freeze_raw %>% select(species, LT50)


#' *Species Info*
#' 
#' 
#' We will pull in a list of species names and families from archived data associated
#' with prior publication from this experiment.
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
species_dat <- read.csv(path)

#' View species data structure
str(species_dat)

#' Subset species info
#' 
#' The dataset contains species names and families (columns 1:3) and many 
#' seed and seedling traits (columns 4:28). We only need the species names 
#' and families, so we'll subset those here.
#' 
names_wp <- species_dat %>% 
  select(species_name, family)




#'###  
#' 
#'  **Create phylogenetic tree for species**
#'
#'
#'  Creating a phylogeny with PhlyoMaker V2
#'   
#'   https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9363651/
#'   
#'   Jin Y, Qian H. V.PhyloMaker2: An updated and enlarged R package that can generate 
#'   very large phylogenies for vascular plants. Plant Divers. 2022 May 27;44(4):335-339. 
#'   doi: 10.1016/j.pld.2022.05.005. 
#'   
#'   Erratum in: Plant Divers. 2022 Nov 26;45(1):122. PMID: 35967255; PMCID: PMC9363651.
#'   
#' 
#' Based on the botanical nomenclature of the World Plants (WP) database (https://www.worldplants.de)
#' Ver. 18.2  Dec. 5th, 2023
#' 

#' View species names
#' 
names_wp$species_name

#' Species names must match the World Plants database, so we will add a column 
#' for updated names that match (names checked manually, via the website)
#' 
names_wp$species_name_WP <- names_wp$species_name
names_wp$species_name_WP <- dplyr::recode (names_wp$species_name_WP, 
                                           'Achnatherum hymenoides' = 'Eriocoma hymenoides',
                                           'Pleuraphis jamesii'='Hilaria jamesii' ,
                                           'Cleome serrulata'='Cleomella serrulata' ,
                                           'Cleome lutea'  ='Cleomella lutea' ,
                                           'Machaeranthera canescens'='Dieteria canescens')
#' Add a genus column
names_wp$genus <-sub(" .*", "", names_wp$species_name_WP) 

#' Order and rename columns
names_wp <- names_wp %>% select(species_name_WP, genus, family)
colnames(names_wp)[1] <- "species"
colnames(names_wp)[2] <- "genus"
str(names_wp)



#' Generate a phylogeny for the species list
#' 
#'    Options are set so tree is based on World Plant database names for the tree and nodes arguments
#'    We use 'S3' for 'scenarios' argument, which is the most commonly used option among S1, S2, or S3.
#'    
#'    Option 'S3' is described in Qin and Jiang 2019:
#'    Jin, Y., & Qian, H. (2019). V.PhyloMaker: an R package that can generate very large phylogenies 
#'    for vascular plants. Ecography (Copenhagen), 42(8), 1353-1359. https://doi.org/10.1111/ecog.04434
#'    
#'    Option nodes.info.1 specifies to extract the genus- or family-level largest 
#'          cluster's root and basal node info from the mega tree (rather than separate phylogeny)
#'    
tree_wp <- phylo.maker(sp.list = names_wp, tree = GBOTB.extended.WP, nodes = nodes.info.1.WP, scenarios = "S3")
write.tree(tree_wp$scenario.3, "phylomaker_wp.tree")


#' Plot tree
tree_wp_read <- read.tree("phylomaker_wp.tree")
plot(tree_wp_read, cex = .8)


#'  Update tip labels with chosen names for final study reporting
tree_final <- tree_wp_read
recode <- dplyr::recode
tree_final$tip.label <- recode(tree_final$tip.label, 
                               'Eriocoma_hymenoides' = 'Achnatherum_hymenoides',
                               'Hilaria_jamesii' = 'Pleuraphis_jamesii',
                               'Cleomella_serrulata' = 'Cleome_serrulata',
                               'Cleomella_lutea' = 'Cleome_lutea',
                               'Dieteria_canescens' = 'Machaeranthera_canescens')


#'  Save tip names (with no underscores between genus and species) to create a nicely formatted tree
tree_final2 <-tree_final
tree_final2$tip.label <- gsub("_", " ", tree_final$tip.label)





#'###  
#' 
#'  **Estimate phylogenetic signal**
#'  for freeze and drought tolerance
#'

#' Create single dataframe of response metrics
#' 
species_dat_phylosig <- merge(drought, freeze, by="species", all=T)
#' Rename variables
species_dat_phylosig <- species_dat_phylosig %>%
  rename("Drought Survival" = "days_to_dead_mean",
         "Freezing Lethal Temp." = "LT50")
#' View data 
str(species_dat_phylosig)


#'  Save species names that match the tree into trait dataframe
#'  
tree_names <- data.frame(tree_names=tree_final$tip.label)
tree_names$species_name <- gsub("_", " ", tree_names$tree_names)
trait_names <- species_dat %>%  select(species, species_name)
final_names <- merge(tree_names, trait_names, by="species_name")


#' Save dataframe with species rownames only
dat_phylo_nice <- merge(final_names, species_dat_phylosig, by="species")
rownames(dat_phylo_nice) <- dat_phylo_nice$tree_names
dat_for_signal <- dat_phylo_nice %>% select(-species, -species_name, -tree_names)



#############################################
#' Match traits using tool from Traits textbook's R materials:
#' 
#' The tool "match_data" comes from the book 
#' "Handbook of Trait-based Ecology: From Theory to R Tools" 
#' by de Bello et al (2021).
#' 
source(here::here("C:/Users/larsonju/OneDrive - UW/Documents/Seed_Team/Data/Data_for_R/R_Scripts_and_Input/EcoLetters_Archive_Zenodo/match_data.R"))
match_traits_phylo <- match_data(traits = dat_for_signal, tree = tree_final) #links phylogeny with trait data
#############################################

TraitSignals<- data.frame(matrix(ncol = 3, nrow = 0)) #make blank dataframe to accept blomberg K data
colnames(TraitSignals) <- c('Trait', 'lambda', 'p') # name columns
TraitSignals[,2]<- as.numeric(TraitSignals[,2]) #set data type for lambda

for (i in 1:length(colnames(match_traits_phylo$traits))) { 
  TraitName<-colnames(match_traits_phylo$traits)[i] #trait name 
  traity<- as.matrix(match_traits_phylo$traits)[,i] #extract trait means per species
  #SEofTrait<-as.matrix(traits_phylo_se)[,i] #extract trait SE per species
  sig <- phylosig(match_traits_phylo$tree, traity, method = "lambda", test = TRUE) # Not included, SE code:  se = switch( is.na(SEofTrait[1]) == TRUE , NULL , SEofTrait  ) 
  plot.phylosig(sig)
  TraitSignals[i,] <- c(TraitName, round(as.numeric(sig$lambda),3), round(as.numeric(sig$P),3)) #add to dataframe
}
TraitSignals 




#'###  
#' 
#' **Estimate phylogenetic correlation**
#' between drought and freezing tolerance
#' 


#' Create a copy of our tree and set the nodes of the tree to NULL so we don't get an error about node
#' labels and tip labels being duplicated
tree_pgls <- tree_final
tree_pgls$node.label<-NULL

#' Second, create a copy of our response dataframe and add a species column
#' 
dat_pgls <- dat_phylo_nice

#   Add a species column, needed for pgls setup
dat_pgls$Species <- rownames(dat_pgls) 

#'   Change the order of the dataframe to match the tree
tree_order <- match_traits_phylo$tree$tip.label
dat_pgls <- dat_pgls[match(tree_order, dat_pgls$Species),]

#' Remove columns we don't need
dat_pgls <- dat_pgls %>% select(-species_name,-tree_names)

#' Third,
#' Specify the tree and trait data sources for PGLS models
pgls_setup <- comparative.data(tree_pgls, dat_pgls, Species, vcv=TRUE, vcv.dim=3)


#' Run PGLS model
#' 
#' Note that drought survival has the larger lambda value, so we will use
#' this as the response variable in PGLS.
#' 
trial_mod_1 <- pgls( `Drought Survival` ~ `Freezing Lethal Temp.`, pgls_setup, lambda = 'ML')
trial_mod_1
sum_trial_mod_1 <-summary(trial_mod_1)
sum_trial_mod_1 


#' Estimate phylogenetic correlation
#' 
#' Squareroot of multiple R2 value, multiplied by negative 1 to account for the negative
#' relationship (beta estimate)
pgls_r <- round(sqrt(sum_trial_mod_1$r.squared)*-1,3)

#' Compare to normal linear model
pearson_r <- round(-1*sqrt(summary(lm(`Drought Survival` ~ `Freezing Lethal Temp.`, data=dat_pgls))$r.squared),3)

pgls_r
pearson_r




#'###  
#'  
#'  **Phylogenetic heatmap**
#'  for freeze and drought tolerance *
#'
#' Species tree
tree_fig <- ggtree(tree_final2, ladderize=F) + 
  geom_tiplab(size=3, fontface=3) + 
  ggplot2::xlim(0,200)
tree_fig


#' Order species by phylogenetic tree
dat_for_signal$species <- row.names(dat_for_signal)
dat_for_signal <- dat_for_signal[match(tree_final$tip.label, dat_for_signal$species),]


#' Heatmap for drought
#' 
dat_drought <- dat_for_signal %>% select(`Drought Survival`, species)
dat_drought$trait <- "Drought Survival"
dat_drought$species <- factor(dat_drought$species, levels=(dat_drought$species))
#
#
# Use red gradient for drought
ggplot(dat_drought, aes( x=trait, y=species, fill = `Drought Survival`)) +
  geom_tile(colour="black", size=0.05) + 
  scale_x_discrete(position = "top") +
  scale_fill_gradientn(colors = c("ivory","#67001F"),
                       values = scales::rescale(c(9, max(dat_drought[1]))))+
  labs(x='Drought', fill="Drought\nSurvival\n(days)", y="Species") + 
  theme_classic() +
  theme(axis.title=element_text(size=6), axis.text.x = element_blank(), legend.text = element_text(size=6), legend.title = element_text(size=6)) 



#' Heatmap for freezing
#'
dat_freezing <- dat_for_signal %>% select(`Freezing Lethal Temp.`,species)
dat_freezing$trait <- "Freezing Lethal Temp."
dat_freezing$species <- factor(dat_freezing$species, levels=(dat_freezing$species))
#
#
#' Use blue gradient for freezing
ggplot(dat_freezing, aes( x=trait, y=species, fill = `Freezing Lethal Temp.`)) +
  geom_tile(colour="black", size=0.05) + 
  scale_x_discrete(position = "top") +
  scale_fill_gradientn(colors = c("azure","#053061"),
                       values = scales::rescale(c(-4, min(dat_freezing[1],na.rm=T) )),
                       na.value='white')+
  labs(x='Freezing', fill='Freezing\nLethal\nTemp\n(°C)', y="Species") + 
  theme_classic() +
  theme(axis.title=element_text(size=6), axis.text.x = element_blank(), legend.text = element_text(size=6), legend.title = element_text(size=6)) 





#'###  
#' 
#'  **Generate PCoA axes based on species' relatedness**
#'
#' We use the approach described in Rosbakh et al. (2022), using ordination to generate
#' orthogonal components that capture species' relatedness based on a phylogenetic tree.
#' 
#' Citation:
#' RosBakh, S., M. Pichler, P. Poschlod. 2022. Machine-learning algorithms predict 
#' soil seed bank persistence from easily available traits. Applied Vegetaiton Science 25:e12660
#' 
#' Paper:
#' https://doi.org/10.1111/avsc.12660
#' 
#' R Code:
#' https://github.com/MaximilianPi/Rosbakh-Pichler-Poschlod-2021/blob/main/Code/phylogeny.R
#' 
#' The approach utilizes the package PVR:
#' https://cran.r-project.org/web/packages/PVR/PVR.pdf#page=9.08
#' 

#' *Prepare Tree*
#' 
#' View phylogenetic tree
tree_final2

#' Solve multichotomies (only one among three Dalea species)
tree_eigen <- multi2di(tree_final2)
ggtree(tree_eigen, ladderize=F) + 
  geom_tiplab(size=3, fontface=3) + 
  ggplot2::xlim(0,200)



#' *PVR approach from Rosbackh et al. 2022*
#' 
pvr1 <- PVRdecomp(tree_eigen)

#' Extract eigenvalues and estimate percent variance explained
#'    
#'    Component 1 explains 33.6%
#'    Componenet 2 explains 14.5%
#'    
phylo_eigvals <- pvr1@Eigen$values
phylo_eigperc <- phylo_eigvals/(sum(phylo_eigvals)) 
phylo_eigperc

#' View the cumulative variance explained
phylo_eigperc_cumulative <- t(cumsum(phylo_eigperc)) 
phylo_eigperc_cumulative 


#' Create a dataframe of species scores to merge with other data
#' 
#' Extract species names
phylo_names <- (pvr1@phylo$tip.label)
#' Extract species scores (eigenvectors)
phylo_scores <- as.data.frame(pvr1@Eigen$vectors)
#' Add a column for species' names
phylo_scores$species_name <- phylo_names 
#' Select only the first 4 components
phylo_scores_final <- phylo_scores %>% select(species_name, c1:c4)
#' Make sure format of species name column is 'character'
phylo_scores_final$species_name <- as.character(phylo_scores_final$species_name)


#' Merge scores with full trait dataframe
species_dat_merge <- merge(species_dat, phylo_scores_final, by='species_name')



#' *Save phylogenetic metrics*
#' 
str(species_dat_merge)
phylo_metrics <- species_dat_merge %>% select(species, c1, c2) %>%
  rename("phylo_pcoa1" = "c1","phylo_pcoa2" = "c2")



