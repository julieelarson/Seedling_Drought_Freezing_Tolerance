# Seedling_Drought_Freezing_Tolerance
This repository contains code and data to conduct analyses and create figures included  in the following manuscript in preparation: Larson, J.E., B.J. Butterfield, S.M. Munson, D. Neuhaus, S.M. Copeland. (In Prep) Functional traits can explain coordinated seedling responses to drought and freezing stress

This study analyzes how traits, phylogenetic relatendness, and climate niche explain variation in seedling drought and freezing tolerance. 

The respostiory contains R scripts that sepearately generate several species-level metrics (e.g., climate niche, phyogenetic relatedness scores, drought and freezing tolerance) as well as a cross-data analysis that unites these species'-level metrics to generate the analysese and figures in the main text.

*Source of response variables:* 
1) Lethal freezing temperatures: Generated in freezing_tolerance.R

   -Input: Freeze_Temps.csv

   -Input: Freezing_Mortality.csv

   -Output: freeze_metric.csv
   
2) Drought tolerance: Generated in drought.tolerance.R
    -Input: Drought_Mortality.csv

   -Input:Drought_Soil_Moisture_At_Death.csv

   -Input:Drought_WaterPots.csv

   -Output: drought_metrics.csv


*Source of explanatory variables for analysis:*

1) Seed and seedling traits: Linked from published repository (https://doi.org/10.5061/dryad.b5mkkwhpt)

2) Phylogenetic PCoA scores: Generated in phylogenetic_metrics&analysis.R

   -Note that this script also includes all standalone phylogenetic analyses (e.g., phylogenetic signal estimation)

   -Input: drought_metrics.csv

   -Input: freeze_metric.csv

   -Input: Archeived seed and seedling trait data

   -Output: phylo_metrics.csv

3) Climate niche metrics: Generated in climate_niche.R

   -Input: species_list.csv

   -Output: climate_niche_metrics.csv


