# Seedling_Drought_Freezing_Tolerance
TThis repository contains code and data to conduct analyses and create figures included in the following manuscript:

Functional traits can explain coordinated seedling responses to drought and freezing stress.

Author Names: Julie E. Larson*1,2, Bradley J. Butterfield3, Seth M. Munson4, Dylan F. Neuhaus1,5, Stella M. Copeland1

Author Institutions and Addresses: 
1USDA-Agricultural Research Service, Eastern Oregon Agricultural Research Center, Burns, OR, 97720, USA 
2Present address: School of Environmental and Forest Sciences, University of Washington, Seattle, WA, 98195, USA 
3Department of Biological Sciences, Northern Arizona University, Flagstaff, AZ, 86011-5640, USA 
4U.S. Geological Survey, Southwest Biological Science Center, Flagstaff, AZ 86001, USA 
5Present address: University of Central Arkansas, Conway, AR, 72035, USA 
*Corresponding Author: [larsonju@uw.edu](mailto:larsonju@uw.edu)

This study analyzes how traits, phylogenetic relatedness, and climate niche explain variation in seedling drought and freezing tolerance.

The repository contains R scripts that separately generate several species-level metrics (e.g., climate niche, phylogenetic relatedness scores, drought and freezing tolerance) as well as a cross-data analysis that unites these species-level metrics to generate the analyses and figures in the main text.


*Source of response variables for analysis:* 
1) Lethal freezing temperatures: Generated in freezing_tolerance.R

   -Input: Freeze_Temps.csv

   -Input: Freezing_Mortality.csv

   -Output: freeze_metric.csv
   
2) Drought tolerance: Generated in drought.tolerance.R
   -Input: Drought_Mortality.csv

   -Input: Drought_Soil_Moisture_At_Death.csv

   -Input: Drought_WaterPots.csv

   -Output: drought_metrics.csv


*Source of explanatory variables for analysis:*

1) Seed and seedling traits: Linked from published repository (https://doi.org/10.5061/dryad.b5mkkwhpt)

2) Phylogenetic PCoA scores: Generated in phylogenetic_metrics&analysis.R

   -Note that this script also includes all standalone phylogenetic analyses (e.g., phylogenetic signal estimation)

   -Input: drought_metrics.csv

   -Input: freeze_metric.csv

   -Input: Archived seed and seedling trait data

   -Output: phylo_metrics.csv

3) Climate niche metrics: Generated in climate_niche.R

   -Input: species_list.csv

   -Output: climate_niche_metrics.csv


