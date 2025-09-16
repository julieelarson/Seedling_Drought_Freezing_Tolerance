#' ---
#' title: "Generating climate niche metrics"
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


#' **This code: Generating climate niche metrics**
#' 
#' In this code, we extract species' occurrence records from the GBIF(Global Biodiversity 
#' Information Facility), then cross these with gridded climate data separately sourced 
#' from WorldClim v2.0 (1-km resolution climate norms for the period 1980-2010; Fick & Hijmans, 2017).
#' 
#' For each species and climate metric, we then 5th, 50th, and 95th percentile from
#' the total distribution of values for a species' range.
#' 
#' For analyses, we ultimately use two of these metrics, which were selected a priori
#' based on the analytical goals of explaining seedling drought and freezing tolerance:
#' -- the 5th percentile of MAT (i.e. species' coldest occurrences, MAT_05)
#' -- the 5th percentile of MAP (i.e. species' driest occurrences, MAP_05)
#' 
#' 
#' Reference:
#' Fick, S. E., and R. J. Hijmans. 2017. "WorldClim 2: new 1-km spatial resolution climate 
#' surfaces for global land areas." International Journal of Climatology 37 (12): 4302-4315. 
#' https://doi.org/10.1002/joc.5086.



########## load packages and set arguments ##############
library(dismo)
library(terra)#NOTE: this code was run during the transition from {raster} to {terra}. The raster functions were still supported by terra, but are now deprecated

args=commandArgs(TRUE)#For running a job array in a HPC environment

spNum=args[1]#Running all of the code below for each species


########## selecting species and extent method #################################################
spList=read.csv("species_list.csv",na.strings='',header=T)#A vector of latin binomials, e.g. 'Agoseris glauca'
sp.i=as.vector(spList[spNum,1])#For running in a HPC environment, or parallel computing. Alternative is to loop this and run for each species

#using dismo::gbif to extract species occurrences
if (spNum == 2){#specifying sub-species for two species widespread species
  occ.i <- gbif('Achillea','millefolium occidentalis')
}else{
  if (spNum == 3){
    occ.i <- gbif('Artemisia','tridentata wyomingensis')
  }else{
  occ.i=gbif(strsplit(sp.i,' ')[[1]][1],strsplit(sp.i,' ')[[1]][2])
  }
}

occ.i=cbind(occ.i$lon,occ.i$lat,rep(1,nrow(occ.i)))#pulling coordinates
occ.i=occ.i[complete.cases(occ.i),]#removing any occurrences without complete lat/lon records

### reading in climate data
envStack=brick(".../climDat.grd")
envUnstack=unstack(envStack)
envStack=stack(envUnstack[[1]],envUnstack[[5]])
rm(envUnstack)

### removing duplicates and occurrences without climate data
occ_filtered.i=rasterize(occ.i[,1:2],envStack,field=occ.i[,3],fun=function(x,...)max(x))

respXY=cbind(xyFromCell(occ_filtered.i,cell=which(occ_filtered.i[]==1)),c(rep(1)))
climDat=extract(envStack,respXY[,1:2])
respXY=respXY[complete.cases(respXY),]
climDat=climDat[complete.cases(climDat),]

# determining climate background
# this helps to eliminate occurrence records with extreme climate values that would unduly influence climate niche estimates
# this excludes values that are >2sd
SD=2
   
climLim=vector()
for (i in 1:ncol(climDat)){
  lim.i=vector()
  for (j in 1:999){
    dist.i=rnorm(3*nrow(climDat),mean=mean(climDat[,i]),sd=sd(climDat[,i])*SD)
    lim.i=cbind(lim.i,c(min(dist.i),max(dist.i)))
  }
  climLim=cbind(climLim,rowMeans(lim.i))
}
colnames(climLim)=colnames(climDat)
rownames(climLim)=c("min","max")

is.between=function(x,a,b){
  (x-a)*(b-x)>0
}

climLimRast=is.between(envStack[[1]],climLim[1,1],climLim[2,1])
for (i in 2:ncol(climLim)){
  climLimRast=climLimRast*is.between(envStack[[i]],climLim[1,i],climLim[2,i])
}

#removing climatic outlier presences
presRast=rasterize(respXY[,1:2],climLimRast,field=respXY[,3],fun=function(x,...)max(x))
presRast=mask(presRast,climLimRast,maskvalue=0)
presRast=mask(presRast,climLimRast)
respXY=xyFromCell(presRast,cell=which(!is.na(presRast[])))
respXY=cbind(respXY,rep(1,nrow(respXY)))

write.csv(respXY,"respXY.csv")

#extracting climate niche values
pct05=function(x){
  x=x[!is.na(x)]
  x=sort(x)
  return(x[floor(length(x)*.05)])
}

pct95=function(x){
  x=x[!is.na(x)]
  x=sort(x)
  return(x[ceiling(length(x)*.95)])
}

climDat=extract(envStack,respXY[,1:2])

#calculating climatic niche values for each species
sppMed=apply(climDat,2,median,na.rm=T)
spp05=apply(climDat,2,pct05)
spp95=apply(climDat,2,pct95)