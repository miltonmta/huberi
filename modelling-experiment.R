# --- Packages                                      ----
# install.packages(c("tidyverse", "raster", "rgdal", "abind", "spThin", "vegan", "maps", "mask","pcych", "kernlab", "dismo", "rJava")
require(tidyverse)
require(raster)
require(rgdal)
require(abind)
require(spThin)
require(vegan)
require(maps)
require(mask)
require(kernlab)
require(dismo)
require(rJava)

# --- List of improvements to the scritp            ----

# 1. Implement validation by the checkerboards method.
# 2. Implement occurrence filtering at the ambiental space and compare with spThin (geographical space)
# 3. Transform maps in frequencies instead of suitabilities.
# 4. Impelement multi cores for running several models simultaneously.
# 5. Reduce code by implementing subfunctions, lopps.
# 6. Rewrite the code using tidy

# --- Maxent/rJava troubleshoot                     ----

# MaxEnt is available as a standalone Java program. Dismo has a function 'maxent' that communicates with this program. To use it you must first download the program from:
# browseURL("http://www.cs.princeton.edu/~schapire/maxent/")
# Put the file 'maxent.jar' in the 'java' folder of the 'dismo' package. That is the folder returned by system.file("java", package="dismo"). You need MaxEnt version 3.3.3b or higher.

# While loagind the package rJava, if you have the error image not found, follow the steps:
# 1. Check if JDK is installed in you machine. If not:
# browseURL("http://www.oracle.com/technetwork/java/javase/downloads/jdk8-downloads-2133151.html")

# 2. Checking if the jar file is present. 
# jar <- paste(system.file(package="dismo"), "/java/maxent.jar", sep='')
# file.exists(jar)

# 3. Use dyn.load for ataching the libjvm.dylib file before running the `rJava` Package.
# dyn.load('/Library/Java/JavaVirtualMachines/jdk1.8.0_131.jdk/Contents/Home/jre/lib/server/libjvm.dylib')

# 4. If you are in MacOS, a permanent solution is to link libjvm.dylib to /usr/local.lib at terminal:
# sudo ln -f -s $(/usr/libexec/java_home)/jre/lib/server/libjvm.dylib /usr/local/lib
# browseURL ("https://stackoverflow.com/questions/30738974/rjava-load-error-in-rstudio-r-after-upgrading-to-osx-yosemite")

# --- General Info                                  ----

# This script has an index table. If you are in RStudio go to Code > Show Document Outline (shift + command / clrt + o)

# The directories with the data utilized and the ones outputted here can be downloaded from the following OneDrive repositorium:
browseURL("https://1drv.ms/f/s!ApJZaitgpPr7gZtfS9n9mU9DDzXQMg")

# The models projected for 2070 (average for 2061-2080) were obtain at "worldclim.com" by the spatial resolution of 2.5min (0.04º or ≈ 4.4km). We selected the variables that appear simultaneously in all the representative concentration pathways scnarios (RCP26, RCP45, RCP60, RCP80). The codes of the 11 GCMs utilized are: bc, cc, gs, hd, he, ip, mi, mr, mc, mg, no.
# browseURL("http://www.worldclim.org/cmip5_2.5m")


# ***************************************************************************************
## 01. Read aogcms models                           ----

### 1.1 current ----
# Worldclim variables for current conditions (~1960-1990), by the spatial resolution of 2.5min. 
# browseURL("http://www.worldclim.org/current")

###.............. processing the current climatic model
read_current <- function (dir)                                
{
  model_raw <- stack(list.files(dir,  pattern = ".bil$", full.names = TRUE))
  e <- extent(-122, -18, -56, 14) 
  model_e <- crop(model_raw, e)                              
  val <- getValues(model_e)                                  
  coord <- xyFromCell(model_e, 1:ncell(model_e))             
  model <- cbind(coord, val)                                 
  model <- na.omit(model)                                    
  return(model)
}

###.............. Creating a matrix
current <- read_current(dir = "./data/climatic_vars/current")

###.............. Creating a RasterBrick
current_spatial <- rasterFromXYZ(current) 


### 1.2 rcp  ----
# Worldclim variables for the RCPs, projected to 2070(average for 2061-2080), by the spatial resolution of 2.5min. 
# browseURL("http://www.worldclim.org/cmip5_2.5m")

###.............. wolrdclim GCM codes

# BCC-CSM1-1	      BC
# CCSM4	            CC
# GISS-E2-R	        GS
# HadGEM2-AO	      HD
# HadGEM2-ES        HE	
# IPSL-CM5A-LR	    IP
# MIROC5            MC
# MRI-CGCM3	        MG
# MIROC-ESM-CHEM    MI
# MIROC-ESM    	    MR
# NorESM1-M	        NO

###.............. selected AOGCMs
# based on the cluster analysis we'll import only the selected variables at each RCP scenario.

# RCP 26: CCSM4(CC), HADGEM2-AO(HD),   IPSL-CMSA-LR(IP), MIROC-ESM(MR) 
# RCP 45: CCSM4(CC), HADGEM2-AO(HD),   IPSL-CMSA-LR(IP), MIROC-ESM(MR) 
# RCP 60: CCSM4(CC), IPSL-CMSA-LR(IP), MRI-CGCM3(MG),    MIROC-ESM(MR)
# RCP 85: CCSM4(CC), IPSL-CMSA-LR(IP), MIROC5(MC),       MIROC-ESM(MR)

# for maintanining comparability... 

# RCP 26: CCSM4(CC),                   IPSL-CMSA-LR(IP), MIROC-ESM(MR) 
# RCP 45: CCSM4(CC),                   IPSL-CMSA-LR(IP), MIROC-ESM(MR) 
# RCP 60: CCSM4(CC), IPSL-CMSA-LR(IP),                   MIROC-ESM(MR)
# RCP 85: CCSM4(CC), IPSL-CMSA-LR(IP),                   MIROC-ESM(MR)

# naming must be in the exact order of the origin directory.
model_names <- c("CCSM4", "IPSL-CM5A-LR", "MIROC-ESM")

###.............. processing rcps
read_rcp <- function (x)
{
  directories <- list.dirs( x, full.names = TRUE)[-1]
  e <- extent(-122, -18, -56, 14)
  models_list <- list()
  rcp <- NULL
  for (i in 1:length(directories))
  {
    models_raw       <- stack(list.files(directories[i],pattern = ".tif$", full.names = TRUE))
    models_e         <- crop( models_raw , e )
    val              <- values (models_e)
    coord            <- xyFromCell(models_e, 1:ncell(models_e))
    models           <- cbind(coord, val)
    models           <- na.omit(models)
    models_list[[i]] <- models_e
    rcp              <- abind (rcp, models, along = 3)
  }
  return(list("array" = rcp, "rasters" = models_list))
}

rcp26_list <- read_rcp( x = "./data/climatic_vars/selected/26bi70/")
rcp45_list <- read_rcp( x = "./data/climatic_vars/selected/45bi70/")
rcp60_list <- read_rcp( x = "./data/climatic_vars/selected/60bi70/")
rcp85_list <- read_rcp( x = "./data/climatic_vars/selected/85bi70/")

###.............. object type: array
rcp26 <- rcp26_list[["array"]]
rcp45 <- rcp45_list[["array"]]
rcp60 <- rcp60_list[["array"]]
rcp85 <- rcp85_list[["array"]]

###.............. object type: RasterStack
rcp26_spatial <- rcp26_list[["rasters"]]
rcp26_spatial <- stack(rcp26_spatial[[1]], rcp26_spatial[[2]], rcp26_spatial[[3]])

rcp45_spatial <- rcp45_list[["rasters"]]
rcp45_spatial <- stack(rcp45_spatial[[1]], rcp45_spatial[[2]], rcp45_spatial[[3]])

rcp60_spatial <- rcp60_list[["rasters"]]
rcp60_spatial <- stack(rcp60_spatial[[1]], rcp60_spatial[[2]], rcp60_spatial[[3]])

rcp85_spatial <- rcp85_list[["rasters"]]
rcp85_spatial <- stack(rcp85_spatial[[1]], rcp85_spatial[[2]], rcp85_spatial[[3]])


rm(rcp26_list, rcp45_list, rcp60_list, rcp85_list)

# ***************************************************************************************
## 02. Variable selection                           ----
###.............. by exploratory factor analysis

fa.parallel(current[ , -c(1:2)], fa = 'fa') #scree plot
current_fa <- fa(current[ , -c(1:2)], nfactors = 5, rotate = 'varimax')
loadings <- loadings(current_fa)
write.table(loadings, "./data/climatic_vars/selected/varimax_loadings.txt")

###.............. Selected variables 
# bio02, bio03, bio10, bio14, bio16.

###.............. bioclimatic variables descriptions

#BIO01 = Annual Mean Temperature
#BIO02 = Mean Diurnal Range (Mean of monthly (max temp - min temp))
#BIO03 = Isothermality (BIO2/BIO7) (* 100)
#BIO04 = Temperature Seasonality (standard deviation *100)
#BIO05 = Max Temperature of Warmest Month
#BIO06 = Min Temperature of Coldest Month
#BIO07 = Temperature Annual Range (BIO5-BIO6)
#BIO08 = Mean Temperature of Wettest Quarter
#BIO09 = Mean Temperature of Driest Quarter
#BIO10 = Mean Temperature of Warmest Quarter
#BIO11 = Mean Temperature of Coldest Quarter
#BIO12 = Annual Precipitation
#BIO13 = Precipitation of Wettest Month
#BIO14 = Precipitation of Driest Month
#BIO15 = Precipitation Seasonality (Coefficient of Variation)
#BIO16 = Precipitation of Wettest Quarter
#BIO17 = Precipitation of Driest Quarter
#BIO18 = Precipitation of Warmest Quarter
#BIO19 = Precipitation of Coldest Quarter

# ***************************************************************************************
## 03. Saving selected variables                    ----
### 3.1a. current - saving as table ----

write.table(current[,c("x", "y", "bio02", "bio03", "bio10", "bio14", "bio16")], "./data/climatic_vars/selected/current/current-select.txt", row.names = F, sep = " ")  

### 3.1b. current - saving as raster ----

variables <- as.factor(c("bio02", "bio03", "bio10", "bio14", "bio16"))
for (i in 1:length(variables))
{
  writeRaster (current_spatial[[i]], filename = paste0("./data/climatic_vars/selected/current/current-", variables[i], ".grd"), format = "raster")
}
rm(variables)

current_select <- stack(list.files("./data/climatic_vars/selected/current/",  pattern = ".grd$", full.names = TRUE))

#### RCPs

### 3.2a. rcp - saving table ----

write.table(rcp26 [ ,c("x", "y", "bio02", "bio03", "bio10", "bio14", "bio16" ), ], "./data/climatic_vars/selected/rcp26/rcp26-select.txt", row.names = F, sep = "	")
write.table(rcp45 [ ,c("x", "y", "bio02", "bio03", "bio10", "bio14", "bio16" ), ], "./data/climatic_vars/selected/rcp45/rcp45-select.txt", row.names = F, sep = "	")
write.table(rcp60 [ ,c("x", "y", "bio02", "bio03", "bio10", "bio14", "bio16" ), ], "./data/climatic_vars/selected/rcp60/rcp60-select.txt", row.names = F, sep = "	")
write.table(rcp85 [ ,c("x", "y", "bio02", "bio03", "bio10", "bio14", "bio16" ), ], "./data/climatic_vars/selected/rcp85/rcp85-select.txt", row.names = F, sep = "	")


### 3.2b. rcp - saving as raster ----

###.............. Creating rasters with the selected variables from each aogcm

variables <- as.factor(c("bio02.1", "bio03.1", "bio10.1", "bio14.1", "bio16.1", "bio02.2", "bio03.2", "bio10.2", "bio14.2", "bio16.2", "bio02.3", "bio03.3", "bio10.3", "bio14.3", "bio16.3"))
## RCP26
for (i in 1:length(variables))
{
  writeRaster (rcp26_spatial[[i]], filename = paste0("./data/climatic_vars/selected/rcp26/rcp26-", variables[i], ".grd"), format = "raster")
}

## RCP45
for (i in 1:length(variables))
{
  writeRaster (rcp45_spatial[[i]], filename = paste0("./data/climatic_vars/selected/rcp45/rcp45-", variables[i], ".grd"), format = "raster")
}

## RCP60
for (i in 1:length(variables))
{
  writeRaster (rcp60_spatial[[i]], filename = paste0("./data/climatic_vars/selected/rcp60/rcp60-", variables[i], ".grd"), format = "raster")
}

## RCP85
for (i in 1:length(variables))
{
  writeRaster (rcp85_spatial[[i]], filename = paste0("./data/climatic_vars/selected/rcp85/rcp85-", variables[i], ".grd"), format = "raster")
}

rm(variables)


# ***************************************************************************************
## 04. Occurrences data                             ----

###.............. Reading data
# if you have more than one species, all of them should be in the file raw_data.
# the names of columms in this file must be: "SPEC", "LONG", "LAT".
occur_raw <- read.table("./data/occurrences/raw_data.txt", h = T)
 
###.............. Filtering occurrences in the geographical space
sp <- gsub("C[1-9]","", occur_raw$SPEC)
sp_names <- unique(sp)
occur_thinned <- NULL
for(i in 1:length(sp_names))
{
  sp <- occur_raw[occur_raw[, 1] == sp_names[i], ]
  occur <- thin(sp, thin.par = 20, reps = 100, max.files = 1, locs.thinned.list.return = T)
  summaryThin(occur)
  occur_thinned <- rbind(occur_thinned, occur)
}
write.table(occur_thinned, "./data/occurrences/occur_thinned.txt", sep = ";", row.names = FALSE)

# occur_thinned <- read.csv("./data/occurrences/occur_thinned.csv", sep = ",", h = T)

# ***************************************************************************************
## 05. Extracting variables                         ----

create_var <- function(sp,name)
{
  sp_cell <- cellFromXY(current_select, sp[, -1])
  duplicated(sp_cell)
  sp_cell <- unique(sp_cell)
  sp_coord <- xyFromCell(current_select, sp_cell)
  sp_var <- raster::extract(current_select, sp_cell)
  sp_var <- na.omit(cbind(sp_coord, sp_var))
  write.table(sp_var, file = paste0("./data/occurrences/var_", as.factor(name), ".txt"), row.names = F, sep = ";")
  return(sp_var)
}

# Creating and saving the object "var" for each studied species
for(i in 1:length(sp_names))
{
  var <- create_var(occur_thinned[occur_thinned[, 1] == sp_names[i], ], sp_names[i])
}

# ***************************************************************************************
## 06. Background Sampling                          ----

###.............. Creating and saving the object "back" for each studied species
create_back <- function(sp, name)
{
  coord <- rasterToPoints(current_select)[, 1:2]
  back_id <- sample(1:nrow(coord), nrow(sp))
  back <- extract(current_select, coord[back_id, ])
  back <- cbind(coord [back_id, ], back)
  write.table(back, paste0("./data/occurrences/back_", as.factor(name), ".txt"), row.names = F, sep = ";") 
  return(back)
}

var_files <- list.files("./data/occurrences/", pattern = "var", full.names = TRUE)
for(i in 1:length(var_files))
{
  var_file <- read.table(var_files[i], h = T, sep = ";")
  create_back(var_file, sp_names[i])
}

###.............. Plotting occurrences and background
sp_names
plot(current_select$bio01)
points(occur_thinned[-(occur_thinned[, 1] == "Lithurgus_huberi"), ][, -1], pch = "*", col = "red")
points(occur_thinned[occur_thinned[, 1] == "Lithurgus_huberi", ][,-1], pch = "*", col = "blue")
# points(back[, "x"], back[, "y"], pch = "*", col = 'black')
# points(back[, "x"], back[, "y"], pch = "*", col = 'magenta')

# ***************************************************************************************
## 07. Modelling Predictions                        ----
rm(list = ls())

species_model <- function(occurrence, 
                          background, 
                          biovar_current,
                          biovar_rcp26,
                          biovar_rcp45,
                          biovar_rcp60,
                          biovar_rcp85 ,
                          cross_validation)
{

  ###.............. Reading the selected bioclimatic variables
  current_select <- stack(list.files(biovar_current,  pattern = ".grd$", full.names = TRUE))
  
  ###.............. Reading ocurrence and background
  occur <- read.table(occurrence, sep = ";", h = T)
  back  <- read.table(background, sep = ";", h = T)
  
  ###.............. Creating objects for saving the results
  ## Predictive models
  bioclim_c <- gower_c  <- maxent_c <- SVM_c <- stack()
  bioclim_rcp26 <- gower_rcp26 <- maxent_rcp26 <- SVM_rcp26 <- stack() 
  bioclim_rcp45 <- gower_rcp45 <- maxent_rcp45 <- SVM_rcp45 <- stack()
  bioclim_rcp60 <- gower_rcp60 <- maxent_rcp60 <- SVM_rcp60 <- stack()
  bioclim_rcp85 <- gower_rcp85 <- maxent_rcp85 <- SVM_rcp85 <- stack()
  
  ## Evaluation
  bioclim_e <- gower_e <- maxent_e <- SVM_e <- NULL # TRF
  bioclim_t <- gower_t <- maxent_t <- SVM_t <- NULL # Threshold type comission
  bioclim_d <- gower_d <- maxent_d <- SVM_d <- NULL # Area predicted as presence
  
  ## Partial results
  bioclim_Pout_c <- gower_Pout_c <- maxent_Pout_c <- SVM_Pout_c <- NULL
  bioclim_Pout_rcp26 <- gower_Pout_rcp26 <- maxent_Pout_rcp26 <- SVM_Pout_rcp26 <- NULL
  bioclim_Pout_rcp45 <- gower_Pout_rcp45 <- maxent_Pout_rcp45 <- SVM_Pout_rcp45 <- NULL
  bioclim_Pout_rcp60 <- gower_Pout_rcp60 <- maxent_Pout_rcp60 <- SVM_Pout_rcp60 <- NULL
  bioclim_Pout_rcp85 <- gower_Pout_rcp85 <- maxent_Pout_rcp85 <- SVM_Pout_rcp85 <- NULL
  
  ## Predictions results
  output_current <- output_rcp26 <- output_rcp45 <- output_rcp60 <- output_rcp85 <- NULL
  
  
  for (j in 1:cross_validation)
  {
    ### OPEN "j" ----
    
    ###.............. Creating trainning-testing subsets
    sample_occur <- sample(1:nrow(occur), round(0.75 * nrow(occur), 0))
    trainning <- prepareData(x = current_select, p = occur[sample_occur,  1:2], b = back[sample_occur,  1:2]) 
    testing   <- prepareData(x = current_select, p = occur[-sample_occur, 1:2], b = back[-sample_occur, 1:2])
    
    
    # ***************************************************************************************
    ###.............. Predictive models
    
    ### Bioclim -----------------------------------------------------------------------------
    
    ## Adjusting models
    bioclim_model <- bioclim(trainning[trainning[, "pb"] == 1, -1])
  
    ## Predicting
    bioclim_c <- stack(bioclim_c, predict(object = bioclim_model, x = current_select))
    
    ## Evaluating models
    bioclim_eval <- evaluate(p = testing[testing[, "pb"] == 1, -1], a = testing[testing[, "pb"] == 0, -1], model = bioclim_model)
    
    ## Saving evaluation metrics
    thr <- threshold(bioclim_eval, "no_omission") # the highest threshold at which there is no omission
    TPR <- bioclim_eval@TPR[which( bioclim_eval@t == thr)]
    bioclim_e <- c(bioclim_e, TPR) # TPR - True positive rate 
    bioclim_t <- c(bioclim_t, thr) 
    
    ## Study area predicted as presence
    n_cells <- nrow(na.omit(values(current_select))) 
    pi <- sum(values(bioclim_c >= thr), na.rm = T) / n_cells # predicted area proportion.
    bioclim_d <- TPR * (1 - pi)
    
    
    ### Gower -------------------------------------------------------------------------------
   
    ## Adjusting models
    gower_model <- domain(trainning[trainning[,"pb"] == 1, -1]) 
    
    ## Predicting
    gower_c <- stack(gower_c, predict(object = gower_model, x = current_select))
    
    ## Evaluating models
    gower_eval <- evaluate(p = testing[testing[, "pb"] == 1, -1], a = testing[testing[, "pb"] == 0, -1], model = gower_model)
    
    ## Saving evaluation metrics
    thr <- threshold(gower_eval, "no_omission")
    TPR <- gower_eval@TPR[which( gower_eval@t == thr)]
    gower_e <- c(gower_e, TPR)
    gower_t <- c(gower_t, thr)
    
    ## Study area predicted as presence
    n_cells <- nrow(na.omit(values(current_select))) 
    pi <- sum(values(gower_c >= thr), na.rm = T) / n_cells
    gower_d <- TPR * (1 - pi)
    
    
    ### Maxent -------------------------------------------------------------------------------
    
    ## Adjusting models
    Sys.setenv(NOAWT = TRUE)
    maxent_model <- maxent(x = trainning[, -1], p = trainning[, 1])
    
    ## Predicting
    maxent_c <- stack(maxent_c, predict(object = maxent_model, x = current_select))
    
    ## Evaluating models
    maxent_eval <- evaluate(p = testing[testing[, "pb"] == 1, -1], a = testing[testing[, "pb"] == 0, -1], model = maxent_model, tr = 0.5414856)
    
    ## Saving evaluation metrics
    thr <- threshold(maxent_eval, "no_omission")
    TPR <- maxent_eval@TPR[which( maxent_eval@t == thr)]
    maxent_e <- c(maxent_e, TPR) 
    maxent_t <- c(maxent_t, thr)
    
    ## Study area predicted as presence
    n_cells <- nrow(na.omit(values(current_select))) 
    pi <- sum(values(maxent_c >= thr), na.rm = T) / n_cells
    TPR <- maxent_e
    maxent_d <- TPR * (1 - pi)
    
    
    ### SVM ----------------------------------------------------------------------------------
    
    ## Adjusting models
    SVM_model <- ksvm(pb ~ ., data = trainning)
    
    ## Predicting
    SVM_c <- stack(SVM_c, predict(model = SVM_model, object = current_select)) 
    
    ## Evaluating models
    SVM_eval <- evaluate(p = testing[testing[, "pb"] == 1, -1], a = testing[testing[, "pb"] == 0, -1], model = SVM_model, tr = 0.5566005)
    
    ## Saving evaluation metrics
    thr <- threshold(SVM_eval, "no_omission")
    TPR <- SVM_eval@TPR[which( SVM_eval@t == thr)]
    SVM_e <- c(SVM_e, TPR)
    SVM_t <- c(SVM_t, thr)
    
    ## Study area predicted as presence 
    n_cells <- nrow(na.omit(values(current_select))) 
    pi <- sum(values(SVM_c >= thr), na.rm = T) / n_cells
    TPR <- SVM_e
    SVM_d <- TPR * (1 - pi)
    
    
    # ***************************************************************************************
    
    ###.............. Saving partial outputs for current model
    bioclim_Pout_c <- cbind(bioclim_Pout_c, values(bioclim_c))
    gower_Pout_c   <- cbind(gower_Pout_c,   values(gower_c))
    maxent_Pout_c  <- cbind(maxent_Pout_c,  values(maxent_c))
    SVM_Pout_c     <- cbind(SVM_Pout_c,     values(SVM_c))
    
    
    ###.............. Making predictions for the RCPs
    
    AOGCMs <- c("CCSM4", "IPSL-CM5A-LR", "MIROC-ESM")
    for (i in AOGCMs)
    {
      # OPEN "i" ----
      ### Reading the selected aogcm models
      
      rcp26_select  <- stack(list.files(biovar_rcp26, pattern = ".grd$", full.names = TRUE))
      rcp45_select  <- stack(list.files(biovar_rcp45, pattern = ".grd$", full.names = TRUE))
      rcp60_select  <- stack(list.files(biovar_rcp60, pattern = ".grd$", full.names = TRUE))
      rcp85_select  <- stack(list.files(biovar_rcp85, pattern = ".grd$", full.names = TRUE))
      
      ### Predicting
      bioclim_rcp26 <- stack(bioclim_rcp26, predict(object = bioclim_model, x = rcp26_select))
      bioclim_rcp45 <- stack(bioclim_rcp45, predict(object = bioclim_model, x = rcp45_select))
      bioclim_rcp60 <- stack(bioclim_rcp60, predict(object = bioclim_model, x = rcp60_select))
      bioclim_rcp85 <- stack(bioclim_rcp85, predict(object = bioclim_model, x = rcp85_select))
      
      gower_rcp26   <- stack(gower_rcp26,   predict(object = gower_model, x = rcp26_select))
      gower_rcp45   <- stack(gower_rcp45,   predict(object = gower_model, x = rcp45_select))
      gower_rcp60   <- stack(gower_rcp60,   predict(object = gower_model, x = rcp60_select))
      gower_rcp85   <- stack(gower_rcp85,   predict(object = gower_model, x = rcp85_select))
      
      maxent_rcp26  <- stack(maxent_rcp26,  predict(object = maxent_model, x = rcp26_select))
      maxent_rcp45  <- stack(maxent_rcp45,  predict(object = maxent_model, x = rcp45_select))
      maxent_rcp60  <- stack(maxent_rcp60,  predict(object = maxent_model, x = rcp60_select))
      maxent_rcp85  <- stack(maxent_rcp85,  predict(object = maxent_model, x = rcp85_select))
      
      SVM_rcp26     <- stack(SVM_rcp26,     predict(model = SVM_model, object = rcp26_select))
      SVM_rcp45     <- stack(SVM_rcp45,     predict(model = SVM_model, object = rcp45_select))
      SVM_rcp60     <- stack(SVM_rcp60,     predict(model = SVM_model, object = rcp60_select))
      SVM_rcp85     <- stack(SVM_rcp85,     predict(model = SVM_model, object = rcp85_select))
      
      
      ### Saving partial outputs for the RCPs
      
      bioclim_Pout_rcp26 <- cbind(bioclim_Pout_rcp26, values(bioclim_rcp26))
      bioclim_Pout_rcp45 <- cbind(bioclim_Pout_rcp45, values(bioclim_rcp45))
      bioclim_Pout_rcp60 <- cbind(bioclim_Pout_rcp60, values(bioclim_rcp60))
      bioclim_Pout_rcp85 <- cbind(bioclim_Pout_rcp85, values(bioclim_rcp85))
      
      gower_Pout_rcp26   <- cbind(gower_Pout_rcp26,   values(gower_rcp26))
      gower_Pout_rcp45   <- cbind(gower_Pout_rcp45,   values(gower_rcp45))
      gower_Pout_rcp60   <- cbind(gower_Pout_rcp60,   values(gower_rcp60))
      gower_Pout_rcp85   <- cbind(gower_Pout_rcp85,   values(gower_rcp85))
      
      maxent_Pout_rcp26  <- cbind(maxent_Pout_rcp26,  values(maxent_rcp26))
      maxent_Pout_rcp45  <- cbind(maxent_Pout_rcp45,  values(maxent_rcp45))
      maxent_Pout_rcp60  <- cbind(maxent_Pout_rcp60,  values(maxent_rcp60))
      maxent_Pout_rcp85  <- cbind(maxent_Pout_rcp85,  values(maxent_rcp85))
      
      SVM_Pout_rcp45     <- cbind(SVM_Pout_rcp45,     values(SVM_rcp45))
      SVM_Pout_rcp26     <- cbind(SVM_Pout_rcp26,     values(SVM_rcp26))
      SVM_Pout_rcp60     <- cbind(SVM_Pout_rcp60,     values(SVM_rcp60))
      SVM_Pout_rcp85     <- cbind(SVM_Pout_rcp85,     values(SVM_rcp85))
     
      # CLOSE "i" ----
      # AOGCMs
    }
    
    # CLOSE "j"  ---- 
    # cross-validation
  }   
  
  ###.............. Calculating mean of partial models outputs
  bioclim_Pout_c_mean     <- apply(bioclim_Pout_c, 1,     mean)
  bioclim_Pout_rcp26_mean <- apply(bioclim_Pout_rcp26, 1, mean)
  bioclim_Pout_rcp45_mean <- apply(bioclim_Pout_rcp45, 1, mean)
  bioclim_Pout_rcp60_mean <- apply(bioclim_Pout_rcp60, 1, mean)
  bioclim_Pout_rcp85_mean <- apply(bioclim_Pout_rcp85, 1, mean)
  
  gower_Pout_c_mean       <- apply(gower_Pout_c, 1,     mean)
  gower_Pout_rcp26_mean   <- apply(gower_Pout_rcp26, 1, mean)
  gower_Pout_rcp45_mean   <- apply(gower_Pout_rcp45, 1, mean)
  gower_Pout_rcp60_mean   <- apply(gower_Pout_rcp60, 1, mean)
  gower_Pout_rcp85_mean   <- apply(gower_Pout_rcp85, 1, mean)
  
  maxent_Pout_c_mean      <- apply(maxent_Pout_c, 1,     mean)
  maxent_Pout_rcp26_mean  <- apply(maxent_Pout_rcp26, 1, mean)
  maxent_Pout_rcp45_mean  <- apply(maxent_Pout_rcp45, 1, mean)
  maxent_Pout_rcp60_mean  <- apply(maxent_Pout_rcp60, 1, mean)
  maxent_Pout_rcp85_mean  <- apply(maxent_Pout_rcp85, 1, mean)
  
  SVM_Pout_c_mean         <- apply(SVM_Pout_c, 1,     mean)
  SVM_Pout_rcp26_mean     <- apply(SVM_Pout_rcp26, 1, mean)
  SVM_Pout_rcp45_mean     <- apply(SVM_Pout_rcp45, 1, mean)
  SVM_Pout_rcp60_mean     <- apply(SVM_Pout_rcp60, 1, mean)
  SVM_Pout_rcp85_mean     <- apply(SVM_Pout_rcp85, 1, mean)
  
  
  ###.............. Saving data into the output objects
  output_current <- cbind(output_current, Bioclim = bioclim_Pout_c_mean, Gower = gower_Pout_c_mean, Maxent = maxent_Pout_c_mean, SVM = SVM_Pout_c_mean )
  
  output_rcp26   <- cbind(output_rcp26, Bioclim = bioclim_Pout_rcp26_mean, Gower = gower_Pout_rcp26_mean, Maxent = maxent_Pout_rcp26_mean, SVM = SVM_Pout_rcp26_mean )
  
  output_rcp45   <- cbind(output_rcp45, Bioclim = bioclim_Pout_rcp45_mean, Gower = gower_Pout_rcp45_mean, Maxent = maxent_Pout_rcp45_mean, SVM = SVM_Pout_rcp45_mean )
  
  output_rcp60   <- cbind(output_rcp60, Bioclim = bioclim_Pout_rcp60_mean, Gower = gower_Pout_rcp60_mean, Maxent = maxent_Pout_rcp60_mean, SVM = SVM_Pout_rcp60_mean )
  
  output_rcp85   <- cbind(output_rcp85, Bioclim = bioclim_Pout_rcp85_mean, Gower = gower_Pout_rcp85_mean, Maxent = maxent_Pout_rcp85_mean, SVM = SVM_Pout_rcp85_mean )
  
  
  ###.............. Inserting coords to outputs
  coords <- xyFromCell(current_select, 1:ncell(current_select))
  
  output_current <- cbind(coords,output_current)
  output_rcp26   <- cbind(coords,output_rcp26)
  output_rcp45   <- cbind(coords,output_rcp45)
  output_rcp60   <- cbind(coords,output_rcp60)
  output_rcp85   <- cbind(coords,output_rcp85)
  
  
  ###.............. Excluding NAs from outputs 
  output_current <- na.omit(output_current)
  output_rcp26   <- na.omit(output_rcp26)
  output_rcp45   <- na.omit(output_rcp45)
  output_rcp60   <- na.omit(output_rcp60)
  output_rcp85   <- na.omit(output_rcp85)
  
  
  ###.............. Evaluation data
  models_e <- data.frame("bioclim" = bioclim_e , "gower" = gower_e, "maxent" = maxent_e, "SVM" = SVM_e)

  models_t <- data.frame("bioclim" = bioclim_t , "gower" = gower_t, "maxent" = maxent_t, "SVM" = SVM_t)

  models_d <- data.frame("bioclim" = bioclim_d , "gower" = gower_d, "maxent" = maxent_d, "SVM" = SVM_d)
  
  
  
  #\o/\o/\o/\o/\o/\o/\o/\o/\o/\o/\o/\o/\o/\o/\o/\o/\o/\o/\o/\o/\o/\o/\o/\o/\o/\o/\o/\o/\o/
  return(list(c("output_current" = output_current, 
                "output_rcp26"   = output_rcp26, 
                "output_rcp45"   = output_rcp45, 
                "output_rcp60"   = output_rcp60, 
                "output_rcp85"   = output_rcp85, 
                "TPR"            = models_e, 
                "Threshold"      = models_t, 
                "Predicted_area" = models_d)))
  
} # CLOSE "species_model"----

# ***************************************************************************************
## 08. Running our model                            ----

sp_names
for (i in 1:length(sp_names))
{
  ###.............. Running the modedelling experiment
  result <- species_model(occurrence       = paste("./data/occurrences/var_", sp_names[i], ".txt"),
                          background       = paste("./data/occurrences/back_", sp_names[i],".txt"),
                          biovar_current   = "./data/climatic_vars/selected/current/",
                          biovar_rcp26     = "./data/climatic_vars/selected/rcp26/",
                          biovar_rcp45     = "./data/climatic_vars/selected/rcp45/",
                          biovar_rcp60     = "./data/climatic_vars/selected/rcp60/",
                          biovar_rcp85     = "./data/climatic_vars/selected/rcp85/",
                          cross_validation = 20)
  
  ###.............. Saving predictions
  writeRaster(result[["output_current"]], paste0("./data/outputs/", sp_names[i],"_current.bil"), format = "EHdr")
  writeRaster(result[["output_rcp26"]],   paste0("./data/outputs/", sp_names[i],"_rcp26.bil"), format = "EHdr")
  writeRaster(result[["output_rcp45"]],   paste0("./data/outputs/", sp_names[i],"_rcp45.bil"), format = "EHdr")
  writeRaster(result[["output_rcp60"]],   paste0("./data/outputs/", sp_names[i],"_rcp60.bil"), format = "EHdr")
  writeRaster(result[["output_rcp85"]],   paste0("./data/outputs/", sp_names[i],"_rcp85.bil"), format = "EHdr")
  
  ###.............. Saving evaluation data
  write.table(result[["TPR"]], paste0("./data/outputs/", sp_names[i], "_TPR.txt"), sep = "\t", row.names = F)
  
  write.table(result[["Threshold"]], paste0("./data/outputs/", sp_names[i], "_t.txt"), sep = "\t", row.names = F)
  
  write.table(result[["Predicted_area"]], paste0("./data/outputs/", sp_names[i], "_d.txt"), sep = "\t", row.names = F)
}

# ***************************************************************************************
## 09. Selecting models (TRP)                       ----

  
###.............. asarifolia


###.............. bahiensis


###.............. cairica


###.............. indica


###.............. nil


###.............. purpurea


###.............. huberi
TPR_huberi <- read.table("./data/results/huberi_TPR.txt", h = T)

## Selecting models with TPR values ≥ 0.7
huberi_c     <- raster("./data/results/huberi_current.bil")[[which(TPR_huberi[,] >= 0.7)]]
huberi_rcp26 <- raster("./data/results/huberi_rcp26.bil")[[which(TPR_huberi[,] >= 0.7)]]
huberi_rcp45 <- raster("./data/results/huberi_rcp45.bil")[[which(TPR_huberi[,] >= 0.7)]]
huberi_rcp60 <- raster("./data/results/huberi_rcp60.bil")[[which(TPR_huberi[,] >= 0.7)]]
huberi_rcp85 <- raster("./data/results/huberi_rcp85.bil")[[which(TPR_huberi[,] >= 0.7)]]

all_huberi <- stack(huberi_c, huberi_rcp26, huberi_rcp45, huberi_rcp60, huberi_rcp85)
TPR_huberi <- TPR_h[which(TPR_h)]

##  Creating ANOVA Factors to be used at Uncertainty Evaluation
# huberi_metodo <- rep("bioclim", nlayers(all_huberi)) #???
# huberi_tempo  <- rep(c("pres","fut"), each = nlayers(all_huberi))


### aegyptia



# ***************************************************************************************
## 10. Standardize suitabilities (suit)             ----
# If having error messages, try "range" instead of "standardize"
## Huberi

bioclim_c_h_val <- values(bioclim_c_h)
bioclim_rcp26_h_val <- values(bioclim_rcp26_h)
bioclim_rcp45_h_val <- values(bioclim_rcp45_h)
bioclim_rcp60_h_val <- values(bioclim_rcp60_h)
bioclim_rcp85_h_val <- values(bioclim_rcp85_h)

bioclim_h_val <- rbind(bioclim_c_h_val, bioclim_rcp26_h_val, bioclim_rcp45_h_val, bioclim_rcp60_h_val, bioclim_rcp85_h_val)
bioclim_h_stand <- decostand(bioclim_h_val, "standardize", 2)


gower_c_h_val <- values(gower_c_h)
gower_rcp26_h_val <- values(gower_rcp26_h)
gower_rcp45_h_val <- values(gower_rcp45_h)
gower_rcp60_h_val <- values(gower_rcp60_h)
gower_rcp85_h_val <- values(gower_rcp85_h)

gower_h_val <- rbind(gower_c_h_val, gower_rcp26_h_val, gower_rcp45_h_val, gower_rcp60_h_val, gower_rcp85_h_val)
gower_h_stand <- decostand(gower_h_val, "standardize", 2)


maha_c_h_val <- values(maha_c_h)
maha_rcp26_h_val <- values(maha_rcp26_h)
maha_rcp45_h_val <- values(maha_rcp45_h)
maha_rcp60_h_val <- values(maha_rcp60_h)
maha_rcp85_h_val <- values(maha_rcp85_h)

maha_h_val <- rbind(maha_c_h_val, maha_rcp26_h_val, maha_rcp45_h_val, maha_rcp60_h_val, maha_rcp85_h_val)
maha_h_stand <- decostand(maha_h_val, "standardize", 2)


maxent_c_h_val <- values(maxent_c_h)
maxent_rcp26_h_val <- values(maxent_rcp26_h)
maxent_rcp45_h_val <- values(maxent_rcp45_h)
maxent_rcp60_h_val <- values(maxent_rcp60_h)
maxent_rcp85_h_val <- values(maxent_rcp85_h)

maxent_h_val <- rbind(maxent_c_h_val, maxent_rcp26_h_val, maxent_rcp45_h_val, maxent_rcp60_h_val, maxent_rcp85_h_val)
maxent_h_stand <- decostand(maxent_h_val, "standardize", 2)


SVM_c_h_val <- values(SVM_c_h)
SVM_rcp26_h_val <- values(SVM_rcp26_h)
SVM_rcp45_h_val <- values(SVM_rcp45_h)
SVM_rcp60_h_val <- values(SVM_rcp60_h)
SVM_rcp85_h_val <- values(SVM_rcp85_h)

SVM_h_val <- rbind(SVM_c_h_val, SVM_rcp26_h_val, SVM_rcp45_h_val, SVM_rcp60_h_val, SVM_rcp85_h_val)
SVM_h_stand <- decostand(SVM_h_val, "standardize", 2)


GLM_c_h_val <- values(GLM_c_h)
GLM_rcp26_h_val <- values(GLM_rcp26_h)
GLM_rcp45_h_val <- values(GLM_rcp45_h)
GLM_rcp60_h_val <- values(GLM_rcp60_h)
GLM_rcp85_h_val <- values(GLM_rcp85_h)

GLM_h_val <- rbind(GLM_c_h_val, GLM_rcp26_h_val, GLM_rcp45_h_val, GLM_rcp60_h_val, GLM_rcp85_h_val)
GLM_h_stand <- decostand(GLM_h_val, "standardize", 2)


## Host Plants


bioclim_c_p_val <- values(bioclim_c_p)
bioclim_rcp26_p_val <- values(bioclim_rcp26_p)
bioclim_rcp45_p_val <- values(bioclim_rcp45_p)
bioclim_rcp60_p_val <- values(bioclim_rcp60_p)
bioclim_rcp85_p_val <- values(bioclim_rcp85_p)

bioclim_p_val <- rbind(bioclim_c_p_val, bioclim_rcp26_p_val, bioclim_rcp45_p_val, bioclim_rcp60_p_val, bioclim_rcp85_p_val)
bioclim_p_stand <- decostand(bioclim_p_val, "standardize", 2)


gower_c_p_val <- values(gower_c_p)
gower_rcp26_p_val <- values(gower_rcp26_p)
gower_rcp45_p_val <- values(gower_rcp45_p)
gower_rcp60_p_val <- values(gower_rcp60_p)
gower_rcp85_p_val <- values(gower_rcp85_p)

gower_p_val <- rbind(gower_c_p_val, gower_rcp26_p_val, gower_rcp45_p_val, gower_rcp60_p_val, gower_rcp85_p_val)
gower_p_stand <- decostand(gower_p_val, "standardize", 2)


maha_c_p_val <- values(maha_c_p)
maha_rcp26_p_val <- values(maha_rcp26_p)
maha_rcp45_p_val <- values(maha_rcp45_p)
maha_rcp60_p_val <- values(maha_rcp60_p)
maha_rcp85_p_val <- values(maha_rcp85_p)

maha_p_val <- rbind(maha_c_p_val, maha_rcp26_p_val, maha_rcp45_p_val, maha_rcp60_p_val, maha_rcp85_p_val)
maha_p_stand <- decostand(maha_p_val, "standardize", 2)


maxent_c_p_val <- values(maxent_c_p)
maxent_rcp26_p_val <- values(maxent_rcp26_p)
maxent_rcp45_p_val <- values(maxent_rcp45_p)
maxent_rcp60_p_val <- values(maxent_rcp60_p)
maxent_rcp85_p_val <- values(maxent_rcp85_p)

maxent_p_val <- rbind(maxent_c_p_val, maxent_rcp26_p_val, maxent_rcp45_p_val, maxent_rcp60_p_val, maxent_rcp85_p_val)
maxent_p_stand <- decostand(maxent_p_val, "standardize", 2)


SVM_c_p_val <- values(SVM_c_p)
SVM_rcp26_p_val <- values(SVM_rcp26_p)
SVM_rcp45_p_val <- values(SVM_rcp45_p)
SVM_rcp60_p_val <- values(SVM_rcp60_p)
SVM_rcp85_p_val <- values(SVM_rcp85_p)

SVM_p_val <- rbind(SVM_c_p_val, SVM_rcp26_p_val, SVM_rcp45_p_val, SVM_rcp60_p_val, SVM_rcp85_p_val)
SVM_p_stand <- decostand(SVM_p_val, "standardize", 2)


GLM_c_p_val <- values(GLM_c_p)
GLM_rcp26_p_val <- values(GLM_rcp26_p)
GLM_rcp45_p_val <- values(GLM_rcp45_p)
GLM_rcp60_p_val <- values(GLM_rcp60_p)
GLM_rcp85_p_val <- values(GLM_rcp85_p)

GLM_p_val <- rbind(GLM_c_p_val, GLM_rcp26_p_val, GLM_rcp45_p_val, GLM_rcp60_p_val, GLM_rcp85_p_val)
GLM_p_stand <- decostand(GLM_p_val, "standardize", 2)

# ***************************************************************************************
## 11. Ensemble                                     ----

## huberi
suit <- data.frame(bioclim_h_stand, gower_h_stand, maha_h_stand, maxent_h_stand, SVM_h_stand, GLM_h_stand)
auc  <- c(bioclim_auc_h, gower_auc_h, maha_auc_h, maxent_auc_h, SVM_auc_h, GLM_auc_h)

# loop Ensemble for huberi

# Expected of ensembles for huberi:
# huberi current;
# huberi rcp26;
# huberi rcp45;
# huberi rcp60;
# huberi rcp85.

 
## host plants
suit <- data.frame(bioclim_p_stand, gower_p_stand, maha_p_stand, maxent_p_stand, SVM_p_stand, GLM_p_stand)
auc  <- c(bioclim_auc_p, gower_auc_p, maha_auc_p, maxent_auc_p, SVM_auc_p, GLM_auc_p)

# loop Ensemble for host plants

# Expected of ensembles for host plants:
# host plants current;
# host plants rcp26;
# host plants rcp45;
# host plants rcp60;
# host plants rcp85.


# ***************************************************************************************
## 12. Uncertainty Evaluation                       ----

## huberi
data   <- values(stack(bioclim_h, gower_h, maha_h, maxent_h, SVM_h, GLM_h))
method <- c(bioclim_method_h, gower_method_h, maha_method_h, maxent_method_h, SVM_method_h, GLM_method_h)
period <- c(bioclim_period_h, gower_period_h, maha_period_h, maxent_period_h, SVM_period_h, GLM_period_h)

# loop the ANOVA

## host plants
data   <- values(stack(bioclim_p, gower_p, maha_p, maxent_p, SVM_p, GLM_p))
method <- c(bioclim_method_p, gower_method_p, maha_method_p, maxent_method_p, SVM_method_p, GLM_method_p)
period <- c(bioclim_period_p, gower_period_p, maha_period_p, maxent_period_p, SVM_period_p, GLM_period_p)

# loop the ANOVA





