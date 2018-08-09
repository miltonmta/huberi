# Packages ######
# install.packages(c("tidyverse", "raster", "rgdal", "abind", "vegan", "maps", "mask","pcych", "kernlab", "dismo", "rJava")
require(tidyverse)
require(raster)
require(rgdal)
require(abind)
require(vegan)
require(maps)
require(mask)
require(kernlab)
require(dismo)
require(rJava)

# --- List of improvements to the scritp            ####

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
# browseURL("https://1drv.ms/f/s!ApJZaitgpPr7gZtfS9n9mU9DDzXQMg")

# The models projected for 2070 (average for 2061-2080) were obtain at "worldclim.com" by the spatial resolution of 2.5min (0.04º or ≈ 4.4km). We selected the variables that appear simultaneously in all the representative concentration pathways scnarios (RCP26, RCP45, RCP60, RCP80). The codes of the 11 GCMs utilized are: bc, cc, gs, hd, he, ip, mi, mr, mc, mg, no.
# browseURL("http://www.worldclim.org/cmip5_2.5m")


# 01. Read aogcms models ###########################################################################

### current ----
# Worldclim variables for current conditions (~1960-1990), by the spatial resolution of 2.5min. 
# browseURL("http://www.worldclim.org/current")

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

## Creating a matrix
current <- read_current(dir = "./data/climatic_vars/current")


## Creating a RasterBrick
current_spatial <- rasterFromXYZ(current) 


## rcp  ----
# Worldclim variables for the RCPs, projected to 2070(average for 2061-2080), by the spatial resolution of 2.5min. 
# browseURL("http://www.worldclim.org/cmip5_2.5m")

## wolrdclim GCM codes

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


# naming must be in the exact order of the origin directory. Best practice would be renaming the origin folders with the corresponding model.
model_names <- c("CCSM4", "IPSL-CM5A-LR", "MIROC-ESM")


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

# object type: array
rcp26 <- rcp26_list[["array"]]
rcp45 <- rcp45_list[["array"]]
rcp60 <- rcp60_list[["array"]]
rcp85 <- rcp85_list[["array"]]

# object type: RasterStack
rcp26_spatial <- rcp26_list[["rasters"]]
rcp26_spatial <- stack(rcp26_spatial[[1]], rcp26_spatial[[2]], rcp26_spatial[[3]])

rcp45_spatial <- rcp45_list[["rasters"]]
rcp45_spatial <- stack(rcp45_spatial[[1]], rcp45_spatial[[2]], rcp45_spatial[[3]])

rcp60_spatial <- rcp60_list[["rasters"]]
rcp60_spatial <- stack(rcp60_spatial[[1]], rcp60_spatial[[2]], rcp60_spatial[[3]])

rcp85_spatial <- rcp85_list[["rasters"]]
rcp85_spatial <- stack(rcp85_spatial[[1]], rcp85_spatial[[2]], rcp85_spatial[[3]])


rm(rcp26_list, rcp45_list, rcp60_list, rcp85_list)

# 02. Variable selection #########################################################################

### by varimax rotation ----

# # install.packages(c("psych", "GPArotation"), dependencies = TRUE)
# require(psych)
# # require(GPArotation)
# 
# fa.parallel(current[ , -c(1:2)], fa = 'fa') #scree plot
# current_fa <- fa(current[ , -c(1:2)], nfactors = 5, rotate = 'varimax')
# loadings <- loadings(current_fa)


### by Pearson's correlation ----

correlacao <- cor(current)
write.csv(correlacao, "./data/climatic_vars/selected/correlacao.csv")

# bio1, bio2, bio3, bio4, bio12, bio14, bio15, bio18, bio19.
# Selected variables ----

## bioclimatic variables descriptions

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


# 03. Saving selected variables ####################################################################

#### current

## 3.1a. current - saving as table ----

write.table(current[,c("x", "y", "bio01", "bio02", "bio03", "bio04", "bio12", "bio14", "bio15", "bio18", "bio19" )], "./data/climatic_vars/selected/current/current-select.txt", row.names = F, sep = " ")  

## 3.1b. current - saving as raster ----

variables <- as.factor(c("bio01", "bio02", "bio03", "bio04", "bio12", "bio14", "bio15", "bio18", "bio19"))
for (i in 1:length(variables))
{
  writeRaster (current_spatial[[i]], filename = paste0("./data/climatic_vars/selected/current/current-", variables[i], ".grd"), format = "raster")
}
rm(variables)

#>> current_select----
current_select <- stack(list.files("./data/climatic_vars/selected/current/",  pattern = ".grd$", full.names = TRUE))
# plot(current_select)

#### RCPs

### 3.2a. rcp - saving table ----

write.table(rcp26 [ ,c("x", "y", "bio01", "bio02", "bio03", "bio04", "bio12", "bio14", "bio15", "bio18", "bio19" ), ], "./data/climatic_vars/selected/rcp26/rcp26-select.txt", row.names = F, sep = "	")
write.table(rcp45 [ ,c("x", "y", "bio01", "bio02", "bio03", "bio04", "bio12", "bio14", "bio15", "bio18", "bio19" ), ], "./data/climatic_vars/selected/rcp45/rcp45-select.txt", row.names = F, sep = "	")
write.table(rcp60 [ ,c("x", "y", "bio01", "bio02", "bio03", "bio04", "bio12", "bio14", "bio15", "bio18", "bio19" ), ], "./data/climatic_vars/selected/rcp60/rcp60-select.txt", row.names = F, sep = "	")
write.table(rcp85 [ ,c("x", "y", "bio01", "bio02", "bio03", "bio04", "bio12", "bio14", "bio15", "bio18", "bio19" ), ], "./data/climatic_vars/selected/rcp85/rcp85-select.txt", row.names = F, sep = "	")


### 3.2b. rcp - saving as raster ----
# Creating rasters with the selected variables from each aogcm
variables <- as.factor(c("bio01.1", "bio02.1", "bio03.1", "bio04.1", "bio12.1", "bio14.1", "bio15.1", "bio18.1", "bio19.1", "bio01.2", "bio02.2", "bio03.2", "bio04.2", "bio12.2", "bio14.2", "bio15.2", "bio18.2", "bio19.2", "bio01.3", "bio02.3", "bio03.3", "bio04.3", "bio12.3", "bio14.3", "bio15.3", "bio18.3", "bio19.3"))

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
# Saving selected variables from CCSM4
for (i in 1:length(variables))
{
  writeRaster (rcp85_spatial[[i]], filename = paste0("./data/climatic_vars/selected/rcp85/rcp85-", variables[i], ".grd"), format = "raster")
}

rm(variables)


### Creating objects for each rcp

#>> rcpxx_select----
rcp26_select <- stack(list.files("./data/climatic_vars/selected/rcp26",  pattern = ".grd$", full.names = TRUE))

rcp45_select <- stack(list.files("./data/climatic_vars/selected/rcp45",  pattern = ".grd$", full.names = TRUE))

rcp60_select <- stack(list.files("./data/climatic_vars/selected/rcp60",  pattern = ".grd$", full.names = TRUE))

rcp85_select <- stack(list.files("./data/climatic_vars/selected/rcp85",  pattern = ".grd$", full.names = TRUE)) 

# par (mfrow = c(2,2))
# plot(rcp26_select)
# plot(rcp45_select)
# plot(rcp60_select)
# plot(rcp85_select)
# dev.off()

# 04. Occurrences data ###################################################################################

##  huberi
huberi <- read.table("./data/occurrences/huberi.txt", h = T)
huberi [1:5, ]
tail(huberi)

## host plants
# reading data
host_plants<- read.table("./data/occurrences/host-plants-species.txt", h = T)
host_plants [1:5, ]
str(host_plants) # object type matrix

# species names
sp <- gsub("C[1-9]","", host_plants$species)
sp <- unique(sp)
# [1] "Ipomoea_asarifolia" "Ipomoea_bahiensis" 
# [3] "Ipomoea_cairica"    "Ipomoea_indica"    
# [5] "Ipomoea_nil"        "Ipomoea_purpurea"  
# [7] "Merremia_aegyptia" 

## Plotting 
plot(current_select$bio01)
points(huberi[,-1], pch = "*", col = "blue")
points(host_plants[, -1], pch = "*", col = "red")
# points(host_plants[host_plants[, 1] == "Ipomoea_asarifolia", ],  pch = "*", col = "red")

## extracting variables based on the occurrences data cells
huberi_cell <- cellFromXY(current_select, huberi[, -1])
duplicated(huberi_cell)
huberi_cell <- unique(huberi_cell)
huberi_var <- raster::extract(current_select, huberi_cell)
huberi_var <- na.omit(huberi_var)
#>> huberi_var----


host_plants_cell <- cellFromXY(current_select, host_plants [, -1])
duplicated(host_plants_cell)
host_plants_cell <- unique(host_plants_cell)
host_plants_var <- raster::extract(current_select, host_plants_cell)
host_plants_var <- na.omit(host_plants_var)
#>> host_plants_var----

## saving data
write.table(huberi_var,      "./data/occurrences/huberi-var.txt", row.names = F, sep = " ") 
write.table(host_plants_var, "./data/occurrences/host-plants-var.txt", row.names = F, sep = " ") 

huberi_var <- read.table("./data/occurrences/huberi-var.txt", sep = " ")
host_plants_var <- read.table("./data/occurrences/host-plants-var.txt", sep = " ")

## 05. Background Sampling ##############################################################################

## CURRENT
# # creating the background with the Neotropic study area
# neot <- raster::extract(current_select, 1:ncell(current_select))
# neot.coords <- xyFromCell(current_select, 1:ncell(current_select)) 
# neot <- cbind(neot.coords, cells = 1:ncell(current_select), neot)
# neot <- na.omit(neot)
# 
# # sampling with huberi
# back_id_huberi <- sample(1:nrow(neot), nrow(huberi_var))
# back_huberi <- neot[back_id_huberi, ]
# points(back_huberi[, "x"], back_huberi[, "y"], pch = "*", col = 'black')

#>> back_huberi----
back_id_huberi <- sample(1:nrow(current_select), nrow(huberi_var))
coord_huberi   <- xyFromCell(current_select, back_id_huberi)
back_huberi    <- extract(current_select, back_id_huberi)
back_huberi    <- cbind(coord_huberi,back_huberi)
points(back_huberi[, "x"], back_huberi[, "y"], pch = "*", col = 'black')


# sampling with host plants
back_id_plants <- sample(1:nrow(neot), nrow(host_plants_var))
back_plants <- neot[back_id_plants, ]
points(back_plants[, "x"], back_plants[, "y"], pch = "*", col = 'magenta') # número incorreto de dimensões
#>> back_plants----

# saving background files
write.table(back_huberi, "./data/occurrences/Background-random-huberi.txt", row.names = F, sep = " ") 
write.table(back_plants, "./data/occurrences/Background-random-plants.txt", row.names = F, sep = " ")

rm(list = ls())

# 06. Modelling Predictions ##############################################################

modelling <- function(occurrence_huberi, 
                      occurrence_plants, 
                      background_huberi, 
                      background_plants, 
                      biovar_current,
                      biovar_rcp26,
                      biovar_rcp45,
                      biovar_rcp60,
                      biovar_rcp85 ,
                      cross_validation = ...)
{
  
  output_current <- output_rcp26 <- output_rcp45 <- output_rcp60 <- output_rcp85 <- NULL
  
  # > model_names
  # [1] "CCSM4"        "IPSL-CM5A-LR" "MIROC-ESM"   
  # > 
  AOGCMs <- c("CCSM4", "IPSL-CM5A-LR", "MIROC-ESM")
  
  for (j in AOGCMs)
  {
    ### OPEN "j" ----
    
    ### Reading the selected bioclimatic variables
    
    # creating objects for processing all aogcm models
    AOGCM_CURRENT <- paste("be_biovar_", j, biovar_current, sep="")
    AOGCM_RCP26   <- paste("be_biovar_", j, biovar_rcp26,   sep="")
    AOGCM_RCP45   <- paste("be_biovar_", j, biovar_rcp45,   sep="")
    AOGCM_RCP60   <- paste("be_biovar_", j, biovar_rcp60,   sep="")
    AOGCM_RCP85   <- paste("be_biovar_", j, biovar_rcp85,   sep="")
    
    # vars_current  <- read.table(AOGCM_CURRENT, h = T)
    # vars_rcp26    <- read.table(AOGCM_RCP26,   h = T)
    # vars_rcp45    <- read.table(AOGCM_RCP45,   h = T)
    # vars_rcp60    <- read.table(AOGCM_RCP60,   h = T)
    # vars_rcp85    <- read.table(AOGCM_RCP85,   h = T)
    # 
    # # rasterizing the variables
    # gridded(vars_current) <- ~ x + y
    # gridded(vars_rcp26)   <- ~ x + y
    # gridded(vars_rcp45)   <- ~ x + y
    # gridded(vars_rcp60)   <- ~ x + y
    # gridded(vars_rcp85)   <- ~ x + y
    # 
    # # creating the bioclimatic variables objects
    # current_select <- stack(vars_current)
    # rcp26_select   <- stack(vars_rcp26) 
    # rcp45_select   <- stack(vars_rcp45)
    # rcp60_select   <- stack(vars_rcp60)
    # rcp85_select   <- stack(vars_rcp85)
    
    # creating the bioclimatic variables objects
    current_select <- stack(list.files(AOGCM_CURRENT,  pattern = ".grd$", full.names = TRUE))
    rcp26_select   <- stack(list.files(AOGCM_RCP26,    pattern = ".grd$", full.names = TRUE))
    rcp45_select   <- stack(list.files(AOGCM_RCP45,    pattern = ".grd$", full.names = TRUE))
    rcp60_select   <- stack(list.files(AOGCM_RCP60,    pattern = ".grd$", full.names = TRUE))
    rcp85_select   <- stack(list.files(AOGCM_RCP85,    pattern = ".grd$", full.names = TRUE))
    
    
    ### Reading occurrences data
    
    occur_h <- read.table(occurrence_huberi, h = T)
    occur_p <- read.table(occurrence_plants, h = T)
    back_h  <- read.table(background_huberi, h = T)
    back_p  <- read.table(background_plants, h = T)
    
    occur_h <- read.table("./data/occurrences/huberi-var.txt", h = T)
    occur_p <- read.table("./data/occurrences/host-plants-var.txt", h = T)
    back_h  <- read.table("./data/occurrences/Background-random-huberi.txt", h = T)
    back_p  <- read.table("./data/occurrences/Background-random-plants.txt", h = T)
      
    ### Creating objects for saving partial results for each cross validation loop
    ## huberi
    bioclim_c_h <- gower_c_h <- maha_c_h <- maxent_c_h <- SVM_c_h <- GLM_c_h <- stack()
    bioclim_rcp26_h <- gower_rcp26_h <- maha_rcp26_h <- maxent_rcp26_h <- SVM_rcp26_h <- GLM_rcp26_h <- stack() 
    bioclim_rcp45_h <- gower_rcp45_h <- maha_rcp45_h <- maxent_rcp45_h <- SVM_rcp45_h <- GLM_rcp45_h <- stack()
    bioclim_rcp60_h <- gower_rcp60_h <- maha_rcp60_h <- maxent_rcp60_h <- SVM_rcp60_h <- GLM_rcp60_h <- stack()
    bioclim_rcp85_h <- gower_rcp85_h <- maha_rcp85_h <- maxent_rcp85_h <- SVM_rcp85_h <- GLM_rcp85_h <- stack()
    
    bioclim_e_h <- gower_e_h <- maha_e_h <- maxent_e_h <- SVM_e_h <- GLM_e_h <- NULL
    bioclim_t_h <- gower_t_h <- maha_t_h <- maxent_t_h <- SVM_t_h <- GLM_t_h <- NULL
    
    ## host plants
    bioclim_c_p <- gower_c_p <- maha_c_p <- maxent_c_p <- SVM_c_p <- GLM_c_p <- stack()
    bioclim_rcp26_p <- gower_rcp26_p <- maha_rcp26_p <- maxent_rcp26_p <- SVM_rcp26_p <- GLM_rcp26_p <- stack() 
    bioclim_rcp45_p <- gower_rcp45_p <- maha_rcp45_p <- maxent_rcp45_p <- SVM_rcp45_p <- GLM_rcp45_p <- stack()
    bioclim_rcp60_p <- gower_rcp60_p <- maha_rcp60_p <- maxent_rcp60_p <- SVM_rcp60_p <- GLM_rcp60_p <- stack()
    bioclim_rcp85_p <- gower_rcp85_p <- maha_rcp85_p <- maxent_rcp85_p <- SVM_rcp85_p <- GLM_rcp85_p <- stack()
    
    bioclim_e_p <- gower_e_p <- maha_e_p <- maxent_e_p <- SVM_e_p <- GLM_e_p <- NULL
    bioclim_t_p <- gower_t_p <- maha_t_p <- maxent_t_p <- SVM_t_p <- GLM_t_p <- NULL
    
    
    for (i in 1:cross_validation)
    {
      ### OPEN "i" ----
      # loop - cross-validation
      
      ### creating trainning-testing subsets
    
      ## huberi
      sample_occur_h <- sample(1:nrow(occur_h), round(0.75 * nrow(occur_h), 0))
      # sample_back_h  <- sample(1:nrow(back_h),  round(0.75 * nrow(back_h),  0))
      
      training_h <- prepareData(x = current_select, p = occur_h[sample_occur_h,  1:2], b = back_h[sample_occur_h,  1:2])
      testing_h  <- prepareData(x = current_select, p = occur_h[-sample_occur_h, 1:2], b = back_h[-sample_occur_h, 1:2])
      
      ## host plants
      sample_occur_p <- sample(1:nrow(occur_p), round(0.75 * nrow(occur_p)))
      sample_back_p  <- sample(1:nrow(back_p), round(0.75 * nrow(back_p)))
      
      training_p <- prepareData(x = current_select, p = occur_p[sample_occur_p,  1:2], b = back_p[sample_back_p,  1:2], xy = T)
      testing_p  <- prepareData(x = current_select, p = occur_p[-sample_occur_p, 1:2], b = back_p[-sample_back_p, 1:2], xy = T)
      
      ### Bioclim -------------------
      
      ## huberi
      # ajusting models
      bioclim_model_h <- bioclim(training_h[training_h[, "pb"] == 1, -1])
      
      # making predictions
      bioclim_c_h <- stack(bioclim_c_h, predict(object = bioclim_model_h, x = current_select))
      bioclim_rcp26_h <- stack(bioclim_rcp26_h, predict(object = bioclim_model_h, x = rcp26_select))
      bioclim_rcp45_h <- stack(bioclim_rcp45_h, predict(object = bioclim_model_h, x = rcp45_select))
      bioclim_rcp60_h <- stack(bioclim_rcp60_h, predict(object = bioclim_model_h, x = rcp60_select))
      bioclim_rcp85_h <- stack(bioclim_rcp85_h, predict(object = bioclim_model_h, x = rcp85_select))
      
      # Evaluating models
      bioclim_eval_h <- evaluate(p=testing_h[testing_h[, "pb"] == 1, -1], a = testing_h[testing_h[, "pb"] == 0, -1], model = bioclim_model_h)
      
      bioclim_e_h <- c(bioclim_e_h, bioclim_eval@auc)
      bioclim_t_h <- c(bioclim_t_h, threshold(bioclim_eval_h, "spec_sens"))
      
      ## host plants
      # ajusting models
      bioclim_model_p <- bioclim(training_p[training_p[,"pb"] == 1, -c(1:3)])
      
      # making predictions
      bioclim_c_p <- stack(bioclim_c_p, predict(object = bioclim_model_p, x = current_select))
      bioclim_rcp26_p <- stack(bioclim_rcp26_p, predict(object = bioclim_model_p, x = rcp26_select))
      bioclim_rcp45_p <- stack(bioclim_rcp45_p, predict(object = bioclim_model_p, x = rcp45_select))
      bioclim_rcp60_p <- stack(bioclim_rcp60_p, predict(object = bioclim_model_p, x = rcp60_select))
      bioclim_rcp85_p <- stack(bioclim_rcp85_p, predict(object = bioclim_model_p, x = rcp85_select))
      
      # Evaluating models
      bioclim_eval_p <- evaluate(p=testing_p[testing_p[, "pb"] == 1, -1], a = testing_p[testing_p[, "pb"] == 0, -1], model = bioclim_model_p)
      
      bioclim_e_p <- c(bioclim_e_p, bioclim_eval@auc)
      bioclim_t_p <- c(bioclim_t_p, threshold(bioclim_eval_p, "spec_sens"))
      
      ### Gower ---------------------
      ## huberi
      # ajusting models
      gower_model_h <- domain(training_h[training_h[,"pb"] == 1, -c(1:3)])
      
      # making predictions
      gower_c_h <- stack(gower_c_h, predict(object = gower_model_h, x = current_select))
      gower_rcp26_h <- stack(gower_rcp26_h, predict(object = gower_model_h, x = rcp26_select))
      gower_rcp45_h <- stack(gower_rcp45_h, predict(object = gower_model_h, x = rcp45_select))
      gower_rcp60_h <- stack(gower_rcp60_h, predict(object = gower_model_h, x = rcp60_select))
      gower_rcp85_h <- stack(gower_rcp85_h, predict(object = gower_model_h, x = rcp85_select))
      
      # Evaluating models
      gower_eval_h <- evaluate(p=testing_h[testing_h[, "pb"] == 1, -1], a = testing_h[testing_h[, "pb"] == 0, -1], model = gower_model_h)
      
      gower_e_h <- c(gower_e_h, gower_eval@auc)
      gower_t_h <- c(gower_t_h, threshold(gower_eval_h, "spec_sens"))
      
      ## host plants
      # ajusting models
      gower_model_p <- domain(training_p[training_p[,"pb"] == 1, -c(1:3)])
      
      # making predictions
      gower_c_p <- stack(gower_c_p, predict(object = gower_model_p, x = current_select))
      gower_rcp26_p <- stack(gower_rcp26_p, predict(object = gower_model_p, x = rcp26_select))
      gower_rcp45_p <- stack(gower_rcp45_p, predict(object = gower_model_p, x = rcp45_select))
      gower_rcp60_p <- stack(gower_rcp60_p, predict(object = gower_model_p, x = rcp60_select))
      gower_rcp85_p <- stack(gower_rcp85_p, predict(object = gower_model_p, x = rcp85_select))
      
      # Evaluating models
      gower_eval_p <- evaluate(p=testing_p[testing_p[, "pb"] == 1, -1], a = testing_p[testing_p[, "pb"] == 0, -1], model = gower_model_p)
      
      gower_e_p <- c(gower_e_p, gower_eval@auc)
      gower_t_p <- c(gower_t_p, threshold(gower_eval_p, "spec_sens"))
      
      
      ### Maha ----
      ## huberi
      # ajusting models
      maha_model_h <- mahal(training_h[training_h[,"pb"] == 1, -c(1:3)])
      
      # making predictions
      maha_c_h <- stack(maha_c_h, predict(object = maha_model_h, x = current_select))
      maha_rcp26_h <- stack(maha_rcp26_h, predict(object = maha_model_h, x = rcp26_select))
      maha_rcp45_h <- stack(maha_rcp45_h, predict(object = maha_model_h, x = rcp45_select))
      maha_rcp60_h <- stack(maha_rcp60_h, predict(object = maha_model_h, x = rcp60_select))
      maha_rcp85_h <- stack(maha_rcp85_h, predict(object = maha_model_h, x = rcp85_select))
      
      # Evaluating models
      maha_eval_h <- evaluate(p=testing_h[testing_h[, "pb"] == 1, -1], a = testing_h[testing_h[, "pb"] == 0, -1], model = maha_model_h)
      
      maha_e_h <- c(maha_e_h, maha_eval@auc)
      maha_t_h <- c(maha_t_h, threshold(maha_eval_h, "spec_sens"))
      
      ## host plants
      # ajusting models
      maha_model_p <- mahal(training_p[training_p[,"pb"] == 1, -c(1:3)])
      
      # making predictions
      maha_c_p <- stack(maha_c_p, predict(object = maha_model_p, x = current_select))
      maha_rcp26_p <- stack(maha_rcp26_p, predict(object = maha_model_p, x = rcp26_select))
      maha_rcp45_p <- stack(maha_rcp45_p, predict(object = maha_model_p, x = rcp45_select))
      maha_rcp60_p <- stack(maha_rcp60_p, predict(object = maha_model_p, x = rcp60_select))
      maha_rcp85_p <- stack(maha_rcp85_p, predict(object = maha_model_p, x = rcp85_select))
      
      # Evaluating models
      maha_eval_p <- evaluate(p=testing_p[testing_p[, "pb"] == 1, -1], a = testing_p[testing_p[, "pb"] == 0, -1], model = maha_model_p)
      
      maha_e_p <- c(maha_e_p, maha_eval@auc)
      maha_t_p <- c(maha_t_p, threshold(maha_eval_p, "spec_sens"))
      
      
      
      ### Maxent ----
    
      
      ## huberi
      # ajusting models
      Sys.setenv(NOAWT = TRUE)
      maxent_model_h <- maxent(x = training_h[, -1], p = training_h[, 1])
      
      # making predictions
      maxent_c_h <- stack(maxent_c_h, predict(object = maxent_model_h, x = current_select))
      maxent_rcp26_h <- stack(maxent_rcp26_h, predict(object = maxent_model_h, x = rcp26_select))
      maxent_rcp45_h <- stack(maxent_rcp45_h, predict(object = maxent_model_h, x = rcp45_select))
      maxent_rcp60_h <- stack(maxent_rcp60_h, predict(object = maxent_model_h, x = rcp60_select))
      maxent_rcp85_h <- stack(maxent_rcp85_h, predict(object = maxent_model_h, x = rcp85_select))
      
      # Evaluating models
      maxent_eval_h <- evaluate(p=testing_h[testing_h[, "pb"] == 1, -1], a = testing_h[testing_h[, "pb"] == 0, -1], model = maxent_model_h)
      
      maxent_e_h <- c(maxent_e_h, maxent_eval@auc)
      maxent_t_h <- c(maxent_t_h, threshold(maxent_eval_h, "spec_sens"))
      
      ## host plants
      # ajusting models
      Sys.setenv(NOAWT = TRUE)
      maxent_model_p <- maxent(x = training_p[, -1], p = training_p[, 1])
      
      # making predictions
      maxent_c_p <- stack(maxent_c_p, predict(object = maxent_model_p, x = current_select))
      maxent_rcp26_p <- stack(maxent_rcp26_p, predict(object = maxent_model_p, x = rcp26_select))
      maxent_rcp45_p <- stack(maxent_rcp45_p, predict(object = maxent_model_p, x = rcp45_select))
      maxent_rcp60_p <- stack(maxent_rcp60_p, predict(object = maxent_model_p, x = rcp60_select))
      maxent_rcp85_p <- stack(maxent_rcp85_p, predict(object = maxent_model_p, x = rcp85_select))
      
      # Evaluating models
      maxent_eval_p <- evaluate(p=testing_p[testing_p[, "pb"] == 1, -1], a = testing_p[testing_p[, "pb"] == 0, -1], model = maxent_model_p)
      
      maxent_e_p <- c(maxent_e_p, maxent_eval@auc)
      maxent_t_p <- c(maxent_t_p, threshold(maxent_eval_p, "spec_sens"))
      
      
      ### SVM ----
      # ajusting models
      SVM_model_h <- ksvm(pb ~ bio+bio..., data = training_h)
      
      # making predictions
      SVM_c_h <- stack(SVM_c_h, predict(object = SVM_model_h, x = current_select))
      SVM_rcp26_h <- stack(SVM_rcp26_h, predict(object = SVM_model_h, x = rcp26_select))
      SVM_rcp45_h <- stack(SVM_rcp45_h, predict(object = SVM_model_h, x = rcp45_select))
      SVM_rcp60_h <- stack(SVM_rcp60_h, predict(object = SVM_model_h, x = rcp60_select))
      SVM_rcp85_h <- stack(SVM_rcp85_h, predict(object = SVM_model_h, x = rcp85_select))
      
      # Evaluating models
      SVM_eval_h <- evaluate(p=testing_h[testing_h[, "pb"] == 1, -1], a = testing_h[testing_h[, "pb"] == 0, -1], model = SVM_model_h)
      
      SVM_e_h <- c(SVM_e_h, SVM_eval@auc)
      SVM_t_h <- c(SVM_t_h, threshold(SVM_eval_h, "spec_sens"))
      
      ## host plants
      # ajusting models
      SVM_model_h <- ksvm(pb ~ bio+bio..., data = training_h)
      
      # making predictions
      SVM_c_p <- stack(SVM_c_p, predict(object = SVM_model_p, x = current_select))
      SVM_rcp26_p <- stack(SVM_rcp26_p, predict(object = SVM_model_p, x = rcp26_select))
      SVM_rcp45_p <- stack(SVM_rcp45_p, predict(object = SVM_model_p, x = rcp45_select))
      SVM_rcp60_p <- stack(SVM_rcp60_p, predict(object = SVM_model_p, x = rcp60_select))
      SVM_rcp85_p <- stack(SVM_rcp85_p, predict(object = SVM_model_p, x = rcp85_select))
      
      # Evaluating models
      SVM_eval_p <- evaluate(p=testing_p[testing_p[, "pb"] == 1, -1], a = testing_p[testing_p[, "pb"] == 0, -1], model = SVM_model_p)
      
      SVM_e_p <- c(SVM_e_p, SVM_eval@auc)
      SVM_t_p <- c(SVM_t_p, threshold(SVM_eval_p, "spec_sens"))
      
      
      ### GLM ----
      GLM_model_h <- glm(pb ~ bio+bio..., data = training_h, family = binomial (link = "logit"))
      
      # making predictions
      GLM_c_h <- stack(GLM_c_h, predict(object = GLM_model_h, x = current_select))
      GLM_rcp26_h <- stack(GLM_rcp26_h, predict(object = GLM_model_h, x = rcp26_select))
      GLM_rcp45_h <- stack(GLM_rcp45_h, predict(object = GLM_model_h, x = rcp45_select))
      GLM_rcp60_h <- stack(GLM_rcp60_h, predict(object = GLM_model_h, x = rcp60_select))
      GLM_rcp85_h <- stack(GLM_rcp85_h, predict(object = GLM_model_h, x = rcp85_select))
      
      # Evaluating models
      GLM_eval_h <- evaluate(p=testing_h[testing_h[, "pb"] == 1, -1], a = testing_h[testing_h[, "pb"] == 0, -1], model = GLM_model_h)
      
      GLM_e_h <- c(GLM_e_h, GLM_eval@auc)
      GLM_t_h <- c(GLM_t_h, threshold(GLM_eval_h, "spec_sens"))
      
      ## host plants
      # ajusting models
      GLM_model_h <-glm(pb ~ bio+bio..., data = training_h, family = binomial (link = "logit"))
      
      # making predictions
      GLM_c_p <- stack(GLM_c_p, predict(object = GLM_model_p, x = current_select))
      GLM_rcp26_p <- stack(GLM_rcp26_p, predict(object = GLM_model_p, x = rcp26_select))
      GLM_rcp45_p <- stack(GLM_rcp45_p, predict(object = GLM_model_p, x = rcp45_select))
      GLM_rcp60_p <- stack(GLM_rcp60_p, predict(object = GLM_model_p, x = rcp60_select))
      GLM_rcp85_p <- stack(GLM_rcp85_p, predict(object = GLM_model_p, x = rcp85_select))
      
      # Evaluating models
      GLM_eval_p <- evaluate(p=testing_p[testing_p[, "pb"] == 1, -1], a = testing_p[testing_p[, "pb"] == 0, -1], model = GLM_model_p)
      
      GLM_e_p <- c(GLM_e_p, GLM_eval@auc)
      GLM_t_p <- c(GLM_t_p, threshold(GLM_eval_p, "spec_sens"))
    
    }# CLOSE "i" ----
    
    
    ## writing predictions as raster
    # huberi
    writeRaster(bioclim_c_h, "./data/outputs/bioclim_c_h.bil", format = "EHdr")
    writeRaster(bioclim_rcp26_h, "./data/outputs/bioclim_rcp26_h.bil", format = "EHdr")
    writeRaster(bioclim_rcp45_h, "./data/outputs/bioclim_rcp45_h.bil", format = "EHdr")
    writeRaster(bioclim_rcp60_h, "./data/outputs/bioclim_rcp60_h.bil", format = "EHdr")
    writeRaster(bioclim_rcp85_h, "./data/outputs/bioclim_rcp85_h.bil", format = "EHdr")
    
    writeRaster(gower_c_h, "./data/outputs/gower_c_h.bil", format = "EHdr")
    writeRaster(gower_rcp26_h, "./data/outputs/gower_rcp26_h.bil", format = "EHdr")
    writeRaster(gower_rcp45_h, "./data/outputs/gower_rcp45_h.bil", format = "EHdr")
    writeRaster(gower_rcp60_h, "./data/outputs/gower_rcp60_h.bil", format = "EHdr")
    writeRaster(gower_rcp85_h, "./data/outputs/gower_rcp85_h.bil", format = "EHdr")
    
    writeRaster(maha_c_h, "./data/outputs/maha_c_h.bil", format = "EHdr")
    writeRaster(maha_rcp26_h, "./data/outputs/maha_rcp26_h.bil", format = "EHdr")
    writeRaster(maha_rcp45_h, "./data/outputs/maha_rcp45_h.bil", format = "EHdr")
    writeRaster(maha_rcp60_h, "./data/outputs/maha_rcp60_h.bil", format = "EHdr")
    writeRaster(maha_rcp85_h, "./data/outputs/maha_rcp85_h.bil", format = "EHdr")
    
    writeRaster(maxent_c_h, "./data/outputs/maxent_c_h.bil", format = "EHdr")
    writeRaster(maxent_rcp26_h, "./data/outputs/maxent_rcp26_h.bil", format = "EHdr")
    writeRaster(maxent_rcp45_h, "./data/outputs/maxent_rcp45_h.bil", format = "EHdr")
    writeRaster(maxent_rcp60_h, "./data/outputs/maxent_rcp60_h.bil", format = "EHdr")
    writeRaster(maxent_rcp85_h, "./data/outputs/maxent_rcp85_h.bil", format = "EHdr")
    
    writeRaster(SVM_c_h, "./data/outputs/SVM_c_h.bil", format = "EHdr")
    writeRaster(SVM_rcp26_h, "./data/outputs/SVM_rcp26_h.bil", format = "EHdr")
    writeRaster(SVM_rcp45_h, "./data/outputs/SVM_rcp45_h.bil", format = "EHdr")
    writeRaster(SVM_rcp60_h, "./data/outputs/SVM_rcp60_h.bil", format = "EHdr")
    writeRaster(SVM_rcp85_h, "./data/outputs/SVM_rcp85_h.bil", format = "EHdr")
    
    writeRaster(GLM_c_h, "./data/outputs/GLM_c_h.bil", format = "EHdr")
    writeRaster(GLM_rcp26_h, "./data/outputs/GLM_rcp26_h.bil", format = "EHdr")
    writeRaster(GLM_rcp45_h, "./data/outputs/GLM_rcp45_h.bil", format = "EHdr")
    writeRaster(GLM_rcp60_h, "./data/outputs/GLM_rcp60_h.bil", format = "EHdr")
    writeRaster(GLM_rcp85_h, "./data/outputs/GLM_rcp85_h.bil", format = "EHdr")
    
    write.table(data.frame(bioclim = bioclim_e_h ,gower = gower_e_h, maha = maha_e_h, maxent = maxent_e_h, SVM = SVM_e_h, GLM = GLM_e_h), "./data/outputs/AUC_huberi.txt", sep = "\t", row.names = F)
    
    write.table(data.frame(bioclim = bioclim_t_h ,gower = gower_t_h, maha = maha_t_h, maxent = maxent_t_h, SVM = SVM_t_h, GLM = GLM_t_h), "./data/outputs/Threshold_huberi.txt", sep = "\t", row.names = F)
    
    # Host plants
    writeRaster(bioclim_c_p, "./data/outputs/bioclim_c_p.bil", format = "EHdr")
    writeRaster(bioclim_rcp26_p, "./data/outputs/bioclim_rcp26_p.bil", format = "EHdr")
    writeRaster(bioclim_rcp45_p, "./data/outputs/bioclim_rcp45_p.bil", format = "EHdr")
    writeRaster(bioclim_rcp60_p, "./data/outputs/bioclim_rcp60_p.bil", format = "EHdr")
    writeRaster(bioclim_rcp85_p, "./data/outputs/bioclim_rcp85_p.bil", format = "EHdr")
    
    writeRaster(gower_c_p, "./data/outputs/gower_c_p.bil", format = "EHdr")
    writeRaster(gower_rcp26_p, "./data/outputs/gower_rcp26_p.bil", format = "EHdr")
    writeRaster(gower_rcp45_p, "./data/outputs/gower_rcp45_p.bil", format = "EHdr")
    writeRaster(gower_rcp60_p, "./data/outputs/gower_rcp60_p.bil", format = "EHdr")
    writeRaster(gower_rcp85_p, "./data/outputs/gower_rcp85_p.bil", format = "EHdr")
    
    writeRaster(maha_c_p, "./data/outputs/maha_c_p.bil", format = "EHdr")
    writeRaster(maha_rcp26_p, "./data/outputs/maha_rcp26_p.bil", format = "EHdr")
    writeRaster(maha_rcp45_p, "./data/outputs/maha_rcp45_p.bil", format = "EHdr")
    writeRaster(maha_rcp60_p, "./data/outputs/maha_rcp60_p.bil", format = "EHdr")
    writeRaster(maha_rcp85_p, "./data/outputs/maha_rcp85_p.bil", format = "EHdr")
    
    writeRaster(maxent_c_p, "./data/outputs/maxent_c_p.bil", format = "EHdr")
    writeRaster(maxent_rcp26_p, "./data/outputs/maxent_rcp26_p.bil", format = "EHdr")
    writeRaster(maxent_rcp45_p, "./data/outputs/maxent_rcp45_p.bil", format = "EHdr")
    writeRaster(maxent_rcp60_p, "./data/outputs/maxent_rcp60_p.bil", format = "EHdr")
    writeRaster(maxent_rcp85_p, "./data/outputs/maxent_rcp85_p.bil", format = "EHdr")
    
    writeRaster(SVM_c_p, "./data/outputs/SVM_c_p.bil", format = "EHdr")
    writeRaster(SVM_rcp26_p, "./data/outputs/SVM_rcp26_p.bil", format = "EHdr")
    writeRaster(SVM_rcp45_p, "./data/outputs/SVM_rcp45_p.bil", format = "EHdr")
    writeRaster(SVM_rcp60_p, "./data/outputs/SVM_rcp60_p.bil", format = "EHdr")
    writeRaster(SVM_rcp85_p, "./data/outputs/SVM_rcp85_p.bil", format = "EHdr")
    
    writeRaster(GLM_c_p, "./data/outputs/GLM_c_p.bil", format = "EHdr")
    writeRaster(GLM_rcp26_p, "./data/outputs/GLM_rcp26_p.bil", format = "EHdr")
    writeRaster(GLM_rcp45_p, "./data/outputs/GLM_rcp45_p.bil", format = "EHdr")
    writeRaster(GLM_rcp60_p, "./data/outputs/GLM_rcp60_p.bil", format = "EHdr")
    writeRaster(GLM_rcp85_p, "./data/outputs/GLM_rcp85_p.bil", format = "EHdr")
    
    write.table(data.frame(bioclim = bioclim_e_p ,gower = gower_e_p, maha = maha_e_p, maxent = maxent_e_p, SVM = SVM_e_p, GLM = GLM_e_p), "./data/outputs/AUC_plants.txt", sep = "\t", row.names = F)
    
    write.table(data.frame(bioclim = bioclim_t_p ,gower = gower_t_p, maha = maha_t_p, maxent = maxent_t_p, SVM = SVM_t_p, GLM = GLM_t_p), "./data/outputs/Threshold_plants.txt", sep = "\t", row.names = F)
    
  }# CLOSE "j"  ----
  
  
}# closes the function modelling

# ! running our model ####

modelling (occurrence_huberi = "./data/occurrences/huberi-var.txt", 
           occurrence_plants = "./data/occurrences/host-plants-var.txt", 
           background_huberi = "./data/occurrences/Background-random-huberi.txt", 
           background_plants = "./data/occurrences/Background-random-plants.txt",
           biovar_current    = "./data/climatic_vars/selected/current/current-select.txt",
           biovar_rcp26      = "./data/climatic_vars/selected/rcp26/rcp26-select.txt",
           biovar_rcp45      = "./data/climatic_vars/selected/rcp45/rcp45-select.txt",
           biovar_rcp60      = "./data/climatic_vars/selected/rcp60/rcp60-select.txt",
           biovar_rcp85      = "./data/climatic_vars/selected/rcp85/rcp85-select.txt",
           cross_validation  = 100)

# 07. Selecting models (auc) ############################################################################

### Huberi
auc_h <- read.table("./data/outputs/AUC_huberi.txt", h = T)

## bioclim

#Selecting models with auc values ≥ 0.7
bioclim_c_h <- stack("./data/outputs/bioclim_c_h.bil")[[which(auc_h[,"bioclim"] >= 0.7)]]
bioclim_rcp26_h <- stack("./data/outputs/bioclim_rcp26_h.bil")[[which(auc_h[,"bioclim"] >= 0.7)]]
bioclim_rcp45_h <- stack("./data/outputs/bioclim_rcp45_h.bil")[[which(auc_h[,"bioclim"] >= 0.7)]]
bioclim_rcp60_h <- stack("./data/outputs/bioclim_rcp60_h.bil")[[which(auc_h[,"bioclim"] >= 0.7)]]
bioclim_rcp85_h <- stack("./data/outputs/bioclim_rcp85_h.bil")[[which(auc_h[,"bioclim"] >= 0.7)]]

bioclim_h <- stack(bioclim_c_h, bioclim_rcp26_h, bioclim_rcp45_h, bioclim_rcp60_h, bioclim_rcp85_h)
bioclim_auc_h <- auc_h[which(auc_h[,"bioclim"] >= 0.7),"bioclim"]

#  Creating ANOVA Factors to be used at Uncertainty Evaluation
bioclim_method_h <- rep("bioclim", nlayers(bioclim_h))
bioclim_period_h <- rep(c("pres","fut"), each = nlayers(bioclim_h))

## gower

#Selecting models with auc values ≥ 0.7
gower_c_h <- stack("./data/outputs/gower_c_h.bil")[[which(auc_h[,"gower"] >= 0.7)]]
gower_rcp26_h <- stack("./data/outputs/gower_rcp26_h.bil")[[which(auc_h[,"gower"] >= 0.7)]]
gower_rcp45_h <- stack("./data/outputs/gower_rcp45_h.bil")[[which(auc_h[,"gower"] >= 0.7)]]
gower_rcp60_h <- stack("./data/outputs/gower_rcp60_h.bil")[[which(auc_h[,"gower"] >= 0.7)]]
gower_rcp85_h <- stack("./data/outputs/gower_rcp85_h.bil")[[which(auc_h[,"gower"] >= 0.7)]]

gower_h <- stack(gower_c_h, gower_rcp26_h, gower_rcp45_h, gower_rcp60_h, gower_rcp85_h)
gower_auc_h <- auc_h[which(auc_h[,"gower"] >= 0.7),"gower"]

#  Creating ANOVA Factors to be used at Uncertainty Evaluation
gower_method_h <- rep("gower", nlayers(gower_h))
gower_period_h <- rep(c("pres","fut"), each = nlayers(gower_h))

## maha

#Selecting models with auc values ≥ 0.7
maha_c_h <- stack("./data/outputs/maha_c_h.bil")[[which(auc_h[,"maha"] >= 0.7)]]
maha_rcp26_h <- stack("./data/outputs/maha_rcp26_h.bil")[[which(auc_h[,"maha"] >= 0.7)]]
maha_rcp45_h <- stack("./data/outputs/maha_rcp45_h.bil")[[which(auc_h[,"maha"] >= 0.7)]]
maha_rcp60_h <- stack("./data/outputs/maha_rcp60_h.bil")[[which(auc_h[,"maha"] >= 0.7)]]
maha_rcp85_h <- stack("./data/outputs/maha_rcp85_h.bil")[[which(auc_h[,"maha"] >= 0.7)]]

maha_h <- stack(maha_c_h, maha_rcp26_h, maha_rcp45_h, maha_rcp60_h, maha_rcp85_h)
maha_auc_h <- auc_h[which(auc_h[,"maha"] >= 0.7),"maha"]

#  Creating ANOVA Factors to be used at Uncertainty Evaluation
maha_method_h <- rep("maha", nlayers(maha_h))
maha_period_h <- rep(c("pres","fut"), each = nlayers(maha_h))

## maxent

#Selecting models with auc values ≥ 0.7
maxent_c_h <- stack("./data/outputs/maxent_c_h.bil")[[which(auc_h[,"maxent"] >= 0.7)]]
maxent_rcp26_h <- stack("./data/outputs/maxent_rcp26_h.bil")[[which(auc_h[,"maxent"] >= 0.7)]]
maxent_rcp45_h <- stack("./data/outputs/maxent_rcp45_h.bil")[[which(auc_h[,"maxent"] >= 0.7)]]
maxent_rcp60_h <- stack("./data/outputs/maxent_rcp60_h.bil")[[which(auc_h[,"maxent"] >= 0.7)]]
maxent_rcp85_h <- stack("./data/outputs/maxent_rcp85_h.bil")[[which(auc_h[,"maxent"] >= 0.7)]]

maxent_h <- stack(maxent_c_h, maxent_rcp26_h, maxent_rcp45_h, maxent_rcp60_h, maxent_rcp85_h)
maxent_auc_h <- auc_h[which(auc_h[,"maxent"] >= 0.7),"maxent"]

#  Creating ANOVA Factors to be used at Uncertainty Evaluation
maxent_method_h <- rep("maxent", nlayers(maxent_h))
maxent_period_h <- rep(c("pres","fut"), each = nlayers(maxent_h))

## SVM

#Selecting models with auc values ≥ 0.7
SVM_c_h <- stack("./data/outputs/SVM_c_h.bil")[[which(auc_h[,"SVM"] >= 0.7)]]
SVM_rcp26_h <- stack("./data/outputs/SVM_rcp26_h.bil")[[which(auc_h[,"SVM"] >= 0.7)]]
SVM_rcp45_h <- stack("./data/outputs/SVM_rcp45_h.bil")[[which(auc_h[,"SVM"] >= 0.7)]]
SVM_rcp60_h <- stack("./data/outputs/SVM_rcp60_h.bil")[[which(auc_h[,"SVM"] >= 0.7)]]
SVM_rcp85_h <- stack("./data/outputs/SVM_rcp85_h.bil")[[which(auc_h[,"SVM"] >= 0.7)]]

SVM_h <- stack(SVM_c_h, SVM_rcp26_h, SVM_rcp45_h, SVM_rcp60_h, SVM_rcp85_h)
SVM_auc_h <- auc_h[which(auc_h[,"SVM"] >= 0.7),"SVM"]

#  Creating ANOVA Factors to be used at Uncertainty Evaluation
SVM_method_h <- rep("SVM", nlayers(SVM_h))
SVM_period_h <- rep(c("pres","fut"), each = nlayers(SVM_h))

## GLM

#Selecting models with auc values ≥ 0.7
GLM_c_h <- stack("./data/outputs/GLM_c_h.bil")[[which(auc_h[,"GLM"] >= 0.7)]]
GLM_rcp26_h <- stack("./data/outputs/GLM_rcp26_h.bil")[[which(auc_h[,"GLM"] >= 0.7)]]
GLM_rcp45_h <- stack("./data/outputs/GLM_rcp45_h.bil")[[which(auc_h[,"GLM"] >= 0.7)]]
GLM_rcp60_h <- stack("./data/outputs/GLM_rcp60_h.bil")[[which(auc_h[,"GLM"] >= 0.7)]]
GLM_rcp85_h <- stack("./data/outputs/GLM_rcp85_h.bil")[[which(auc_h[,"GLM"] >= 0.7)]]

GLM_h <- stack(GLM_c_h, GLM_rcp26_h, GLM_rcp45_h, GLM_rcp60_h, GLM_rcp85_h)
GLM_auc_h <- auc_h[which(auc_h[,"GLM"] >= 0.7),"GLM"]

#  Creating ANOVA Factors to be used at Uncertainty Evaluation
GLM_method_h <- rep("GLM", nlayers(GLM_h))
GLM_period_h <- rep(c("pres","fut"), each = nlayers(GLM_h))


### Host plants
auc_p <- read.table("./data/outputs/AUC_plants.txt", h = T)

## bioclim

#Selecting models with auc values ≥ 0.7
bioclim_c_p <- stack("./data/outputs/bioclim_c_p.bil")[[which(auc_p[,"bioclim"] >= 0.7)]]
bioclim_rcp26_p <- stack("./data/outputs/bioclim_rcp26_p.bil")[[which(auc_p[,"bioclim"] >= 0.7)]]
bioclim_rcp45_p <- stack("./data/outputs/bioclim_rcp45_p.bil")[[which(auc_p[,"bioclim"] >= 0.7)]]
bioclim_rcp60_p <- stack("./data/outputs/bioclim_rcp60_p.bil")[[which(auc_p[,"bioclim"] >= 0.7)]]
bioclim_rcp85_p <- stack("./data/outputs/bioclim_rcp85_p.bil")[[which(auc_p[,"bioclim"] >= 0.7)]]

bioclim_p <- stack(bioclim_c_p, bioclim_rcp26_p, bioclim_rcp45_p, bioclim_rcp60_p, bioclim_rcp85_p)
bioclim_auc_p <- auc_p[which(auc_p[,"bioclim"] >= 0.7),"bioclim"]

#  Creating ANOVA Factors to be used at Uncertainty Evaluation
bioclim_method_p <- rep("bioclim", nlayers(bioclim_p))
bioclim_period_p <- rep(c("pres","fut"), each = nlayers(bioclim_p))

## gower

#Selecting models with auc values ≥ 0.7
gower_c_p <- stack("./data/outputs/gower_c_p.bil")[[which(auc_p[,"gower"] >= 0.7)]]
gower_rcp26_p <- stack("./data/outputs/gower_rcp26_p.bil")[[which(auc_p[,"gower"] >= 0.7)]]
gower_rcp45_p <- stack("./data/outputs/gower_rcp45_p.bil")[[which(auc_p[,"gower"] >= 0.7)]]
gower_rcp60_p <- stack("./data/outputs/gower_rcp60_p.bil")[[which(auc_p[,"gower"] >= 0.7)]]
gower_rcp85_p <- stack("./data/outputs/gower_rcp85_p.bil")[[which(auc_p[,"gower"] >= 0.7)]]

gower_p <- stack(gower_c_p, gower_rcp26_p, gower_rcp45_p, gower_rcp60_p, gower_rcp85_p)
gower_auc_p <- auc_p[which(auc_p[,"gower"] >= 0.7),"gower"]

#  Creating ANOVA Factors to be used at Uncertainty Evaluation
gower_method_p <- rep("gower", nlayers(gower_p))
gower_period_p <- rep(c("pres","fut"), each = nlayers(gower_p))

## maha

#Selecting models with auc values ≥ 0.7
maha_c_p <- stack("./data/outputs/maha_c_p.bil")[[which(auc_p[,"maha"] >= 0.7)]]
maha_rcp26_p <- stack("./data/outputs/maha_rcp26_p.bil")[[which(auc_p[,"maha"] >= 0.7)]]
maha_rcp45_p <- stack("./data/outputs/maha_rcp45_p.bil")[[which(auc_p[,"maha"] >= 0.7)]]
maha_rcp60_p <- stack("./data/outputs/maha_rcp60_p.bil")[[which(auc_p[,"maha"] >= 0.7)]]
maha_rcp85_p <- stack("./data/outputs/maha_rcp85_p.bil")[[which(auc_p[,"maha"] >= 0.7)]]

maha_p <- stack(maha_c_p, maha_rcp26_p, maha_rcp45_p, maha_rcp60_p, maha_rcp85_p)
maha_auc_p <- auc_p[which(auc_p[,"maha"] >= 0.7),"maha"]

#  Creating ANOVA Factors to be used at Uncertainty Evaluation
maha_method_p <- rep("maha", nlayers(maha_p))
maha_period_p <- rep(c("pres","fut"), each = nlayers(maha_p))

## maxent

#Selecting models with auc values ≥ 0.7
maxent_c_p <- stack("./data/outputs/maxent_c_p.bil")[[which(auc_p[,"maxent"] >= 0.7)]]
maxent_rcp26_p <- stack("./data/outputs/maxent_rcp26_p.bil")[[which(auc_p[,"maxent"] >= 0.7)]]
maxent_rcp45_p <- stack("./data/outputs/maxent_rcp45_p.bil")[[which(auc_p[,"maxent"] >= 0.7)]]
maxent_rcp60_p <- stack("./data/outputs/maxent_rcp60_p.bil")[[which(auc_p[,"maxent"] >= 0.7)]]
maxent_rcp85_p <- stack("./data/outputs/maxent_rcp85_p.bil")[[which(auc_p[,"maxent"] >= 0.7)]]

maxent_p <- stack(maxent_c_p, maxent_rcp26_p, maxent_rcp45_p, maxent_rcp60_p, maxent_rcp85_p)
maxent_auc_p <- auc_p[which(auc_p[,"maxent"] >= 0.7),"maxent"]

#  Creating ANOVA Factors to be used at Uncertainty Evaluation
maxent_method_p <- rep("maxent", nlayers(maxent_p))
maxent_period_p <- rep(c("pres","fut"), each = nlayers(maxent_p))

## SVM

#Selecting models with auc values ≥ 0.7
SVM_c_p <- stack("./data/outputs/SVM_c_p.bil")[[which(auc_p[,"SVM"] >= 0.7)]]
SVM_rcp26_p <- stack("./data/outputs/SVM_rcp26_p.bil")[[which(auc_p[,"SVM"] >= 0.7)]]
SVM_rcp45_p <- stack("./data/outputs/SVM_rcp45_p.bil")[[which(auc_p[,"SVM"] >= 0.7)]]
SVM_rcp60_p <- stack("./data/outputs/SVM_rcp60_p.bil")[[which(auc_p[,"SVM"] >= 0.7)]]
SVM_rcp85_p <- stack("./data/outputs/SVM_rcp85_p.bil")[[which(auc_p[,"SVM"] >= 0.7)]]

SVM_p <- stack(SVM_c_p, SVM_rcp26_p, SVM_rcp45_p, SVM_rcp60_p, SVM_rcp85_p)
SVM_auc_p <- auc_p[which(auc_p[,"SVM"] >= 0.7),"SVM"]

#  Creating ANOVA Factors to be used at Uncertainty Evaluation
SVM_method_p <- rep("SVM", nlayers(SVM_p))
SVM_period_p <- rep(c("pres","fut"), each = nlayers(SVM_p))

## GLM

#Selecting models with auc values ≥ 0.7
GLM_c_p <- stack("./data/outputs/GLM_c_p.bil")[[which(auc_p[,"GLM"] >= 0.7)]]
GLM_rcp26_p <- stack("./data/outputs/GLM_rcp26_p.bil")[[which(auc_p[,"GLM"] >= 0.7)]]
GLM_rcp45_p <- stack("./data/outputs/GLM_rcp45_p.bil")[[which(auc_p[,"GLM"] >= 0.7)]]
GLM_rcp60_p <- stack("./data/outputs/GLM_rcp60_p.bil")[[which(auc_p[,"GLM"] >= 0.7)]]
GLM_rcp85_p <- stack("./data/outputs/GLM_rcp85_p.bil")[[which(auc_p[,"GLM"] >= 0.7)]]

GLM_p <- stack(GLM_c_p, GLM_rcp26_p, GLM_rcp45_p, GLM_rcp60_p, GLM_rcp85_p)
GLM_auc_p <- auc_p[which(auc_p[,"GLM"] >= 0.7),"GLM"]

#  Creating ANOVA Factors to be used at Uncertainty Evaluation
GLM_method_p <- rep("GLM", nlayers(GLM_p))
GLM_period_p <- rep(c("pres","fut"), each = nlayers(GLM_p))


# 08. Standardize suitabilities (suit) #######################################################
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

# 09. Ensemble ####################################################################################

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


# 10. Uncertainty Evaluation ######################################################################

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





