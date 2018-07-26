# Packages ######
# install.packages(c("raster", "rgdal", "abind", "vegan", "maps", "mask","pcych", "kernlab", "dismo", "rJava")

require(raster)
require(rgdal)
require(abind)
require(vegan)
require(maps)
require(mask)
require(kernlab)

# This script has an index table. If you are in RStudio go to Code > Show Document Outline (shift + command / clrt + o)

# The directories with the data utilized and the ones outputted here can be downloaded from the following OneDrive repositorium:
browseURL("https://1drv.ms/f/s!ApJZaitgpPr7gZtfS9n9mU9DDzXQMg")

# The models projected for 2070 (average for 2061-2080) were obtain at "worldclim.com" by the spatial resolution of 2.5min (0.04º or ≈ 4.4km). We selected the variables that appear simultaneously in all the representative concentration pathways scnarios (RCP26, RCP45, RCP60, RCP80). The codes of the 11 GCMs utilized are: bc, cc, gs, hd, he, ip, mi, mr, mc, mg, no.
browseURL("http://www.worldclim.org/cmip5_2.5m")


# 01. read aogcms models ############################################################################

#?? Here we will import and process only the worldclim variables for current conditions (~1960-1990), by the spatial resolution of 2.5min. We'll submmit them to a varimax selection procedure. Once, selected through the loadings scores, we'll make use of the same variable number in each future model across all the GCM models at the four RCPs scenarios for 2070.
browseURL("http://www.worldclim.org/current")

## Current Model ----
tinker_bell <- function (dir)                                                # 19 .bil files at the directory
{
  model_raw <- stack(list.files(dir,  pattern = ".bil$", full.names = TRUE)) # Object Type RasterStack (241.1KB) nrow 3600
  e <- extent(-122, -18, -56, 14) 
  model_e <- crop(model_raw, e)                                              # Object Type RasterBrick (303.9 MB) nrow 1680
  val <- getValues(model_e)                                                  # Object Type matrix (303.9 MB) nrow 4193280
  coord <- xyFromCell(model_e, 1:ncell(model_e))                             # Object Type matrix (64MB) nrow 4193280
  model <- cbind(coord, val)                                                 # Object Type matrix (671,8 MB) nrow 4193280
  model <- na.omit(model)                                                    # Object Type matrix (171,2 MB) nrow 912559
  return(model)
}

current <- tinker_bell(dir = "./data/climatic_vars/current")                  # Object Type matrix (171,2 MB) nrow 912559

## RCPs Models ----

# Wolrdclim GCM codes

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

#!#####
# for maintanining comparability... 

# RCP 26: CCSM4(CC),                   IPSL-CMSA-LR(IP), MIROC-ESM(MR) 
# RCP 45: CCSM4(CC),                   IPSL-CMSA-LR(IP), MIROC-ESM(MR) 
# RCP 60: CCSM4(CC), IPSL-CMSA-LR(IP),                   MIROC-ESM(MR)
# RCP 85: CCSM4(CC), IPSL-CMSA-LR(IP),                   MIROC-ESM(MR)


# naming must be in the exact order of the origin directory.
model_names <- c("CCSM4", "IPSL-CM5A-LR", "MIROC-ESM")


tooth_fairy <- function (x)
{
  directories <- list.dirs( x, full.names = TRUE)[-1]
  e <- extent(-122, -18, -56, 14)
  rcp <- NULL
  # models <- list()
  
  for (i in 1:length(directories))
  {
    models_raw <- stack(list.files(directories[i],pattern = ".tif$", full.names = TRUE))
    models_e <- crop( models_raw , e ) 
    val <- values (models_e)
    coord <- xyFromCell(models_e, 1:ncell(models_e))
    models <- cbind(coord, val)
    models <- na.omit(models)
    # models <- rasterToPoints(model) # suggested by R. Hijimans at SO
    rcp <- abind (rcp, models, along = 3)
  }
  
  return(rcp) 
}

rcp26_select <- tooth_fairy( x = "./data/climatic_vars/selected/26bi70/")
rcp45_select <- tooth_fairy( x = "./data/climatic_vars/selected/45bi70/")
rcp60_select <- tooth_fairy( x = "./data/climatic_vars/selected/60bi70/")
rcp85_select <- tooth_fairy( x = "./data/climatic_vars/selected/85bi70/")


# 02. Varimax variable selection #########################################################################
# install.packages(c("psych", "GPArotation"), dependencies = TRUE)
require(psych)
# require(GPArotation)

fa.parallel(current[ , -c(1:2)], fa = 'fa') #scree plot
current_fa <- fa(current[ , -c(1:2)], nfactors = 5, rotate = 'varimax')
loadings <- loadings(current_fa)

#!!####
# Matheus Ribeiro is running the factorial selection analysis with the "current" object to provide the loadings table. For now we'll use the following variables for adjusting the scrip:
#bio1, bio2, bio3, bio16, bio17

## Bioclimatic Variables Descriptions #######

#BIO1 = Annual Mean Temperature
#BIO2 = Mean Diurnal Range (Mean of monthly (max temp - min temp))
#BIO3 = Isothermality (BIO2/BIO7) (* 100)
#BIO4 = Temperature Seasonality (standard deviation *100)
#BIO5 = Max Temperature of Warmest Month
#BIO6 = Min Temperature of Coldest Month
#BIO7 = Temperature Annual Range (BIO5-BIO6)
#BIO8 = Mean Temperature of Wettest Quarter
#BIO9 = Mean Temperature of Driest Quarter
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

# 03. Saving selected variables #######################################################################
# Selected variables: bio1, bio2, bio3, bio16, bio17.

### current

## saving as table

write.table(current[,c("x", "y", "bio1", "bio2", "bio3", "bio16", "bio17" )], "./data/climatic_vars/selected/current-select.txt", row.names = F, sep = " ")  

## saving as raster
#??----

current <- rasterFromXYZ(current) # The object current was outputted from the function tinker_bell as a matrix. It needs to be a RasterBrick type of object so we could extract the selected variables from it.

writeRaster(current$bio1, "./data/climatic_vars/selected/bio01_current.grd", format = "raster")
# Error in current$bio1 : $ operator is invalid for atomic vectors
writeRaster(current$bio2, "./data/climatic_vars/selected/bio02_current.grd", format = "raster")
writeRaster(current$bio3, "./data/climatic_vars/selected/bio03_current.grd", format = "raster")
writeRaster(current$bio16, "./data/climatic_vars/selected/bio16_current.grd", format = "raster")
writeRaster(current$bio17, "./data/climatic_vars/selected/bio17_current.grd", format = "raster")


current_select <- stack(list.files("./data/climatic_vars/selected",  pattern = ".grd$", full.names = TRUE))
plot(current_select)


## RCPs


# saving array of aogcm models as table??
write.table(rcp26_select, "./data/climatic_vars/selected/rcp26-select.txt", row.names = F, sep = "	")
write.table(rcp45_select, "./data/climatic_vars/selected/rcp45-select.txt", row.names = F, sep = "	")
write.table(rcp60_select, "./data/climatic_vars/selected/rcp60-select.txt", row.names = F, sep = "	")
write.table(rcp85_select, "./data/climatic_vars/selected/rcp85-select.txt", row.names = F, sep = "	")

# saving array of aogcm models as raster
writeRaster(rcp26_select, "./data/climatic_vars/selected/rcp26-select.grd", format = "raster")
writeRaster(rcp45_select, "./data/climatic_vars/selected/rcp45-select.grd", format = "raster")
writeRaster(rcp60_select, "./data/climatic_vars/selected/rcp60-select.grd", format = "raster")
writeRaster(rcp85_select, "./data/climatic_vars/selected/rcp85-select.grd", format = "raster")


# 04. Occurrencies data ###################################################################################


huberi <- read.table("./data/ocurrencies/huberi.txt", h = T)
huberi [1:5, ]
huberi <- huberi[, -1]

host_plants <- read.table("./data/ocurrencies/host-plants-species.txt", h = T)
host_plants [1:5, ]
host_plants <- host_plants[, -1]

plot(current_select$bio1)
points(huberi[,"long"], huberi[,"lat"], pch = 20)
points(host_plants[,"long"], host_plants[,"lat"], pch = 20)


# extracting variables from ocurrencies data cells
huberi_cell <- cellFromXY(current_select, huberi)
duplicated(huberi_cell)
huberi_cell <- unique(huberi_cell)
huberi_var <- extract(current_select, huberi_cell)

host_plants_cell <- cellFromXY(current_select, host_plants)
duplicated(host_plantss_cell)
host_plants_cell <- unique(host_plants_cell)
host_plants_var <- extract(current_select, host_plants_cell)

write.table(huberi_var, "./data/ocurrencies/huberi-var.txt", row.names = F, sep = " ") 
write.table(host_plants_var, "./data/ocurrencies/host-plants-var.txt", row.names = F, sep = " ") 

## 05. Background Sampling ##############################################################################

## CURRENT
# creating the background with the Neotropic study area
neot_c <- extract(current_select, 1:ncell(current_select))
neot_c.coords <- xyFromCell(current_select, 1:ncell(current_select)) 
neot_c <- cbind(neot_c.coords, cells = 1:ncell(current_select), neot_c)
neot_c <- na.omit(neot_c)

# sampling with huberi
back_id_huberi <- sample(1:nrow(neot_c), nrow(huberi_var))
back_huberi_c <- neot_c[back_id_huberi, ]
points(back_huberi_c[, "x"], back_huberi_c[, "y"], pch = 20, col = 'red')

# sampling with host plants
back_id_plants <- sample(1:nrow(neot_c), nrow(host_plants_var))
back_plants_c <- neot_c[back_id_plants, ]
points(back_plants_c[, "x"], back_plants_c[, "y"], pch = 20, col = 'blue')

# saving background files
write.table(back_huberi, "./data/ocurrencies/Background-random-huberi.txt", row.names = F, sep = " ") 
write.table(back_plants, "./data/ocurrencies/Background-random-plants.txt", row.names = F, sep = " ")


rm(list = ls())

# 06. Modelling Adequability Predictions ##############################################################

# require(raster)
# require(rgdal)
# require(vegan)
# require(abind)
require(dismo)
require(kernlab)
require(maps)
dyn.load('/Library/Java/JavaVirtualMachines/jdk1.8.0_131.jdk/Contents/Home/jre/lib/server/libjvm.dylib')
require(rJava)

# MaxEnt is available as a standalone Java program. Dismo has a function 'maxent' that communicates with this program. To use it you must first download the program from http://www.cs.princeton.edu/~schapire/maxent/. Put the le 'maxent.jar' in the 'java' folder of the 'dismo' package. That is the folder returned by system.file("java", package="dismo"). You need MaxEnt version 3.3.3b or higher.

# 1. Check if JDK is installed going to the path below. If not:
# browseURL("http://www.oracle.com/technetwork/java/javase/downloads/jdk8-downloads-2133151.html")
# browseURL ("https://stackoverflow.com/questions/30738974/rjava-load-error-in-rstudio-r-after-upgrading-to-osx-yosemite")


# 2. checking if the jar file is present. If not, skip this bit
# 	jar <- paste(system.file(package="dismo"), "/java/maxent.jar", sep='')
#	 file.exists(jar)

# 3. Use dyn.load for ataching the libjvm.dylib file before running the `rJava` Package.


fairy_godmother <- function(occurrency_huberi = "...", occurrency_plants = "...", background_huberi = "...", background_plants = "...", cross_validation = ...)
{
  
  # Creating empty objects for considering the difefferent AOGCMs
  output_current <- output_rcp26 <- output_rcp45 <- output_rcp60 <- output_rcp85 <- NULL
  
  # model_names <- c("CCSM4", "IPSL-CM5A-LR", "MIROC-ESM")
  AOGCMs <- model_names
  
  for (j in AOGCMs)
  {
    ### Loop AOGCMs -----
    
    ### Reading the climatic files
    # ???----
    # I know I must creat a loop here for iteranting througt the aogcm  of my selected models at all RCPs.  How?
    
    AOGCM_CURRENT <- # j   # only one model here. 
    AOGCM_RCP26   <- # j   # Three models here.
    AOGCM_RCP45   <- # j   # Three models here.
    AOGCM_RCP60   <- # j   # Three models here.
    AOGCM_RCP85   <- # j   # Three models here.
    
    ### loading occurrency and backgound data.
    occur_h <- read.table(occurrency_huberi, h = T)
    occur_p <- read.table(occurrency_plants, h = T)
    back_h <- read.table(background_huberi, h = T)
    back_p <- read.talbe(background_plants, h = T)  
      
    ### Creating objects for saving partial results for each cross validation loop
    ## huberi
    bioclim_c_h <- gower_c_h <- maha_c_h <- maxent_c_h <- SVM_c_h <- GLM_c_h <- stack()
    bioclim_rcp26_h <- gower_rcp26_h <- maha_rcp26_h <- maxent_rcp26_h <- SVM_rcp26_h <- GLM_rcp26_h <- stack() # Should we run abind instead of stack??
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
      ### Loop Cross-validation ----
      
      ### creating trainning-testing subsets
      
      ## huberi
      sample_occur_h <- sample(1:nrow(occur_h), round(0.75 * nrow(occur_h)))
      sample_back_h  <- sample(1:nrow(back_h), round(0.75 * nrow(back_h)))
      training_h <- prepareData(x = current_select, p = occur_h[sample_occur_h,  1:2], b = back_h[sample_back_h,  1:2], xy = T)
      testing_h  <- prepareData(x = current_select, p = occur_h[-sample_occur_h, 1:2], b = back_h[-sample_back_h, 1:2], xy = T)
      
      ## host plants
      sample_occur_p <- sample(1:nrow(occur_p), round(0.75 * nrow(occur_p)))
      sample_back_p  <- sample(1:nrow(back_p), round(0.75 * nrow(back_p)))
      training_p <- prepareData(x = current_select, p = occur_p[sample_occur_p,  1:2], b = back_p[sample_back_p,  1:2], xy = T)
      testing_p  <- prepareData(x = current_select, p = occur_p[-sample_occur_p, 1:2], b = back_p[-sample_back_p, 1:2], xy = T)
      
      ### Bioclim -------------------
      
      ## huberi
      # ajusting models
      bioclim_model_h <- bioclim(training_h[training_h[,"pb"] == 1, -c(1:3)])
      
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
      
    }# closes the "i" - loop cross validation
    
    
    ## Writing predictions as raster  ################################################
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
    
  }# closes the "j" - loop AOGCMs
  
  
}# closes the function fairy_godmother

# Running my model ####
fairy_godmother (occurrency_huberi = "./data/ocurrencies/huberi-var.txt", occurrency_plants = "./data/ocurrencies/plants-var.txt", background_huberi = "./data/ocurrencies/Background-random-huberi.txt", background_plants = "./data/ocurrencies/Background-random-plants.txt")


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


################################ List of improvements to the scritp ###############################

# 1. Implement occurrency data filtering at the ambiental space!
# 1. Transform maps in frequencies instead of suitabilities.
# 2. Impelement multi cores for running several models simultaneously.




