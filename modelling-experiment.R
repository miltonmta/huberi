require(raster)
require(rgdal)
require(abind)
require(vegan)
require(maps)
require(mask)

# 01. read aogcms models####

#?? Here we will import and process only the worldclim varibles for current conditions (~1960-1990).We'll submmit them to a varimax selection procedure. Once, selected throug the loadings, we'll make use of the same varible number in each future model acroos the RCPs

# Current Model ----
nuala <- function (dir)
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

current <- nuala(dir = "./data/climatic_vars/current")
current [1:5, ]
nrow(current)

# RCPs Models ----

# based on the cluster analysis we'll import only the selected variables at each RCP scenario.

tinker_bell <- function(x)
{
  model_raw <- stack(list.files(x,  pattern = ".tif$", full.names = TRUE))
  e <- extent(-122, -18, -56, 14)
  model_e <- crop(model_raw, e)
  model_val <- getValues(model_e)
  coord_model <- xyFromCell(model_e, 1:ncell(model_e))
  model <- cbind(coord_model, model_val)
  model <- na.omit(model)
  # model <- rasterToPoints(model) # suggested by R. Hijimans at SO
  return(model)
}

# append function
apn <- function(...) abind(..., along = 3) # empty function for setting appending par to main function

# Models from RCP 26
x <- list.dirs("./data/climatic_vars/selected/26bi70/", full.names = TRUE)[-1]
model_list <- lapply(x, tinker_bell)
rcp_26 <- do.call("apn", model_list)

# Models from RCP 45
x <- list.dirs("./data/climatic_vars/selected/45bi70/", full.names = TRUE)[-1]
model_list <- lapply(x, tinker_bell)
rcp_45 <- do.call("apn", model_list)
rm(rcp_45_tinker_bell)

# Models from RCP 60
x <- list.dirs("./data/climatic_vars/selected/60bi70/", full.names = TRUE)[-1]
model_list <- lapply(x, tinker_bell)
rcp_60 <- do.call("apn", model_list)

# Models from RCP 85
x <- list.dirs("./data/climatic_vars/selected/85bi70/", full.names = TRUE)[-1]
model_list <- lapply(x, tinker_bell)
rcp_85 <- do.call("apn", model_list)


## Wolrdclim GCM code----
# At worldclim.com by the resolution of 2.5min we selected  only the variables that apper simultaneously in all the Representative Concetration Pathways projetions (RCP26, RCP45, RCP60, RCP80). The codes of the 11 GCMs utilized are: bc, cc, gs, hd, he, ip, mi, mr, mc, mg, no.
browseURL("http://www.worldclim.org/cmip5_2.5m")

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

model_names <- c("BCC-CSM1-1", "CCSM4", "GISS-EZ-R", "HadGEM2-AO", "HadGEM2-ES", "IPSL-CM5A-LR", "MIROC5", "MRI-CGCM3", "MIROC-ESM-CHEM", "MIROC-ESM", "NorESM1-M")# naming must be in the same directory reading order.

# 02. Varimax variable selection####
require(psych)

fa.parallel(current[ , -c(1:2)], fa = 'fa') #scree plot
# The estimated weights for the factor scores are probably incorrect.  Try a different factor extraction method.
# In factor.scores, the correlation matrix is singular, an approximation is used
# Parallel analysis suggests that the number of factors =  5  and the number of components =  NA 
# Warning messages:
#   1: In cor.smooth(R) : Matrix was not positive definite, smoothing was done
# 2: In cor.smooth(R) : Matrix was not positive definite, smoothing was done
# 3: In cor.smooth(R) : Matrix was not positive definite, smoothing was done
# 4: In cor.smooth(r) : Matrix was not positive definite, smoothing was done
# 5: In cor.smooth(r) : Matrix was not positive definite, smoothing was done
current_fa <- fa(current[ , -c(1:2)], nfactors = 5, rotate = 'varimax')
current_loadings <- loadings(current_fa)


# fa.parallel(clima_now_coord_val[,-c(1:2)], fa='fa') #screen plot
# clima_now_fa <- fa(clima_now_val[,-c(1:2)], nfactors= 5, rotate= 'varimax')
# clima_now_loadings <- loadings(clima_now_fa)
# Where is bio1, bio2, and bio19?
# Loadings:
#   MR1    MR3    MR2    MR4    MR5   
# bio3   0.235  0.819         0.253  0.144
# bio4  -0.215 -0.912        -0.327 -0.102
# bio5   0.976                0.133       
# bio6   0.744  0.596  0.141  0.244 -0.115
# bio7  -0.169 -0.857 -0.270 -0.245  0.135
# bio8   0.907  0.227         0.159  0.240
# bio9   0.732  0.552         0.247 -0.205
# bio10  0.980  0.139         0.136       
# bio11  0.783  0.563         0.261       
# bio12  0.241  0.378  0.526  0.714       
# bio13  0.300  0.425  0.169  0.830       
# bio14         0.178  0.934  0.216  0.106
# bio15         0.123 -0.748         0.195
# bio16  0.296  0.421  0.189  0.835       
# bio17         0.193  0.940  0.248       
# bio18  0.136  0.119  0.497  0.447  0.321
# 
# MR1   MR3   MR2   MR4   MR5
# SS loadings    4.824 3.884 3.009 2.687 0.326
# Proportion Var 0.301 0.243 0.188 0.168 0.020
# Cumulative Var 0.301 0.544 0.732 0.900 0.921

## Bioclimatic Variables Descriptions ----

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

### 03. Saving selected variables####

# Now we have selected the variables bio... from the current GMC (worldclim v 1.4 at 2,5"), we'll take the same variable number from the selected models of the each one of the four RCP scenarious  and save them  at "./data/climatic_vars" as a .grd file.

## current
# saving as table

write.table(current[,c("x", "y" )], "./data/climatic_vars/selected/current.txt", row.names = F, sep = " ")  
# alternatively:
# write.table(current[,c("x", "y")], "clima_AS.csv", row.names=F, sep=",")

# saving as raster
writeRaster(current$bio, "./data/climatic_vars/selected/bio_current.grd", format = "raster")
writeRaster(current$bio, "./data/climatic_vars/selected/bio_current.grd", format = "raster")
writeRaster(current$bio, "./data/climatic_vars/selected/bio_current.grd", format = "raster")
writeRaster(current$bio, "./data/climatic_vars/selected/bio_current.grd", format = "raster")
writeRaster(current$bio, "./data/climatic_vars/selected/bio_current.grd", format = "raster")

# reading selected variables raster files
bio <- raster("./data/climatic_vars/selected/bio_current.grd")
bio <- raster("./data/climatic_vars/selected/bio_current.grd")
bio <- raster("./data/climatic_vars/selected/bio_current.grd")
bio <- raster("./data/climatic_vars/selected/bio_current.grd")
bio <- raster("./data/climatic_vars/selected/bio_current.grd")

# Binding the variables
current_select <- stack(c(bio, bio, bio, bio, bio))
names(clima.AS) <- c("bio","bio", "bio", "bio", "bio")
plot(clima.AS)


## Creating objects for all RCPs with the selected variables 
#??####

# how to extract the selected variables simultaneously form all model at the four RCPs?
# since the variables have already been processed within the function tinker_bell, here we just have to get the selected ones from the arrays.

rcp26_select

rcp45_select

rcp60_select

rcp85_select


write.table(rcp26_select, "./data/climatic_vars/selected/rcp-26-select.txt", row.names = F, sep = "	")
write.table(rcp45_select, "./data/climatic_vars/selected/rcp-26-select.txt", row.names = F, sep = "	")
write.table(rcp60_select, "./data/climatic_vars/selected/rcp-26-select.txt", row.names = F, sep = "	")
write.table(rcp85_select, "./data/climatic_vars/selected/rcp-26-select.txt", row.names = F, sep = "	")


# 04. Occurrencies data####


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


## 05. Background Sampling####

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

# 06. Modelling Adequability Predictions####

require(raster)
require(dismo)
require(kernlab)
require(maps)
require(abind)
dyn.load('/Library/Java/JavaVirtualMachines/jdk1.8.0_131.jdk/Contents/Home/jre/lib/server/libjvm.dylib')
require(rJava)

# Problems loading Rjava ( necessary package for running Maxent)? 
# 1. Check if JDK is installed going to the path below. If not:
# browseURL("http://www.oracle.com/technetwork/java/javase/downloads/jdk8-downloads-2133151.html")
# browseURL ("https://stackoverflow.com/questions/30738974/rjava-load-error-in-rstudio-r-after-upgrading-to-osx-yosemite")


fairy_godmother <- function(occurrency_huberi = "...", occurrency_plants = "...", background_huberi = "...", background_plants = "...", cross_validation = ...)
{
  AOGCMs <- model_names
  
  ## loading occurrency and backgound data.
  occur_h <- read.table(occurrency_huberi, h = T)
  occur_p <- read.table(occurrency_plants, h = T)
  back_h <- read.table(background_huberi, h = T)
  back_p <- read.talbe(background_plants, h = T)
  
  ## Creating empty objects for saving results from all AOGCMs
  output_current <- output_rcp26 <- output_rcp45 <- output_rcp60 <- output_rcp85 <- NULL
  
  for (j in AOGCMs)
  {
    ### Loop AOGCMs ####
    
    ### Reading the climatic files
    # ???----
    # I know I must creat a loop here for iteranting througt the aogcm  of my selected models at all RCPs
    # how?
    AOGCM_CURRENT <- # j # only one model here. 
    AOGCM_RCP26 <- # j
    AOGCM_RCP45 <- # j
    AOGCM_RCP60 <- # j
    AOGCM_RCP85 <- # j
    
    ### Creating objects for saving partial results for each cross validation loop.
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
    bioclim_rcp26_p <- gower_rcp26_p <- maha_rcp26_p <- maxent_rcp26_p <- SVM_rcp26_p <- GLM_rcp26_p <- stack() # Should we run abind instead of stack??
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
      ### Bioclim ----
      
      ### huberi
      # ajusting models
      bioclim_model_h <- bioclim(training_h[training_h[,"pb"] == 1, -c(1:3)])
      
      # making predictions
      bioclim_c_h <- stack (bioclim_c_h, predict(object = bioclim_model_h, x = current_select))
      bioclim_rcp26_h <- stack (bioclim_rcp26_h, predict(object = bioclim_model_h, x = rcp26_select))
      bioclim_rcp45_h <- stack (bioclim_rcp45_h, predict(object = bioclim_model_h, x = rcp45_select))
      bioclim_rcp60_h <- stack (bioclim_rcp60_h, predict(object = bioclim_model_h, x = rcp60_select))
      bioclim_rcp85_h <- stack (bioclim_rcp85_h, predict(object = bioclim_model_h, x = rcp85_select))
      
      # Evaluating models
      bioclim_eval_h <- evaluate(p=testing_h[testing_h[, "pb"] == 1, -1], a = testing_h[testing_h[, "pb"] == 0, -1], model = bioclim_model_h)
      
      bioclim_e_h <- c(bioclim_e_h, bioclim_eval@auc)
      bioclim_t_h <- c(bioclim_t_h, threshold(bioclim_eval_h, "spec_sens"))
      
      ## host plants
      # ajusting models
      bioclim_model_p <- bioclim(training_p[training_p[,"pb"] == 1, -c(1:3)])
      
      # making predictions
      bioclim_c_p <- stack (bioclim_c_p, predict(object = bioclim_model_p, x = current_select))
      bioclim_rcp26_p <- stack (bioclim_rcp26_p, predict(object = bioclim_model_p, x = rcp26_select))
      bioclim_rcp45_p <- stack (bioclim_rcp45_p, predict(object = bioclim_model_p, x = rcp45_select))
      bioclim_rcp60_p <- stack (bioclim_rcp60_p, predict(object = bioclim_model_p, x = rcp60_select))
      bioclim_rcp85_p <- stack (bioclim_rcp85_p, predict(object = bioclim_model_p, x = rcp85_select))
      
      # Evaluating models
      bioclim_eval_p <- evaluate(p=testing_p[testing_p[, "pb"] == 1, -1], a = testing_p[testing_p[, "pb"] == 0, -1], model = bioclim_model_p)
      
      bioclim_e_p <- c(bioclim_e_p, bioclim_eval@auc)
      bioclim_t_p <- c(bioclim_t_p, threshold(bioclim_eval_p, "spec_sens"))
      
      ### Gower ----
      ## huberi
      # ajusting models
      # making predictions
      # Evaluating models
      
      ## host plants
      # ajusting models
      # making predictions
      # Evaluating models
      
      
      ### Maha ----
      ## huberi
      # ajusting models
      # making predictions
      # Evaluating models
      
      ## host plants
      # ajusting models
      # making predictions
      # Evaluating models
      
      
      ### Maxent ----
      ## huberi
      # ajusting models
      # making predictions
      # Evaluating models
      
      ## host plants
      # ajusting models
      # making predictions
      # Evaluating models
      
      
      ### SVM ----
      ## huberi
      # ajusting models
      # making predictions
      # Evaluating models
      
      ## host plants
      # ajusting models
      # making predictions
      # Evaluating models
      
      
      ### GLM ----
      ## huberi
      # ajusting models
      # making predictions
      # Evaluating models
      
      ## host plants
      # ajusting models
      # making predictions
      # Evaluating models
      
    }
    
  }
}

fairy_godmother (occurrency_huberi = "./data/ocurrencies/huberi-var.txt", occurrency_plants = "./data/ocurrencies/plants-var.txt", background_huberi = "./data/ocurrencies/Background-random-huberi.txt", background_plants = "./data/ocurrencies/Background-random-plants.txt")

# 07. Write predictions####

# 08. Selecting models####

# 09. Find/standardize suitabilities (suit) ####

# 10. Ensemble ####

# 11. Uncertainty Evaluation ####


