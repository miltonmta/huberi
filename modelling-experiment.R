require(raster)
require(rgdal)
require(psych)
require(abind)
require(amap)
require(stats)
require(vegan)
require(dismo)
require(kernlab)
require(maps)
require(mask)
#install.packages("rJava)
dyn.load('/Library/Java/JavaVirtualMachines/jdk1.8.0_131.jdk/Contents/Home/jre/lib/server/libjvm.dylib')
require(rJava)

# Problems loading Rjava ( necessary package for running Maxent)? 
# 1. Check if JDK is installed going to the path below. If not:
# browseURL("http://www.oracle.com/technetwork/java/javase/downloads/jdk8-downloads-2133151.html")
# browseURL ("https://stackoverflow.com/questions/30738974/rjava-load-error-in-rstudio-r-after-upgrading-to-osx-yosemite")

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
writeRaster(ccsm.0k.ASr$bio, "./data/climatic_vars/selected/bio_current.grd", format = "raster")
writeRaster(ccsm.0k.ASr$bio, "./data/climatic_vars/selected/bio_current.grd", format = "raster")
writeRaster(ccsm.0k.ASr$bio, "./data/climatic_vars/selected/bio_current.grd", format = "raster")
writeRaster(ccsm.0k.ASr$bio, "./data/climatic_vars/selected/bio_current.grd", format = "raster")
writeRaster(ccsm.0k.ASr$bio, "./data/climatic_vars/selected/bio_current.grd", format = "raster")

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

host_plant <- read.table("./data/ocurrencies/host-plant-species.txt", h = T)
host_plant [1:5, ]
host_plant <- host_plant[, -1]

plot(current_select$bio1)
points(huberi[,"long"], huberi[,"lat"], pch = 20)
points(host_plant[,"long"], host_plant[,"lat"], pch = 20)


# extracting variables from ocurrencies data cells
huberi_cell <- cellFromXY(current_select, huberi)
duplicated(huberi_cell)
huberi_cell <- unique(huberi_cell)
huberi_var <- extract(current_select, huberi_cell)

host_plant_cell <- cellFromXY(current_select, host_plant)
duplicated(host_plant_cell)
host_plant_cell <- unique(host_plant_cell)
host_plant_var <- extract(current_select, host_plant_cell)


# saving as a table
write.table(huberi_var, "./data/ocurrencies/huberi-var.txt", row.names = F, sep = " ") 

write.table(host_plant_var, "./data/ocurrencies/host-plant-var.txt", row.names = F, sep = " ") 


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
back_id_plant <- sample(1:nrow(neot_c), nrow(host_plant_var))
back_plant_c <- neot_c[back_id_plant, ]
points(back_plant_c[, "x"], back_plant_c[, "y"], pch = 20, col = 'blue')

write.table(back_huberi, "./data/ocurrencies/Background-random-huberi-current.txt", row.names = F, sep = " ") 
write.table(back_plant, "./data/ocurrencies/Background-random-plant-current.txt", row.names = F, sep = " ") 

## RCP26
# creating the background with the Neotropic study area
neot_rcp26 <- extract(rcp26_select, 1:ncell(rcp26_select))
neot_rcp26.coords <- xyFromCell(rcp26_select, 1:ncell(rcp26_select)) 
neot_rcp26 <- cbind(neot_rcp26.coords, cells = 1:ncell(rcp26_select), neot_rcp26)
neot_rcp26 <- na.omit(neot_rcp26)

# sampling with huberi
back_id_huberi <- sample(1:nrow(neot_rcp26), nrow(huberi_var))
back_huberi_rcp26 <- neot_rcp26[back_id_huberi, ]
points(back_huberi_rcp26[, "x"], back_huberi_rcp26[, "y"], pch = 20, col = 'red')

# sampling with host plants
back_id_plant <- sample(1:nrow(neot_rcp26), nrow(host_plant_var))
back_plant_rcp26 <- neot_rcp26[back_id_plant, ]
points(back_plant_rcp26[, "x"], back_plant_rcp26[, "y"], pch = 20, col = 'blue')

write.table(back_huberi_rcp26, "./data/ocurrencies/Background-random-huberi-rcp26.txt", row.names = F, sep = " ") 
write.table(back_plant_rcp26, "./data/ocurrencies/Background-random-plant-rcp26.txt", row.names = F, sep = " ") 

## RCP45
# creating the background with the Neotropic study area
neot_rcp45 <- extract(rcp45_select, 1:ncell(rcp45_select))
neot_rcp45.coords <- xyFromCell(rcp45_select, 1:ncell(rcp45_select)) 
neot_rcp45 <- cbind(neot_rcp45.coords, cells = 1:ncell(rcp45_select), neot_rcp45)
neot_rcp45 <- na.omit(neot_rcp45)

# sampling with huberi
back_id_huberi <- sample(1:nrow(neot_rcp45), nrow(huberi_var))
back_huberi_rcp45 <- neot_rcp45[back_id_huberi, ]
points(back_huberi_rcp45[, "x"], back_huberi_rcp45[, "y"], pch = 20, col = 'red')

# sampling with host plants
back_id_plant <- sample(1:nrow(neot_rcp45), nrow(host_plant_var))
back_plant_rcp45 <- neot_rcp45[back_id_plant, ]
points(back_plant_rcp45[, "x"], back_plant_rcp45[, "y"], pch = 20, col = 'blue')

write.table(back_huberi_rcp45, "./data/ocurrencies/Background-random-huberi-rcp45.txt", row.names = F, sep = " ") 
write.table(back_plant_rcp45, "./data/ocurrencies/Background-random-plant-rcp45.txt", row.names = F, sep = " ") 

## RCP60
# creating the background with the Neotropic study area
neot_rcp60 <- extract(rcp60_select, 1:ncell(rcp60_select))
neot_rcp60.coords <- xyFromCell(rcp60_select, 1:ncell(rcp60_select)) 
neot_rcp60 <- cbind(neot_rcp60.coords, cells = 1:ncell(rcp60_select), neot_rcp60)
neot_rcp60 <- na.omit(neot_rcp60)

# sampling with huberi
back_id_huberi <- sample(1:nrow(neot_rcp60), nrow(huberi_var))
back_huberi_rcp60 <- neot_rcp60[back_id_huberi, ]
points(back_huberi_rcp60[, "x"], back_huberi_rcp60[, "y"], pch = 20, col = 'red')

# sampling with host plants
back_id_plant <- sample(1:nrow(neot_rcp60), nrow(host_plant_var))
back_plant_rcp60 <- neot_rcp60[back_id_plant, ]
points(back_plant_rcp60[, "x"], back_plant_rcp60[, "y"], pch = 20, col = 'blue')

write.table(back_huberi_rcp60, "./data/ocurrencies/Background-random-huberi-rcp60.txt", row.names = F, sep = " ") 
write.table(back_plant_rcp60, "./data/ocurrencies/Background-random-plant-rcp60.txt", row.names = F, sep = " ") 


## RCP85
# creating the background with the Neotropic study area
neot_rcp85 <- extract(rcp85_select, 1:ncell(rcp85_select))
neot_rcp85.coords <- xyFromCell(rcp85_select, 1:ncell(rcp85_select)) 
neot_rcp85 <- cbind(neot_rcp85.coords, cells = 1:ncell(rcp85_select), neot_rcp85)
neot_rcp85 <- na.omit(neot_rcp85)

# sampling with huberi
back_id_huberi <- sample(1:nrow(neot_rcp85), nrow(huberi_var))
back_huberi_rcp85 <- neot_rcp85[back_id_huberi, ]
points(back_huberi_rcp85[, "x"], back_huberi_rcp85[, "y"], pch = 20, col = 'red')

# sampling with host plants
back_id_plant <- sample(1:nrow(neot_rcp85), nrow(host_plant_var))
back_plant_rcp85 <- neot_rcp85[back_id_plant, ]
points(back_plant_rcp85[, "x"], back_plant_rcp85[, "y"], pch = 20, col = 'blue')

write.table(back_huberi_rcp85, "./data/ocurrencies/Background-random-huberi-rcp85.txt", row.names = F, sep = " ") 
write.table(back_plant_rcp85, "./data/ocurrencies/Background-random-plant-rcp85.txt", row.names = F, sep = " ") 


# 06. Modelling Adequability Predictions####

output_current <- NULL
output_rcp26 <- NULL
output_rcp45 <- NULL
output_rcp60 <- NULL
output_rcp85 <- NULL
model_names
for (j in model_names)
{
  aogcm_current_huberi <- paste("huberi_var", j, "back_huberi", sep = "")
  aogcm_current_plant <- paste("host_plant_var", j, "back_plant", sep = "")
  
  
}


# 07. Write predictions####

# 08. Selecting models####

# 09. Find/standardize suitabilities (suit) ####

# 10. Ensemble ####

# 11. Uncertainty Evaluation ####


