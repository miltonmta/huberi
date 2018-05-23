# " -----------------------------------------------
# packages####
#rm(list=ls())
require(raster)
require(maps)
#library(psych)
require(vegan)
library(dismo)
require(kernlab)
#install.packages(mask)
dyn.load('/Library/Java/JavaVirtualMachines/jdk1.8.0_131.jdk/Contents/Home/jre/lib/server/libjvm.dylib')
library(rJava)
#"--------------------------------------------------

# 01. Preparing for variables selection####
#Here we will import and prepare only the worldclim varibles for current conditions (~1960-1990) to proceed them for varimax selection. Once, selected throug the loadings, we'll make use of the same varible number in each future GCMs.
#variable rastering and stacking



now_bio_1 <- raster("./data/climatic_vars/vars_present/bio1.asc")
now_bio_2 <- raster("./data/climatic_vars/vars_present/bio2.asc")
now_bio_3 <- raster("./data/climatic_vars/vars_present/bio3.asc")
now_bio_4 <- raster("./data/climatic_vars/vars_present/bio4.asc")
now_bio_5 <- raster("./data/climatic_vars/vars_present/bio5.asc")
now_bio_6 <- raster("./data/climatic_vars/vars_present/bio6.asc")
now_bio_7 <- raster("./data/climatic_vars/vars_present/bio7.asc")
now_bio_8 <- raster("./data/climatic_vars/vars_present/bio8.asc")
now_bio_9 <- raster("./data/climatic_vars/vars_present/bio9.asc")
now_bio_10 <- raster("./data/climatic_vars/vars_present/bio10.asc")
now_bio_11 <- raster("./data/climatic_vars/vars_present/bio11.asc")
now_bio_12 <- raster("./data/climatic_vars/vars_present/bio12.asc")
now_bio_13 <- raster("./data/climatic_vars/vars_present/bio13.asc")
now_bio_14 <- raster("./data/climatic_vars/vars_present/bio14.asc")
now_bio_15 <- raster("./data/climatic_vars/vars_present/bio15.asc")
now_bio_16 <- raster("./data/climatic_vars/vars_present/bio16.asc")
now_bio_17 <- raster("./data/climatic_vars/vars_present/bio17.asc")
now_bio_18 <- raster("./data/climatic_vars/vars_present/bio18.asc")
now_bio_19 <- raster("./data/climatic_vars/vars_present/bio19.asc")


clima_now <- stack(now_bio_1, now_bio_2, now_bio_3, now_bio_4, now_bio_5, now_bio_6, now_bio_7, now_bio_8, now_bio_9, now_bio_10, now_bio_11, now_bio_12, now_bio_13, now_bio_14, now_bio_15, now_bio_16, now_bio_17, now_bio_18, h=T) 

# Here we delimit the working extent for South America (sa) so we do not analyse the whole world extent of the original worldclim variables
e <- extent(-100,-35,-60,25)
clima_now_sa <- crop(clima_now, e)

#extracting raster values
clima_now_val <- values(clima_now_sa)
clima_now_val[1:5,] # bio 19 missing?
nrow(clima_now_val)

coord_sa <- xyFromCell(clima_now_sa, 1:ncell(clima_now_sa))
coord_sa[1:5,]
nrow(coord_sa)

clima_now_coord_val <- cbind(coord_sa, clima_now_val)
clima_now_coord_val[1:5,]
nrow(clima_now_coord_val)
clima_now_coord_val <- na.omit(clima_now_coord_val)
nrow(clima_now_coord_val)

# 02. Varimax variable selection####
library(psych)

fa.parallel(clima_now_coord_val[,-c(1:2)], fa='fa') #screen plot
clima_now_fa <- fa(clima_now_val[,-c(1:2)], nfactors= 5, rotate= 'varimax')
clima_now_loadings <- loadings(clima_now_fa)
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

# 03. Saving selected variables####

# Now we have selected the variables bio... from the current GMC (worldclim v 1.4 at 2,5"), we'll take the same variable number from all 17 future projections at RCP +8.5 and save them (current + 17 future) at "./data/climatic_vars" as a .grd file.

# 04. Background Sampling####

# 05. Occurrences data####

# 06. Modelling Adequability Predictions####

# 07. Write predictions####

# 08. Selecting models####

# 09. Find/standardize suitabilities (suit) ####

# 10. Ensemble ####

# 11. Uncertainty Evaluation ####


