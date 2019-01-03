file.edit("readme.R")
citation(package = "kernlab", lib.loc = NULL)
RStudio.Version()
R.Version()
# ****  Loading functions                            ----
source("./Auxiliary_functions.R")
source("./Multiple_ENMs.R")
source("./ensemble.R")
source("./predict_enfa.R")

# file.edit(c("./Auxiliary_functions.R", "./Multiple_ENMs.R"))

# ****  Loading packages                             ----
load_pak(c("tidyverse", "raster", "rgdal", "abind", "spThin", "dismo", "kernlab", "vegan", "maps", "psych", "rJava", "dendextend", "beepr", "data.table", "adehabitatHS"))

# ***************************************************************************************
## 01.  Read aogcms models                           ----

##### 1.1.............. current                                     
current_list <- read_current(dir = "./data/climatic_vars/current/")

current_mat <- current_list[["matrix"]] # required at Exploratoty Factor Analysis
current     <- current_list[["raster"]] # required for salving selected variables and creating "back"
beep(2)

##### 1.2.............. rcp
#. Following the order of the models in the directory:
# 1 == "CCS4"
# 2 == "IPSL-CMSA-LR"
# 3 == "MIROC-ESM"
# 
rcp26_list <- read_rcp( x = "./data/climatic_vars/selected/26bi70/")
rcp45_list <- read_rcp( x = "./data/climatic_vars/selected/45bi70/")
rcp60_list <- read_rcp( x = "./data/climatic_vars/selected/60bi70/")
rcp85_list <- read_rcp( x = "./data/climatic_vars/selected/85bi70/")


rcp26 <- rcp26_list[["rasters"]]
rcp26 <- stack(rcp26[[1]], rcp26[[2]], rcp26[[3]]) # this stack is adding .1, .2, .3 to the variables as models id. It needs to be removed when running the model. See line 189 from Multiple_ENMs.R
rcp45 <- rcp45_list[["rasters"]]
rcp45 <- stack(rcp45[[1]], rcp45[[2]], rcp45[[3]])

rcp60 <- rcp60_list[["rasters"]]
rcp60 <- stack(rcp60[[1]], rcp60[[2]], rcp60[[3]])

rcp85 <- rcp85_list[["rasters"]]
rcp85 <- stack(rcp85[[1]], rcp85[[2]], rcp85[[3]])

beep(2)
rm(rcp26_list, rcp45_list, rcp60_list, rcp85_list)

# ***************************************************************************************
## 02.  Variable selection                           ----

### By exploratory factor analysis

fa.parallel(current_mat[ , -c(1:2)], fa = 'fa') #scree plot
current_fa <- fa(current_mat[ , -c(1:2)], nfactors = 5, rotate = 'varimax')
loadings <- loadings(current_fa)
write.table(loadings, "./data/climatic_vars/selected/varimax_loadings.txt")
beep(2)

### Selected variables 
# bio02, bio03, bio10, bio14, bio16.

# ***************************************************************************************
## 03.  Saving selected variables                    ----

###................ CURRENT

writeRaster (current[[c("bio02", "bio03", "bio10", "bio14", "bio16")]], filename = "./data/climatic_vars/selected/current/current.grd", format = "raster", bylayer = F)

current_select <- stack("./data/climatic_vars/selected/current/current.grd")

beep(2)

###................ RCPs

variables <- c("bio02.1", "bio03.1", "bio10.1", "bio14.1", "bio16.1", 
               "bio02.2", "bio03.2", "bio10.2", "bio14.2", "bio16.2", 
               "bio02.3", "bio03.3", "bio10.3", "bio14.3", "bio16.3")


for (i in 1:length(variables))
{
  writeRaster(rcp26[[i]], filename = paste0("./data/climatic_vars/selected/rcp26/rcp26-", i, ".grd"), format = "raster")
  writeRaster(rcp45[[i]], filename = paste0("./data/climatic_vars/selected/rcp45/rcp45-", i, ".grd"), format = "raster")
  writeRaster(rcp60[[i]], filename = paste0("./data/climatic_vars/selected/rcp60/rcp60-", i, ".grd"), format = "raster")
  writeRaster(rcp85[[i]], filename = paste0("./data/climatic_vars/selected/rcp85/rcp85-", i, ".grd"), format = "raster")
}

beep(2)
rm(variables)

# ***************************************************************************************
## 04.  Occurrences data                             ----

###.............. Reading data
# the names of columms in the file raw_data must be: "SPEC", "LONG", "LAT".
occur_raw <- read.table("./data/occurrences/raw_data.txt", h = T)
 
###.............. Filtering occurrences in the geographical space
sp <- gsub("C[1-9]","", occur_raw$SPEC)
sp_names <- unique(sp)
occur_thinned <- NULL
for(i in 1:length(sp_names))
{
  sp <- occur_raw[occur_raw[, 1] == sp_names[i], ]
  occur <- thin(sp, thin.par = 20, reps = 100, max.files = 1, locs.thinned.list.return = T)
  occur_thinned <- rbind(occur_thinned, occur)
}
beep(2)

write.table(occur_thinned, "./data/occurrences/occur_thinned.txt", sep = ";", row.names = FALSE)


###.............. Extrancting bio variables based on the ocurrence cells
for(i in 1:length(sp_names))
{
  var <- create_var(occur_thinned[occur_thinned[, 1] == sp_names[i], ], sp_names[i])
}
beep(2)


# ***************************************************************************************
## 05.  Background Sampling                          ----

###.............. Background files 
# Creating and saving the object "back" for each studied species

var_files <- list.files("./data/occurrences/", pattern = "var", full.names = TRUE)
for(i in 1:length(var_files))
{
  var_file <- read.table(var_files[i], h = T, sep = ";")
  create_back(var_file, sp_names[i])
}
beep(2)

###.............. Plotting occurrences and background
back <- list.files("./data/occurrences/", pattern = "back", full.names = TRUE)
back <- bind_rows(lapply(back, fread))
back <- as.data.frame(back)
beep(2)

sp_names
plot(current_select$bio02)
points(occur_thinned[!grepl("Lithurgus_huberi", occur_thinned$SPEC), ][, -1], pch = "*", col = "red")
points(occur_thinned[occur_thinned[, 1] == "Lithurgus_huberi", ][,-1], pch = "*", col = "blue")
points(back[, "x"], back[, "y"], pch = "*", col = 'magenta')


# ***************************************************************************************
## 06.  Creating trainning-testing subsets           ----
###.....................................

#  For variation control, the occurrence that runs in all experiments (from XP1 to XP2),  will be modeled with the same random subsets. (see cross_validation loop in Multiple_ENMs)

occur <- read.table("./data/occurrences/var_Lithurgus_huberi.txt",  sep = ";", h = T)
back  <- read.table("./data/occurrences/back_Lithurgus_huberi.txt", sep = ";", h = T)
cross_validation <- 20
for (i in 1:cross_validation)
{
  sample_occur <- sample(1:nrow(occur), round(0.75 * nrow(occur), 0))
  trainning <- prepareData(x = current, 
                           p = occur[sample_occur,  1:2], 
                           b = back[sample_occur,  1:2]) 
  testing   <- prepareData(x = current, 
                           p = occur[-sample_occur, 1:2], 
                           b = back[-sample_occur, 1:2])
  
  write.table(trainning, paste0("./data/occurrences/subsets/trainning", i, ".txt"), sep = ";")
  write.table(testing,   paste0("./data/occurrences/subsets/testing",   i, ".txt"), sep = ";")
}
rm(occur, back)

# ***************************************************************************************

## 07.  XP1                                          ----
###.....................................
rm(list = ls())
source("./Multiple_ENMs.R")
# Variables	 Abiotic ( 5 vars )
# Input	     9 sps :: bee + 7 plants + "resource" * (summed occurs of all plants)
# Output     9 inputs predictions  + 4 predictive methods * 4 rcps * 3 aogcms
# Ensembles  117 * (9 present + 12 future (4 rcps * 3 aogcms)) 
# NEW VARS   Biotic Predictor Variables PA_SEP, PA_STK, SUIT_SEP, SUIT_SKT, resourceSEP, resourceSUIT.

# We've splitted XP1 in two for running a specifif set of variables (e.g SOIL) for just one group of the input occurrences.

## 07.a ..... Loading species names                  ---- 

## ... bee
occur_thinned <- read.table("./data/occurrences/occur_thinned.txt", sep = ";", h = T)
sp <- occur_thinned[occur_thinned[, 1] == "Lithurgus_huberi", ] # keepinng only the bee specie
sp <- gsub("C[1-9]","", sp$SPEC)
sp_name <- unique(sp)
sp_name

## ... 7 plants species
occur_thinned <- read.table("./data/occurrences/occur_thinned.txt", sep = ";", h = T) 
sp <- occur_thinned[!grepl("Lithurgus_huberi", occur_thinned$SPEC), ] # removing the bee species
sp <- gsub("C[1-9]","", sp$SPEC)
sp_names <- unique(sp)
sp_names

## ... All species
occur_thinned <- read.table("./data/occurrences/occur_thinned.txt", sep = ";", h = T) 
sp <- gsub("C[1-9]","", occur_thinned$SPEC)
ALLsp_names <- unique(sp)
ALLsp_names

rm(sp, occur_thinned)

## 07.b ..... XP1.1 - bee                            -----
#   .....................................................................................
sp_name
for (i in 1:length(sp_name))
{
  ###.............. Running the modedelling experiment
  result <- multiple_ENMs(occurrence       = paste0("./data/occurrences/var_",  sp_name[i], ".txt"),
                          background       = paste0("./data/occurrences/back_", sp_name[i], ".txt"),
                          biovar_current   = "./data/climatic_vars/selected/current/",
                          biovar_rcp26     = "./data/climatic_vars/selected/rcp26/",
                          biovar_rcp45     = "./data/climatic_vars/selected/rcp45/",
                          biovar_rcp60     = "./data/climatic_vars/selected/rcp60/",
                          biovar_rcp85     = "./data/climatic_vars/selected/rcp85/",
                          newvar_current   = 0,
                          newvar_rcp26     = 0,
                          newvar_rcp45     = 0,
                          newvar_rcp60     = 0,
                          newvar_rcp85     = 0,
                          trainning        = "./data/occurrences/subsets/trainning",
                          testing          = "./data/occurrences/subsets/testing",
                          AOGCMs           = c(1, 2, 3),
                          Pout             = paste0("./data/outputs/XP1/Pout/xp1_", sp_name[i], "_"),
                          cross_validation = 10)
  
  ###.............. Saving evaluation data
  write.table(result[["TPR_c"]],        paste0("./data/outputs/XP1/", sp_name[i], "_TPR_current.txt"), sep = "\t", row.names = F)
  write.table(result[["Threshold_c"]],  paste0("./data/outputs/XP1/", sp_name[i], "_t_current.txt"),   sep = "\t", row.names = F)
  write.table(result[["Pred_area_c"]],  paste0("./data/outputs/XP1/", sp_name[i], "_d_current.txt"),   sep = "\t", row.names = F)
  write.table(result[["Threshold_sens"]],  paste0("./data/outputs/XP1/", sp_name[i], "_t_sens_current.txt"),   sep = "\t", row.names = F)
  write.table(result[["AUC"]],  paste0("./data/outputs/XP1/", sp_name[i], "_AUC_current.txt"),   sep = "\t", row.names = F)
  
  rm(result)
  gc()
}
beep(8)

## 07.c ..... XP1.2 - 7 plants                       ----
#   .....................................................................................
# here we could include the soil vars.
sp_names
for (i in 1:length(sp_names))
{
  ###.............. Running the modedelling experiment
  result <- multiple_ENMs(occurrence       = paste0("./data/occurrences/var_",  sp_names[i], ".txt"),
                          background       = paste0("./data/occurrences/back_", sp_names[i], ".txt"),
                          biovar_current   = "./data/climatic_vars/selected/current/",
                          biovar_rcp26     = "./data/climatic_vars/selected/rcp26/",
                          biovar_rcp45     = "./data/climatic_vars/selected/rcp45/",
                          biovar_rcp60     = "./data/climatic_vars/selected/rcp60/",
                          biovar_rcp85     = "./data/climatic_vars/selected/rcp85/",
                          newvar_current   = 0,
                          newvar_rcp26     = 0,
                          newvar_rcp45     = 0,
                          newvar_rcp60     = 0,
                          newvar_rcp85     = 0,
                          trainning        = 0,
                          testing          = 0,
                          AOGCMs           = c(1, 2, 3),
                          Pout             = paste0("./data/outputs/XP1/Pout/xp1_", sp_names[i], "_"),
                          cross_validation = 10)
  
  ###.............. Saving evaluation data
  write.table(result[["TPR_c"]],        paste0("./data/outputs/XP1/", sp_names[i], "_TPR_current.txt"), sep = "\t", row.names = F)
  write.table(result[["Threshold_c"]],  paste0("./data/outputs/XP1/", sp_names[i], "_t_current.txt"),   sep = "\t", row.names = F)
  write.table(result[["Pred_area_c"]],  paste0("./data/outputs/XP1/", sp_names[i], "_d_current.txt"),   sep = "\t", row.names = F)
  write.table(result[["Threshold_sens"]],  paste0("./data/outputs/XP1/", sp_name[i], "_t_sens_current.txt"),   sep = "\t", row.names = F)
  write.table(result[["AUC"]],  paste0("./data/outputs/XP1/", sp_name[i], "_AUC_current.txt"),   sep = "\t", row.names = F)
  
  rm(result)
  gc()
}
beep(8)


## 07.e ..... New Predictor variables (XP2)          ----
#   .....................................................................................

biovar_current <- "./data/climatic_vars/selected/current/"
current <- stack(list.files(biovar_current,  pattern = ".grd$", full.names = TRUE))
current <- na.omit(current)
coords <- na.omit(cbind(xyFromCell(current, 1:ncell(current)), values(current)))[,1:2]
rm(biovar_current, current)

occur_thinned <- read.table("./data/occurrences/occur_thinned.txt", sep = ";", h = T)

sp_names
AOGCMs  <-  c(1, 2, 3)
SPid <- paste0("SP", 1:7)
SPid_stk <- "stk"

##      .......... sep_pa                            ----
## .......... generating vars .............
PApres    <- NULL
PArcp26.1 <- PArcp26.2 <- PArcp26.3 <- NULL
PArcp45.1 <- PArcp45.2 <- PArcp45.3 <- NULL
PArcp60.1 <- PArcp60.2 <- PArcp60.3 <- NULL
PArcp85.1 <- PArcp85.2 <- PArcp85.3 <- NULL
for (i in 1:length(sp_names))
{
  PA_df <- read.table(paste0("./data/outputs/XP1/ensembles/", sp_names[i],"_ENSEMBLES.txt"), sep = "\t", h = T)
  ensemble <- rasterFromXYZ(PA_df[, c(1:3)])
  coords.sp <- occur_thinned[occur_thinned[, 1] == sp_names[i], ][-1]
  thr <- quantile(na.omit(extract(ensemble, coords.sp)), 0.05) # get the no-omission thr from present
  
  PApres    <- cbind( PApres,    ifelse((PA_df[,  3]) >= thr, 1, 0))
  
  PArcp26.1 <- cbind( PArcp26.1, ifelse((PA_df[,  8]) >= thr, 1, 0))
  PArcp26.2 <- cbind( PArcp26.2, ifelse((PA_df[,  9]) >= thr, 1, 0))
  PArcp26.3 <- cbind( PArcp26.3, ifelse((PA_df[, 10]) >= thr, 1, 0))
  
  PArcp45.1 <- cbind( PArcp45.1, ifelse((PA_df[, 11]) >= thr, 1, 0))
  PArcp45.2 <- cbind( PArcp45.2, ifelse((PA_df[, 12]) >= thr, 1, 0))
  PArcp45.3 <- cbind( PArcp45.3, ifelse((PA_df[, 13]) >= thr, 1, 0))
  
  PArcp60.1 <- cbind( PArcp60.1, ifelse((PA_df[, 14]) >= thr, 1, 0))
  PArcp60.2 <- cbind( PArcp60.2, ifelse((PA_df[, 15]) >= thr, 1, 0))
  PArcp60.3 <- cbind( PArcp60.3, ifelse((PA_df[, 16]) >= thr, 1, 0))
  
  PArcp85.1 <- cbind( PArcp85.1, ifelse((PA_df[, 17]) >= thr, 1, 0))
  PArcp85.2 <- cbind( PArcp85.2, ifelse((PA_df[, 18]) >= thr, 1, 0))
  PArcp85.3 <- cbind( PArcp85.3, ifelse((PA_df[, 19]) >= thr, 1, 0))
}

PArcp26 <- list(PArcp26.1, PArcp26.2, PArcp26.3)
PArcp45 <- list(PArcp45.1, PArcp45.2, PArcp45.3)
PArcp60 <- list(PArcp60.1, PArcp60.2, PArcp60.3)
PArcp85 <- list(PArcp85.1, PArcp85.2, PArcp85.3)

rm(PArcp26.1, PArcp26.2, PArcp26.3, 
   PArcp45.1, PArcp45.2, PArcp45.3, 
   PArcp60.1, PArcp60.2, PArcp60.3,
   PArcp85.1, PArcp85.2, PArcp85.3)

## ............ saving vars ...............
PA_c  <- rasterFromXYZ(cbind(coords, PApres))
names(PA_c) <- SPid
writeRaster(PA_c,  "./data/outputs/predictors/XP2/sep_PA_c.tif",  format = "GTiff", overwrite = TRUE)

for (i in AOGCMs)
{
  PA_26 <- rasterFromXYZ(cbind(coords, PArcp26[[i]]))
  PA_45 <- rasterFromXYZ(cbind(coords, PArcp45[[i]]))
  PA_60 <- rasterFromXYZ(cbind(coords, PArcp60[[i]]))
  PA_85 <- rasterFromXYZ(cbind(coords, PArcp85[[i]]))
  
  names(PA_26) <- names(PA_45) <- names(PA_60) <- names(PA_85) <- SPid
  
  writeRaster(PA_26, paste0("./data/outputs/predictors/XP2/rcp26/sep_PA_26_", i,".tif"), format = "GTiff", overwrite = TRUE)
  writeRaster(PA_45, paste0("./data/outputs/predictors/XP2/rcp45/sep_PA_45_", i,".tif"), format = "GTiff", overwrite = TRUE)
  writeRaster(PA_60, paste0("./data/outputs/predictors/XP2/rcp60/sep_PA_60_", i,".tif"), format = "GTiff", overwrite = TRUE) 
  writeRaster(PA_85, paste0("./data/outputs/predictors/XP2/rcp85/sep_PA_85_", i,".tif"), format = "GTiff", overwrite = TRUE)
}

##      .......... sep_suit                          ----
## .......... generating vars .............
SUITpres    <- NULL
SUITrcp26.1 <- SUITrcp26.2 <- SUITrcp26.3 <- NULL
SUITrcp45.1 <- SUITrcp45.2 <- SUITrcp45.3 <- NULL
SUITrcp60.1 <- SUITrcp60.2 <- SUITrcp60.3 <- NULL
SUITrcp85.1 <- SUITrcp85.2 <- SUITrcp85.3 <- NULL

for (i in 1:length(sp_names))
{
  suit_df <- read.table(paste0("./data/outputs/XP1/ensembles/", sp_names[i],"_ENSEMBLES.txt"), sep = "\t", h = T)
  SUITpres  <- cbind( SUITpres,  suit_df[, 3]) 
  
  SUITrcp26.1 <- cbind( SUITrcp26.1, suit_df[,  8]) 
  SUITrcp26.2 <- cbind( SUITrcp26.2, suit_df[,  9])
  SUITrcp26.3 <- cbind( SUITrcp26.3, suit_df[, 10]) 
  
  SUITrcp45.1 <- cbind( SUITrcp45.1, suit_df[, 11])
  SUITrcp45.2 <- cbind( SUITrcp45.2, suit_df[, 12])
  SUITrcp45.3 <- cbind( SUITrcp45.3, suit_df[, 13])
  
  SUITrcp60.1 <- cbind( SUITrcp60.1, suit_df[, 14])    
  SUITrcp60.2 <- cbind( SUITrcp60.2, suit_df[, 15]) 
  SUITrcp60.3 <- cbind( SUITrcp60.3, suit_df[, 16]) 
  
  SUITrcp85.1 <- cbind( SUITrcp85.1, suit_df[, 17])
  SUITrcp85.2 <- cbind( SUITrcp85.2, suit_df[, 18])
  SUITrcp85.3 <- cbind( SUITrcp85.3, suit_df[, 19])
}

SUITrcp26 <- list(SUITrcp26.1, SUITrcp26.2, SUITrcp26.3)
SUITrcp45 <- list(SUITrcp45.1, SUITrcp45.2, SUITrcp45.3)
SUITrcp60 <- list(SUITrcp60.1, SUITrcp60.2, SUITrcp60.3)
SUITrcp85 <- list(SUITrcp85.1, SUITrcp85.2, SUITrcp85.3)

rm(SUITrcp26.1, SUITrcp26.2, SUITrcp26.3,
   SUITrcp45.1, SUITrcp45.2, SUITrcp45.3,
   SUITrcp60.1, SUITrcp60.2, SUITrcp60.3,
   SUITrcp85.1, SUITrcp85.2, SUITrcp85.3)

## ............ saving vars ...............
SUIT_c  <- rasterFromXYZ(cbind(coords, SUITpres))
names(SUIT_c) <- SPid
writeRaster(SUIT_c,  "./data/outputs/predictors/XP3/sep_suit_c.tif",  format = "GTiff", overwrite = TRUE)

for (i in AOGCMs)
{
  SUIT_26 <- rasterFromXYZ(cbind(coords, SUITrcp26[[i]]))
  SUIT_45 <- rasterFromXYZ(cbind(coords, SUITrcp45[[i]]))
  SUIT_60 <- rasterFromXYZ(cbind(coords, SUITrcp60[[i]]))
  SUIT_85 <- rasterFromXYZ(cbind(coords, SUITrcp85[[i]]))
  
  names(SUIT_26) <- names(SUIT_45) <- names(SUIT_60) <- names(SUIT_85) <- SPid
  
  writeRaster(SUIT_26, paste0("./data/outputs/predictors/XP3/rcp26/sep_suit_26_", i,".tif"), format = "GTiff", overwrite = TRUE)
  writeRaster(SUIT_45, paste0("./data/outputs/predictors/XP3/rcp45/sep_suit_45_", i,".tif"), format = "GTiff", overwrite = TRUE)
  writeRaster(SUIT_60, paste0("./data/outputs/predictors/XP3/rcp60/sep_suit_60_", i,".tif"), format = "GTiff", overwrite = TRUE)
  writeRaster(SUIT_85, paste0("./data/outputs/predictors/XP3/rcp85/sep_suit_85_", i,".tif"), format = "GTiff", overwrite = TRUE)
}

##      .......... stk_pa                            ----

pres_STK  <- rowSums(PApres)
PA_c_STK  <- rasterFromXYZ(cbind(coords, pres_STK))
names(PA_c_STK) <- SPid_stk
writeRaster(PA_c_STK,  "./data/outputs/predictors/XP4/stk_PA_c.tif",  format = "GTiff", overwrite = TRUE)

for (i in AOGCMs)
{
  rcp26_STK <- rowSums(PArcp26[[i]]) 
  rcp45_STK <- rowSums(PArcp45[[i]]) 
  rcp60_STK <- rowSums(PArcp60[[i]]) 
  rcp85_STK <- rowSums(PArcp85[[i]])
  
  PA_26_STK <- rasterFromXYZ(cbind(coords, rcp26_STK))
  PA_45_STK <- rasterFromXYZ(cbind(coords, rcp45_STK))
  PA_60_STK <- rasterFromXYZ(cbind(coords, rcp60_STK))
  PA_85_STK <- rasterFromXYZ(cbind(coords, rcp85_STK))
  
  names(PA_26_STK) <- names(PA_45_STK) <- names(PA_60_STK) <- names(PA_85_STK) <- SPid_stk
  
  writeRaster(PA_26_STK, paste0("./data/outputs/predictors/XP4/rcp26/stk_PA_26_", i,".tif"), format = "GTiff", overwrite = TRUE)
  writeRaster(PA_45_STK, paste0("./data/outputs/predictors/XP4/rcp45/stk_PA_45_", i,".tif"), format = "GTiff", overwrite = TRUE)
  writeRaster(PA_60_STK, paste0("./data/outputs/predictors/XP4/rcp60/stk_PA_60_", i,".tif"), format = "GTiff", overwrite = TRUE) 
  writeRaster(PA_85_STK, paste0("./data/outputs/predictors/XP4/rcp85/stk_PA_85_", i,".tif"), format = "GTiff", overwrite = TRUE)
}

##      .......... stk_suit                          ----

# SUITrcp26.mat <- do.call(rbind, lapply(SUITrcp26, matrix, ncol = 21, byrow = TRUE))
pres_STK  <- rowMeans(SUITpres)  / ncol(SUITpres)
SUIT_c_STK  <- rasterFromXYZ(cbind(coords, pres_STK))
names(SUIT_c_STK) <- SPid_stk
writeRaster(SUIT_c_STK,  "./data/outputs/predictors/XP5/stk_SUIT_c.tif",  format = "GTiff")

for (i in AOGCMs)
{
  rcp26_STK <- rowMeans(SUITrcp26[[i]]) / ncol(SUITrcp26[[i]]) 
  rcp45_STK <- rowMeans(SUITrcp26[[i]]) / ncol(SUITrcp45[[i]])
  rcp60_STK <- rowMeans(SUITrcp26[[i]]) / ncol(SUITrcp60[[i]]) 
  rcp85_STK <- rowMeans(SUITrcp26[[i]]) / ncol(SUITrcp85[[i]])
  
  SUIT_26_STK <- rasterFromXYZ(cbind(coords, rcp26_STK))
  SUIT_45_STK <- rasterFromXYZ(cbind(coords, rcp45_STK))
  SUIT_60_STK <- rasterFromXYZ(cbind(coords, rcp60_STK))
  SUIT_85_STK <- rasterFromXYZ(cbind(coords, rcp85_STK))
  
  names(SUIT_26_STK) <- names(SUIT_45_STK) <- names(SUIT_60_STK) <- names(SUIT_85_STK) <- SPid_stk
  
  writeRaster(SUIT_26_STK, paste0("./data/outputs/predictors/XP5/rcp26/stk_SUIT_26_", i,".tif"), format = "GTiff", overwrite = TRUE)
  writeRaster(SUIT_45_STK, paste0("./data/outputs/predictors/XP5/rcp45/stk_SUIT_45_", i,".tif"), format = "GTiff", overwrite = TRUE)
  writeRaster(SUIT_60_STK, paste0("./data/outputs/predictors/XP5/rcp60/stk_SUIT_60_", i,".tif"), format = "GTiff", overwrite = TRUE)
  writeRaster(SUIT_85_STK, paste0("./data/outputs/predictors/XP5/rcp85/stk_SUIT_85_", i,".tif"), format = "GTiff", overwrite = TRUE)
}

# ***************************************************************************************
## 08.  SEP_PA                                       ----
###.....................................
# Biotic Predictor	 SEP/PA - Plants XP1.2
# Variables       	 abiotic + SEP/PA -  (12 vars = 5  + 7 )
# Input           	 bee

###.............. Running the modedelling experiment
sp_name
for (i in 1:length(sp_name))
{
  result <- multiple_ENMs(occurrence       = paste0("./data/occurrences/var_",  sp_name[i], ".txt"),
                          background       = paste0("./data/occurrences/back_", sp_name[i], ".txt"),
                          biovar_current   = "./data/climatic_vars/selected/current/",
                          biovar_rcp26     = "./data/climatic_vars/selected/rcp26/",
                          biovar_rcp45     = "./data/climatic_vars/selected/rcp45/",
                          biovar_rcp60     = "./data/climatic_vars/selected/rcp60/",
                          biovar_rcp85     = "./data/climatic_vars/selected/rcp85/",
                          newvar_current   = "./data/outputs/predictors/XP2/",
                          newvar_rcp26     = "./data/outputs/predictors/XP2/rcp26/",
                          newvar_rcp45     = "./data/outputs/predictors/XP2/rcp45/",
                          newvar_rcp60     = "./data/outputs/predictors/XP2/rcp60/",
                          newvar_rcp85     = "./data/outputs/predictors/XP2/rcp85/",
                          trainning        = "./data/occurrences/subsets/trainning",
                          testing          = "./data/occurrences/subsets/testing",
                          AOGCMs           = c(1, 2, 3),
                          Pout             = paste0("./data/outputs/XP2/Pout/xp2_", sp_name[i], "_"),
                          cross_validation = 10)
  
  ###.............. Saving evaluation data
  write.table(result[["TPR_c"]],       paste0("./data/outputs/XP2/", sp_name[i], "_TPR_current.txt"), sep = "\t", row.names = F)
  write.table(result[["Threshold_c"]], paste0("./data/outputs/XP2/", sp_name[i], "_t_current.txt"),   sep = "\t", row.names = F)
  write.table(result[["Pred_area_c"]], paste0("./data/outputs/XP2/", sp_name[i], "_d_current.txt"),   sep = "\t", row.names = F)
  write.table(result[["Threshold_sens"]],  paste0("./data/outputs/XP2/", sp_name[i], "_t_sens_current.txt"),   sep = "\t", row.names = F)
  write.table(result[["AUC"]],  paste0("./data/outputs/XP2/", sp_name[i], "_AUC_current.txt"),   sep = "\t", row.names = F)
  
  rm(result)
  gc()
  
}
beep(8)


# ***************************************************************************************

## 09.  SEP_SUIT                                     ----
###.....................................
# Biotic Predictor	 SEP/SUIT - Plants XP1.2
# Variables       	 abiotic + SEP/SUIT -  (12 vars = 5  + 7 )
# Input           	 bee

###.............. Running the modedelling experiment
sp_name
for (i in 1:length(sp_name))
{
  result <- multiple_ENMs(occurrence       = paste0("./data/occurrences/var_",  sp_name[i], ".txt"),
                          background       = paste0("./data/occurrences/back_", sp_name[i], ".txt"),
                          biovar_current   = "./data/climatic_vars/selected/current/",
                          biovar_rcp26     = "./data/climatic_vars/selected/rcp26/",
                          biovar_rcp45     = "./data/climatic_vars/selected/rcp45/",
                          biovar_rcp60     = "./data/climatic_vars/selected/rcp60/",
                          biovar_rcp85     = "./data/climatic_vars/selected/rcp85/",
                          newvar_current   = "./data/outputs/predictors/XP3/",
                          newvar_rcp26     = "./data/outputs/predictors/XP3/rcp26/",
                          newvar_rcp45     = "./data/outputs/predictors/XP3/rcp45/",
                          newvar_rcp60     = "./data/outputs/predictors/XP3/rcp60/",
                          newvar_rcp85     = "./data/outputs/predictors/XP3/rcp85/",
                          trainning        = "./data/occurrences/subsets/trainning",
                          testing          = "./data/occurrences/subsets/testing",
                          AOGCMs           = c(1, 2, 3),
                          Pout             = paste0("./data/outputs/XP3/Pout/xp3_", sp_name[i], "_"),
                          cross_validation = 10)
  
  ###.............. Saving evaluation data
  write.table(result[["TPR_c"]],       paste0("./data/outputs/XP3/", sp_name[i], "_TPR_current.txt"), sep = "\t", row.names = F)
  write.table(result[["Threshold_c"]], paste0("./data/outputs/XP3/", sp_name[i], "_t_current.txt"),   sep = "\t", row.names = F)
  write.table(result[["Pred_area_c"]], paste0("./data/outputs/XP3/", sp_name[i], "_d_current.txt"),   sep = "\t", row.names = F)
  write.table(result[["Threshold_sens"]],  paste0("./data/outputs/XP3/", sp_name[i], "_t_sens_current.txt"),   sep = "\t", row.names = F)
  write.table(result[["AUC"]],  paste0("./data/outputs/XP3/", sp_name[i], "_AUC_current.txt"),   sep = "\t", row.names = F)
  
  rm(result)
  gc()
}
beep(8)

# ***************************************************************************************

## 10.  STK_PA                                       ----
###.....................................
# Biotic Predictor	 STK/PA - Plants XP1.2
# Variables       	 abiotic + STK/PA -  (6 vars = 5  + 1 )
# Input           	 bee

###.............. Running the modedelling experiment
sp_name
for (i in 1:length(sp_name))
{
  result <- multiple_ENMs(occurrence       = paste0("./data/occurrences/var_",  sp_name[i], ".txt"),
                          background       = paste0("./data/occurrences/back_", sp_name[i], ".txt"),
                          biovar_current   = "./data/climatic_vars/selected/current/",
                          biovar_rcp26     = "./data/climatic_vars/selected/rcp26/",
                          biovar_rcp45     = "./data/climatic_vars/selected/rcp45/",
                          biovar_rcp60     = "./data/climatic_vars/selected/rcp60/",
                          biovar_rcp85     = "./data/climatic_vars/selected/rcp85/",
                          newvar_current   = "./data/outputs/predictors/XP4/",
                          newvar_rcp26     = "./data/outputs/predictors/XP4/rcp26/",
                          newvar_rcp45     = "./data/outputs/predictors/XP4/rcp26/",
                          newvar_rcp60     = "./data/outputs/predictors/XP4/rcp26/",
                          newvar_rcp85     = "./data/outputs/predictors/XP4/rcp26/",
                          trainning        = "./data/occurrences/subsets/trainning",
                          testing          = "./data/occurrences/subsets/testing",
                          AOGCMs           = c(1, 2, 3),
                          Pout             = paste0("./data/outputs/XP4/Pout/xp4_", sp_name[i], "_"),
                          cross_validation = 10)
  
  
  ###.............. Saving evaluation data
  write.table(result[["TPR_c"]],       paste0("./data/outputs/XP4/", sp_name[i], "_TPR_current.txt"), sep = "\t", row.names = F)
  write.table(result[["Threshold_c"]], paste0("./data/outputs/XP4/", sp_name[i], "_t_current.txt"),   sep = "\t", row.names = F)
  write.table(result[["Pred_area_c"]], paste0("./data/outputs/XP4/", sp_name[i], "_d_current.txt"),   sep = "\t", row.names = F)
  write.table(result[["Threshold_sens"]],  paste0("./data/outputs/XP4/", sp_name[i], "_t_sens_current.txt"),   sep = "\t", row.names = F)
  write.table(result[["AUC"]],  paste0("./data/outputs/XP4/", sp_name[i], "_AUC_current.txt"),   sep = "\t", row.names = F)
  
  rm(result)
  gc()
}
beep(8)


# ***************************************************************************************

## 11.  STK_SUIT                                     ----
###.....................................
# Biotic Predictor	 STK/SUIT - Plants XP1.2
# Variables       	 abiotic + STK/SUIT -  (6 vars = 5  + 1 )
# Input           	 bee

###.............. Running the modedelling experiment
sp_name
for (i in 1:length(sp_name))
{
  result <- multiple_ENMs(occurrence       = paste0("./data/occurrences/var_",  sp_name[i], ".txt"),
                          background       = paste0("./data/occurrences/back_", sp_name[i], ".txt"),
                          biovar_current   = "./data/climatic_vars/selected/current/",
                          biovar_rcp26     = "./data/climatic_vars/selected/rcp26/",
                          biovar_rcp45     = "./data/climatic_vars/selected/rcp45/",
                          biovar_rcp60     = "./data/climatic_vars/selected/rcp60/",
                          biovar_rcp85     = "./data/climatic_vars/selected/rcp85/",
                          newvar_current   = "./data/outputs/predictors/XP5/",
                          newvar_rcp26     = "./data/outputs/predictors/XP5/rcp26/",
                          newvar_rcp45     = "./data/outputs/predictors/XP5/rcp45/",
                          newvar_rcp60     = "./data/outputs/predictors/XP5/rcp60/",
                          newvar_rcp85     = "./data/outputs/predictors/XP5/rcp85/",
                          trainning        = "./data/occurrences/subsets/trainning",
                          testing          = "./data/occurrences/subsets/testing",
                          AOGCMs           = c(1, 2, 3),
                          Pout             = paste0("./data/outputs/XP5/Pout/xp5_", sp_name[i], "_"),
                          cross_validation = 10)
  
  ###.............. Saving evaluation data
  write.table(result[["TPR_c"]],       paste0("./data/outputs/XP5/", sp_name[i], "_TPR_current.txt"), sep = "\t", row.names = F)
  write.table(result[["Threshold_c"]], paste0("./data/outputs/XP5/", sp_name[i], "_t_current.txt"),   sep = "\t", row.names = F)
  write.table(result[["Pred_area_c"]], paste0("./data/outputs/XP5/", sp_name[i], "_d_current.txt"),   sep = "\t", row.names = F)
  write.table(result[["Threshold_sens"]],  paste0("./data/outputs/XP5/", sp_name[i], "_t_sens_current.txt"),   sep = "\t", row.names = F)
  write.table(result[["AUC"]],  paste0("./data/outputs/XP5/", sp_name[i], "_AUC_current.txt"),   sep = "\t", row.names = F)
  
  rm(result)
  gc()
}
beep(8)


# ***************************************************************************************
## 12.  Ensembles and Final Outputs                  ----
#   .....................................................................................
source("./ensemble.R")

#   ...............XP1
ALLsp_names
for (i in 1:length(ALLsp_names))
{
  result <- ensemble(Pout           =  "./data/outputs/XP1/Pout/",
                     Alld           =  paste0("./data/outputs/XP1/", ALLsp_names[i], "_d_current.txt"),
                     sp             =  ALLsp_names[i],
                     AOGCMs         =  c(1, 2, 3),
                     biovar_current =  "./data/climatic_vars/selected/current/")
  
### Saving predictions
  
### Saving Ensembles
  write.table(result[["Ensemble"]],   paste0("./data/outputs/ensembles/xp1_", ALLsp_names[i], "_ENSEMBLES.txt"), sep = "\t", row.names = F)
  
  rm(result)
  gc()
}
beep(8)

### Plotting Ensembles
PlotEnsemble(occur     = "./data/occurrences/occur_thinned.txt",
             ensb_data = "./data/outputs/ensembles/xp1_",
             sp        = ALLsp_names,
             output    = "./data/outputs/ensembles/xp1_plot_ind_")

#   ...............XP2

sp_name
for (i in 1:length(sp_name))
{
  result <- ensemble(Pout           = "./data/outputs/XP2/Pout/",
                     Alld           =  paste0("./data/outputs/XP2/", sp_name[i], "_d_current.txt"),
                     sp             = sp_name,
                     AOGCMs         = c(1, 2, 3),
                     biovar_current = "./data/climatic_vars/selected/current/")
  
### Saving Ensembles
  write.table(result[["Ensemble"]], paste0("./data/outputs/ensembles/xp2_", sp_name[i], "_ENSEMBLES.txt"), sep = "\t", row.names = F)
  
  rm(result)
  gc()
}
beep(8)

### Plotting Ensembles
PlotEnsemble(occur     = "./data/occurrences/occur_thinned.txt",
             ensb_data = "./data/outputs/ensembles/xp2_",
             sp        = sp_name,
             output    = "./data/outputs/ensembles/xp2_plot_ind")

#   ...............XP3

sp_name
for (i in 1:length(sp_name))
{
  result <- ensemble(Pout           = "./data/outputs/XP3/Pout/",
                     Alld           =  paste0("./data/outputs/XP3/", sp_name[i], "_d_current.txt"),
                     sp             = sp_name[i],
                     AOGCMs         = c(1, 2, 3),
                     biovar_current = "./data/climatic_vars/selected/current/")

### Saving Ensembles
  write.table(result[["Ensemble"]], paste0("./data/outputs/ensembles/xp3_", sp_name[i], "_ENSEMBLES.txt"), sep = "\t", row.names = F)
  
  rm(result)
  gc()
}
beep(8)

### Plotting Ensembles
PlotEnsemble(occur     = "./data/occurrences/occur_thinned.txt",
             ensb_data = "./data/outputs/ensembles/xp3_",
             sp        = sp_name,
             output    = "./data/outputs/ensembles/xp3_plot_ind")


#   ...............XP4

sp_name
for (i in 1:length(sp_name))
{
  result <- ensemble(Pout           = "./data/outputs/XP4/Pout/",
                     Alld           =  paste0("./data/outputs/XP4/", sp_name[i], "_d_current.txt"),
                     sp             = sp_name[i],
                     AOGCMs         = c(1, 2, 3),
                     biovar_current = "./data/climatic_vars/selected/current/")
  
### Saving Ensembles
  write.table(result[["Ensemble"]], paste0("./data/outputs/ensembles/xp4_", sp_name[i], "_ENSEMBLES.txt"), sep = "\t", row.names = F)
  
  rm(result)
  gc()
}
beep(8)

### Plotting Ensembles
PlotEnsemble(occur     = "./data/occurrences/occur_thinned.txt",
             ensb_data = "./data/outputs/ensembles/xp4_",
             sp        = sp_name,
             output    = "./data/outputs/ensembles/xp4_plot_ind")

#   ...............XP5

sp_name
for (i in 1:length(sp_name))
{
  result <- ensemble(Pout           = "./data/outputs/XP5/Pout/",
                     Alld           =  paste0("./data/outputs/XP5/", sp_name[i], "_d_current.txt"),
                     sp             = sp_name[i],
                     AOGCMs         = c(1, 2, 3),
                     biovar_current = "./data/climatic_vars/selected/current/")
  
### Saving predictions
  # writeRaster(result[["output_current"]], paste0("./data/outputs/XP5/", sp_name[i],"_current.tif"), format = "GTiff")
  # writeRaster(result[["output_rcp26"]],   paste0("./data/outputs/XP5/", sp_name[i],"_rcp26.tif"), format = "GTiff")
  # writeRaster(result[["output_rcp45"]],   paste0("./data/outputs/XP5/", sp_name[i],"_rcp45.tif"), format = "GTiff")
  # writeRaster(result[["output_rcp60"]],   paste0("./data/outputs/XP5/", sp_name[i],"_rcp60.tif"), format = "GTiff")
  # writeRaster(result[["output_rcp85"]],   paste0("./data/outputs/XP5/", sp_name[i],"_rcp85.tif"), format = "GTiff")
  
### Saving Ensembles
  write.table(result[["Ensemble"]], paste0("./data/outputs/ensembles/xp5_", sp_name[i], "_ENSEMBLES.txt"), sep = "\t", row.names = F)
  
  rm(result)
  gc()
}
beep(8)

### Plotting Ensembles
PlotEnsemble(occur     = "./data/occurrences/occur_thinned.txt",
             ensb_data = "./data/outputs/ensembles/xp5_",
             sp        = sp_name,
             output    = "./data/outputs/ensembles/xp5_plot_ind")

# ***************************************************************************************
## 12.  Selecting treatment                          ----
##      .......... by d                              ----
# Anova de Medidas Repetidas (mesmos subconjuntos de XP2 a XP5)	
# preditor:  XP1:XP5	
# resposta: estatística "d" para cada método	

# XP1 <- read.table("./data/outputs/XP1/Lithurgus_huberi_d_current.txt", sep = "\t", h = T)
# XP2 <- read.table("./data/outputs/XP2/Lithurgus_huberi_d_current.txt", sep = "\t", h = T)
# XP3 <- read.table("./data/outputs/XP3/Lithurgus_huberi_d_current.txt", sep = "\t", h = T)
# XP4 <- read.table("./data/outputs/XP4/Lithurgus_huberi_d_current.txt", sep = "\t", h = T)
# XP5 <- read.table("./data/outputs/XP5/Lithurgus_huberi_d_current.txt", sep = "\t", h = T)
# 
# CV_only  <- rep("CV-only", 10)
# sep_pa   <- rep("sep_pa", 10)
# sep_suit <- rep("sep_suit", 10)
# stk_pa   <- rep("stk_pa", 10)
# stk_suit <- rep("stk_suit", 10)
# 
# CV_only  <- cbind(CV_only, XP1)
# sep_pa   <- cbind(sep_pa, XP2)
# sep_suit <- cbind(sep_suit, XP3)
# stk_pa   <- cbind(stk_pa, XP4)
# stk_suit <- cbind(stk_suit, XP5)
# 
# names(CV_only) <- names(sep_pa) <- names(sep_suit) <- names(stk_pa) <- names(stk_suit) <- c("id","bioclim", "maxent", "SVM")
# 
# 
# data <- as.data.frame(rbind(CV_only, sep_pa, sep_suit, stk_pa, stk_suit))
# head(data)
# data_stack <- data.frame(data[1], stack(data[2:4]))
# head(data_stack)
# tail(data_stack)
# 
# write.table(data_stack, "Lithurgus_anorep_suit.txt", sep = "\t")
# 
# if(!require(lme4)){install.packages("lme4")}
# if(!require(lmerTest)){install.packages("lmerTest")}
# library(car)
# factordata <- as.fac
# 
# anova <- aov(values~id, data = data_stack)
# 
#  rmaov <- lmer(values ~ id + (1|ind), data = data_stack)
# anova(rmaov)
# 
# data_stack$values <- factor(data_stack$values)
# data_stack$id <- factor(data_stack$id)
# 
# boxplot(id ~ values, data_stack)

# .................................................


##      .......... by range                          ----

# predictor: XP1:XP5	
# resposta: tamanho de range

#........... calculating range.
source("./get_range.R")

#...........
XP1 <- get.range(biovar_current = "./data/climatic_vars/selected/current/",
                 sp_name        = "Lithurgus_huberi_",
                 pout           = "./data/outputs/XP1/Pout/xp1_",
                 thr            = "./data/outputs/XP1/",
                 cross_val      = 10,
                 AOGCMs         = 3)

write.table(XP1[["pres"]], "./data/outputs/range/xp1_Lithurgus_huberi_current.txt", sep = "\t", row.names = F)
write.table(XP1[["fut26"]], "./data/outputs/range/xp1_Lithurgus_huberi_rcp26.txt", sep = "\t", row.names = F)
write.table(XP1[["fut45"]], "./data/outputs/range/xp1_Lithurgus_huberi_rcp45.txt", sep = "\t", row.names = F)
write.table(XP1[["fut60"]], "./data/outputs/range/xp1_Lithurgus_huberi_rcp60.txt", sep = "\t", row.names = F)
write.table(XP1[["fut85"]], "./data/outputs/range/xp1_Lithurgus_huberi_rcp85.txt", sep = "\t", row.names = F)

#...........
XP2 <- get.range(biovar_current = "./data/climatic_vars/selected/current",
                 sp_name        = "Lithurgus_huberi_",
                 pout           = "./data/outputs/XP2/Pout/xp2_",
                 thr            = "./data/outputs/XP2/",
                 cross_val      = 10,
                 AOGCMs         = 3)

write.table(XP2[["pres"]], "./data/outputs/range/xp2_Lithurgus_huberi_current_sumgrid.txt", sep = "\t", row.names = F)
write.table(XP2[["fut26"]], "./data/outputs/range/xp2_Lithurgus_huberi_rcp26_sumgrid.txt", sep = "\t", row.names = F)
write.table(XP2[["fut45"]], "./data/outputs/range/xp2_Lithurgus_huberi_rcp45_sumgrid.txt", sep = "\t", row.names = F)
write.table(XP2[["fut60"]], "./data/outputs/range/xp2_Lithurgus_huberi_rcp60_sumgrid.txt", sep = "\t", row.names = F)
write.table(XP2[["fut85"]], "./data/outputs/range/xp2_Lithurgus_huberi_rcp85_sumgrid.txt", sep = "\t", row.names = F)

#...........
XP3 <- get.range(biovar_current = "./data/climatic_vars/selected/current/",
                 sp_name        = "Lithurgus_huberi_",
                 pout           = "./data/outputs/XP3/Pout/xp3_",
                 thr            = "./data/outputs/XP3/",
                 cross_val      = 10,
                 AOGCMs         = 3)

write.table(XP3[["pres"]], "./data/outputs/range/xp3_Lithurgus_huberi_current.txt", sep = "\t", row.names = F)
write.table(XP3[["fut26"]], "./data/outputs/range/xp3_Lithurgus_huberi_rcp26.txt", sep = "\t", row.names = F)
write.table(XP3[["fut45"]], "./data/outputs/range/xp3_Lithurgus_huberi_rcp45.txt", sep = "\t", row.names = F)
write.table(XP3[["fut60"]], "./data/outputs/range/xp3_Lithurgus_huberi_rcp60.txt", sep = "\t", row.names = F)
write.table(XP3[["fut85"]], "./data/outputs/range/xp3_Lithurgus_huberi_rcp85.txt", sep = "\t", row.names = F)

#...........
XP4 <- get.range(biovar_current = "./data/climatic_vars/selected/current/",
                 sp_name        = "Lithurgus_huberi_",
                 pout           = "./data/outputs/XP4/Pout/xp4_",
                 thr            = "./data/outputs/XP4/",
                 cross_val      = 10,
                 AOGCMs         = 3)

write.table(XP4[["pres"]], "./data/outputs/range/xp4_Lithurgus_huberi_current.txt", sep = "\t", row.names = F)
write.table(XP4[["fut26"]], "./data/outputs/range/xp4_Lithurgus_huberi_rcp26.txt", sep = "\t", row.names = F)
write.table(XP4[["fut45"]], "./data/outputs/range/xp4_Lithurgus_huberi_rcp45.txt", sep = "\t", row.names = F)
write.table(XP4[["fut60"]], "./data/outputs/range/xp4_Lithurgus_huberi_rcp60.txt", sep = "\t", row.names = F)
write.table(XP4[["fut85"]], "./data/outputs/range/xp4_Lithurgus_huberi_rcp85.txt", sep = "\t", row.names = F)

#...........
XP5 <- get.range(biovar_current = "./data/climatic_vars/selected/current/",
                 sp_name        = "Lithurgus_huberi_",
                 pout           = "./data/outputs/XP5/Pout/xp5_",
                 thr            = "./data/outputs/XP5/",
                 cross_val      = 10,
                 AOGCMs         = 3)

write.table(XP5[["pres"]], "./data/outputs/range/xp5_Lithurgus_huberi_current.txt", sep = "\t", row.names = F)
write.table(XP5[["fut26"]], "./data/outputs/range/xp5_Lithurgus_huberi_rcp26.txt", sep = "\t", row.names = F)
write.table(XP5[["fut45"]], "./data/outputs/range/xp5_Lithurgus_huberi_rcp45.txt", sep = "\t", row.names = F)
write.table(XP5[["fut60"]], "./data/outputs/range/xp5_Lithurgus_huberi_rcp60.txt", sep = "\t", row.names = F)
write.table(XP5[["fut85"]], "./data/outputs/range/xp5_Lithurgus_huberi_rcp85.txt", sep = "\t", row.names = F)


## 13.  Preparing analysis factors                   ----
###.....................................


back <- list.files("./data/occurrences/", pattern = "back", full.names = TRUE)
back <- bind_rows(lapply(back, fread))
back <- as.data.frame(back)

### Reading predictions data
all_outputs  <- laaply (list.files("./data/outputs/", patern = ".bil$", full.names = TRUE), stack)

all_d        <- laaply (list.files("./data/outputs/", patern = "_d_",   full.names = TRUE), fread)

all_TRP      <- laaply (list.files("./data/outputs/", patern = "_TRP_", full.names = TRUE), fread)

all_t        <- laaply (list.files("./data/outputs/", patern = "_t_",   full.names = TRUE), fread)

### Standardizing suitabilities
all_val <- values(all_output)
all_pad <- deconstante(all_val, "standardize", 2)

# ***************************************************************************************
## 14.  Uncertainty Evaluation                       ----

all_huberi <- stack(huberi_c, huberi_rcp26, huberi_rcp45, huberi_rcp60, huberi_rcp85)
TPR_huberi <- TPR_h[which(TPR_h)]

##  Creating ANOVA Factors to be used at Uncertainty Evaluation
# huberi_metodo <- rep("bioclim", nlayers(all_huberi)) #???
# huberi_tempo  <- rep(c("pres","fut"), each = nlayers(all_huberi))
## huberi

data   <- values(stack(bioclim_h, gower_h, maha_h, maxent_h, SVM_h, GLM_h))
method <- c(bioclim_method_h, gower_method_h, maha_method_h, maxent_method_h, SVM_method_h, GLM_method_h)
period <- c(bioclim_period_h, gower_period_h, maha_period_h, maxent_period_h, SVM_period_h, GLM_period_h)

# loop the ANOVA



# ****  List of improvements to the scritp           ----

# 1. Implement validation by the checkerboards method.
# 2. Implement occurrence filtering at the ambiental space and compare with the geographical space one (spThin).
# 3. Transform maps in frequencies instead of suitabilities.
# 4. Impelement multi cores for running several models simultaneously.
# 5. Reduce code by implementing subfunctions, lopps.
# 6. Rewrite the code using tidy



