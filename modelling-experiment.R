# *** Loading functions                             ----

source("./Auxiliary_functions.R")
source("./Multiple_ENMs.R")

# *** Loading packages                              ----

# install.packages(c("tidyverse", "raster", "rgdal", "abind", "spThin", "vegan", "maps", "pcych", "kernlab", "dismo", "rJava", "dendextend"))

load_pak(c("tidyverse", "raster", "rgdal", "abind", "spThin", "dismo", "kernlab", "vegan", "maps", "psych", "rJava", "dendextend"))

# ***************************************************************************************
## 01. Read aogcms models                           ----
##### 1.1.............. current                                     

current <- read_current(dir = "./data/climatic_vars/current/")
current_spatial <- rasterFromXYZ(current) 

##### 1.2.............. rcp

rcp26_list <- read_rcp( x = "./data/climatic_vars/selected/26bi70/")
rcp45_list <- read_rcp( x = "./data/climatic_vars/selected/45bi70/")
rcp60_list <- read_rcp( x = "./data/climatic_vars/selected/60bi70/")
rcp85_list <- read_rcp( x = "./data/climatic_vars/selected/85bi70/")

# ### Creating array
# rcp26 <- rcp26_list[["array"]]
# rcp45 <- rcp45_list[["array"]]
# rcp60 <- rcp60_list[["array"]]
# rcp85 <- rcp85_list[["array"]]

### Creating RasterStack
#. Following the order of the models in the directory:
# 1 == "CCS4"
# 2 == "IPSL-CMSA-LR"
# 3 == "MIROC-ESM"

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
### By exploratory factor analysis

fa.parallel(current[ , -c(1:2)], fa = 'fa') #scree plot
current_fa <- fa(current[ , -c(1:2)], nfactors = 5, rotate = 'varimax')
loadings <- loadings(current_fa)
write.table(loadings, "./data/climatic_vars/selected/varimax_loadings.txt")

### Selected variables 
# bio02, bio03, bio10, bio14, bio16.

# ***************************************************************************************
## 03. Saving selected variables                    ----

###................ CURRENT

# ### 3.1a current - saving as table 
# 
# write.table(current[,c("x", "y", "bio02", "bio03", "bio10", "bio14", "bio16")], "./data/climatic_vars/selected/current/current-select.txt", row.names = F, sep = " ")  

### 3.1b current - saving as raster

variables <- as.factor(c("bio02", "bio03", "bio10", "bio14", "bio16"))
for (i in 1:length(variables))
{
  writeRaster (current_spatial[[i]], filename = paste0("./data/climatic_vars/selected/current/current-", variables[i], ".grd"), format = "raster")
}
rm(variables)

current_select <- stack(list.files("./data/climatic_vars/selected/current/",  pattern = ".grd$", full.names = TRUE))

###................ RCP

# ### 3.2a rcp - saving table
# 
# write.table(rcp26 [ ,c("x", "y", "bio02", "bio03", "bio10", "bio14", "bio16" ), ], "./data/climatic_vars/selected/rcp26/rcp26-select.txt", row.names = F, sep = "	")
# write.table(rcp45 [ ,c("x", "y", "bio02", "bio03", "bio10", "bio14", "bio16" ), ], "./data/climatic_vars/selected/rcp45/rcp45-select.txt", row.names = F, sep = "	")
# write.table(rcp60 [ ,c("x", "y", "bio02", "bio03", "bio10", "bio14", "bio16" ), ], "./data/climatic_vars/selected/rcp60/rcp60-select.txt", row.names = F, sep = "	")
# write.table(rcp85 [ ,c("x", "y", "bio02", "bio03", "bio10", "bio14", "bio16" ), ], "./data/climatic_vars/selected/rcp85/rcp85-select.txt", row.names = F, sep = "	")


### 3.2b rcp - saving as raster

variables <- as.factor(c("bio02.1", "bio03.1", "bio10.1", "bio14.1", "bio16.1", 
                         "bio02.2", "bio03.2", "bio10.2", "bio14.2", "bio16.2", 
                         "bio02.3", "bio03.3", "bio10.3", "bio14.3", "bio16.3"))

## RCP26
for (i in variables)
{
  writeRaster(rcp26_spatial[[which(names(rcp26_spatial) %in% i)]], filename = paste0("./data/climatic_vars/selected/rcp26/rcp26-", i, ".grd"), format = "raster")
}

## RCP45
for (i in variables)
{
  writeRaster (rcp45_spatial[[which(names(rcp45_spatial) %in% i)]], filename = paste0("./data/climatic_vars/selected/rcp45/rcp45-", i, ".grd"), format = "raster")
}

## RCP60
for (i in variables)
{
  writeRaster (rcp60_spatial[[which(names(rcp60_spatial) %in% i)]], filename = paste0("./data/climatic_vars/selected/rcp60/rcp60-", i, ".grd"), format = "raster")
}

## RCP85
for (i in variables)
{
  writeRaster (rcp85_spatial[[which(names(rcp85_spatial) %in% i)]], filename = paste0("./data/climatic_vars/selected/rcp85/rcp85-", i, ".grd"), format = "raster")
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
  occur_thinned <- rbind(occur_thinned, occur)
}
write.table(occur_thinned, "./data/occurrences/occur_thinned.txt", sep = ";", row.names = FALSE)

# occur_thinned <- read.csv("./data/occurrences/occur_thinned.csv", sep = ",", h = T)

# ***************************************************************************************
## 05. Extracting variables                         ----

# Creating and saving the object "var" for each studied species
for(i in 1:length(sp_names))
{
  var <- create_var(occur_thinned[occur_thinned[, 1] == sp_names[i], ], sp_names[i])
}

# ***************************************************************************************
## 06. Background Sampling                          ----

###.............. Creating and saving the object "back" for each studied species

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
## 07. Running our model                            ----
rm(list = ls())

occur_thinned <- read.csv("./data/occurrences/raw_data.txt", sep = "")
sp <- gsub("C[1-9]","", occur_thinned$SPEC)
sp_names <- unique(sp)
sp_names
for (i in 1:length(sp_names))
{
  ###.............. Running the modedelling experiment
  result <- multiple_ENMs(occurrence       = paste0("./data/occurrences/var_", sp_names[i], ".txt"),
                          background       = paste0("./data/occurrences/back_", sp_names[i],".txt"),
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
  #..... current
  write.table(result[["TPR_c"]], paste0("./data/outputs/", sp_names[i], "_TPR_current.txt"), sep = "\t", row.names = F)
  write.table(result[["Threshold_c"]], paste0("./data/outputs/", sp_names[i], "_t_current.txt"), sep = "\t", row.names = F)
  write.table(result[["Predicted_area_c"]], paste0("./data/outputs/", sp_names[i], "_d_current.txt"), sep = "\t", row.names = F)
  
  #..... rcp26
  write.table(result[["TPR_rcp26"]], paste0("./data/outputs/", sp_names[i], "_TPR_rcp26.txt"), sep = "\t", row.names = F)
  write.table(result[["Threshold_rcp26"]], paste0("./data/outputs/", sp_names[i], "_t_rcp26.txt"), sep = "\t", row.names = F)
  write.table(result[["Predicted_area_rcp26"]], paste0("./data/outputs/", sp_names[i], "_d_rcp26.txt"), sep = "\t", row.names = F)
  
  #..... rcp45
  write.table(result[["TPR_rcp45"]], paste0("./data/outputs/", sp_names[i], "_TPR_rcp45.txt"), sep = "\t", row.names = F)
  write.table(result[["Threshold_rcp45"]], paste0("./data/outputs/", sp_names[i], "_t_rcp45.txt"), sep = "\t", row.names = F)
  write.table(result[["Predicted_area_rcp45"]], paste0("./data/outputs/", sp_names[i], "_d_rcp45.txt"), sep = "\t", row.names = F)
  
  #..... rcp60
  write.table(result[["TPR_rcp60"]], paste0("./data/outputs/", sp_names[i], "_TPR_rcp60.txt"), sep = "\t", row.names = F)
  write.table(result[["Threshold_rcp60"]], paste0("./data/outputs/", sp_names[i], "_t_rcp60.txt"), sep = "\t", row.names = F)
  write.table(result[["Predicted_area_rcp60"]], paste0("./data/outputs/", sp_names[i], "_d_rcp60.txt"), sep = "\t", row.names = F)
  
  #..... rcp85
  write.table(result[["TPR_rcp85"]], paste0("./data/outputs/", sp_names[i], "_TPR_rcp85.txt"), sep = "\t", row.names = F)
  write.table(result[["Threshold_rcp85"]], paste0("./data/outputs/", sp_names[i], "_t_rcp85.txt"), sep = "\t", row.names = F)
  write.table(result[["Predicted_area_rcp85"]], paste0("./data/outputs/", sp_names[i], "_d_rcp85.txt"), sep = "\t", row.names = F)
}
  
# ***************************************************************************************
## 08. Standardize suitabilities (suit)             ----


# ***************************************************************************************
## 09. Ensemble                                     ----

## huberi


# loop Ensemble for huberi

# Expected of ensembles for huberi:
# huberi current;
# huberi rcp26;
# huberi rcp45;
# huberi rcp60;
# huberi rcp85.


# Expected of ensembles for each host plants:
# host plants current;
# host plants rcp26;
# host plants rcp45;
# host plants rcp60;
# host plants rcp85.


# ***************************************************************************************
## 10. Uncertainty Evaluation                       ----

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

## host plants
data   <- values(stack(bioclim_p, gower_p, maha_p, maxent_p, SVM_p, GLM_p))
method <- c(bioclim_method_p, gower_method_p, maha_method_p, maxent_method_p, SVM_method_p, GLM_method_p)
period <- c(bioclim_period_p, gower_period_p, maha_period_p, maxent_period_p, SVM_period_p, GLM_period_p)

# loop the ANOVA






#  *** List of improvements to the scritp           ----

# 1. Implement validation by the checkerboards method.
# 2. Implement occurrence filtering at the ambiental space and compare with spThin (geographical space)
# 3. Transform maps in frequencies instead of suitabilities.
# 4. Impelement multi cores for running several models simultaneously.
# 5. Reduce code by implementing subfunctions, lopps.
# 6. Rewrite the code using tidy

