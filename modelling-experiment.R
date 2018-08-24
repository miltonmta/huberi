file.edit("readme.R")
# **** Loading functions                            ----

source("./Auxiliary_functions.R")
source("./Multiple_ENMs.R")
# file.edit(c("./Auxiliary_functions.R", "./Multiple_ENMs.R"))

# **** Loading packages                             ----

# install.packages(c("tidyverse", "raster", "rgdal", "abind", "spThin", "vegan", "maps", "pcych", "kernlab", "dismo", "rJava", "dendextend", "beepr"))

load_pak(c("tidyverse", "raster", "rgdal", "abind", "spThin", "dismo", "kernlab", "vegan", "maps", "psych", "rJava", "dendextend", "beepr", "data.table"))

# ***************************************************************************************
## 01. Read aogcms models                           ----
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
## 02. Variable selection                           ----
### By exploratory factor analysis

fa.parallel(current_mat[ , -c(1:2)], fa = 'fa') #scree plot
current_fa <- fa(current_mat[ , -c(1:2)], nfactors = 5, rotate = 'varimax')
loadings <- loadings(current_fa)
write.table(loadings, "./data/climatic_vars/selected/varimax_loadings.txt")
beep(2)

### Selected variables 
# bio02, bio03, bio10, bio14, bio16.

# ***************************************************************************************
## 03. Saving selected variables                    ----

###................ CURRENT

variables <- as.factor(c("bio02", "bio03", "bio10", "bio14", "bio16"))
for (i in variables)
{
  writeRaster (current[[which(names(current) %in% i)]], filename = paste0("./data/climatic_vars/selected/current/current-", i, ".grd"), format = "raster")
}
rm(variables)

current_select <- stack(list.files("./data/climatic_vars/selected/current",  pattern = ".grd$", full.names = TRUE))
beep(2)

###................ RCPs

variables <- as.factor(c("bio02.1", "bio03.1", "bio10.1", "bio14.1", "bio16.1", 
                         "bio02.2", "bio03.2", "bio10.2", "bio14.2", "bio16.2", 
                         "bio02.3", "bio03.3", "bio10.3", "bio14.3", "bio16.3"))

## RCP26
for (i in variables)
{
  writeRaster(rcp26[[which(names(rcp26) %in% i)]], filename = paste0("./data/climatic_vars/selected/rcp26/rcp26-", i, ".grd"), format = "raster")
}

## RCP45
for (i in variables)
{
  writeRaster (rcp45[[which(names(rcp45) %in% i)]], filename = paste0("./data/climatic_vars/selected/rcp45/rcp45-", i, ".grd"), format = "raster")
}

## RCP60
for (i in variables)
{
  writeRaster (rcp60[[which(names(rcp60) %in% i)]], filename = paste0("./data/climatic_vars/selected/rcp60/rcp60-", i, ".grd"), format = "raster")
}

## RCP85
for (i in variables)
{
  writeRaster (rcp85[[which(names(rcp85) %in% i)]], filename = paste0("./data/climatic_vars/selected/rcp85/rcp85-", i, ".grd"), format = "raster")
}
beep(2)

rm(variables)

# ***************************************************************************************
## 04. Occurrences data                             ----

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

# occur_thinned <- read.csv("./data/occurrences/occur_thinned.csv", sep = ",", h = T)

###.............. Extrancting bio variables based on the ocurrence cells

## Creating and saving the object "var" for the stacked occurrences but the response one.
resource <- occur_thinned[!grepl("Lithurgus_huberi", occur_thinned$SPEC), ]
write.table(resource, "./data/occurrences/resource.txt", sep = ";", row.names = FALSE)

# resource <- read.table("./data/occurrences/resource.txt", sep = ";", h = T)

var <- create_var(resource, "resource")

## Creating and saving the object "var" for each studied species
for(i in 1:length(sp_names))
{
  var <- create_var(occur_thinned[occur_thinned[, 1] == sp_names[i], ], sp_names[i])
}
beep(2)


# ***************************************************************************************
## 05. Background Sampling                          ----

###.............. Background files 
# Creating and saving the object "back" for the stacked species but the response one

var_file <- read.table("./data/occurrences/var_resource.txt", h = T, sep = ";")
create_back(var_file, "resource")
beep(2)

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

# rm(list = ls())
# ***************************************************************************************
## 06. XP1                                          ----
# Scenario   Current + Future - 4 RCPs
# Variables	 Abiótic ( 5 vars ) # include SOIL only for host plants
# Input	     9 sps :: bee + 7 plants + "resource" * (summed occurs of all plants)
# Output     9 prediction  x 4 predictive models.
# Result  	 9 Ensembles + Biotic Predictor Variables PA_SEP, PA_STK, SUIT_SEP, SUIT_SKT, resourceSEP, resourceSUIT.

### XP1.1 - bee
#   .....................................................................................

occur_thinned <- read.csv("./data/occurrences/occur_thinned.csv", sep = ",")
sp <- occur_thinned[occur_thinned[, 1] == "Lithurgus_huberi", ]
sp <- gsub("C[1-9]","", sp$SPEC)
sp_name <- unique(sp)
sp_name
for (i in 1:length(sp_name))
{
  ###.............. Running the modedelling experiment
  result <- multiple_ENMs(occurrence       = paste0("./data/occurrences/var_",  sp_name[i], ".txt"),
                          background       = paste0("./data/occurrences/back_", sp_name[i], ".txt"),
                          predictor        = , # no aditional predictor for this occur.
                          biovar_current   = "./data/climatic_vars/selected/current/",
                          biovar_rcp26     = "./data/climatic_vars/selected/rcp26/",
                          biovar_rcp45     = "./data/climatic_vars/selected/rcp45/",
                          biovar_rcp60     = "./data/climatic_vars/selected/rcp60/",
                          biovar_rcp85     = "./data/climatic_vars/selected/rcp85/",
                          cross_validation = 20)
  
  ###.............. Saving predictions
  writeRaster(result[["output_current"]], paste0("./data/outputs/XP1", sp_name[i],"_current.bil"), format = "EHdr")
  writeRaster(result[["output_rcp26"]],   paste0("./data/outputs/XP1", sp_name[i],"_rcp26.bil"), format = "EHdr")
  writeRaster(result[["output_rcp45"]],   paste0("./data/outputs/XP1", sp_name[i],"_rcp45.bil"), format = "EHdr")
  writeRaster(result[["output_rcp60"]],   paste0("./data/outputs/XP1", sp_name[i],"_rcp60.bil"), format = "EHdr")
  writeRaster(result[["output_rcp85"]],   paste0("./data/outputs/XP1", sp_name[i],"_rcp85.bil"), format = "EHdr")

  
  ###.............. Saving evaluation data
  write.table(result[["TPR_c"]],       paste0("./data/outputs/XP1", sp_name[i], "_TPR_current.txt"), sep = "\t", row.names = F)
  write.table(result[["Threshold_c"]], paste0("./data/outputs/XP1", sp_name[i], "_t_current.txt"),   sep = "\t", row.names = F)
  write.table(result[["Pred_area_c"]], paste0("./data/outputs/XP1", sp_name[i], "_d_current.txt"),   sep = "\t", row.names = )
  
}
beep(8)

### XP1.2 - 7 plants + resource
#   .....................................................................................

occur_thinned <- read.csv("./data/occurrences/occur_thinned.csv", sep = ",") # all species but "resource"
sp <- occur_thinned[!grepl("Lithurgus_huberi", occur_thinned$SPEC), ]
sp <- gsub("C[1-9]","", sp$SPEC)
sp_names <- unique(sp)
sp_names <- c(sp_names, "resource")
sp_names
for (i in 1:length(sp_names))
{
  ###.............. Running the modedelling experiment
  result <- multiple_ENMs(occurrence       = paste0("./data/occurrences/var_",  sp_names[i], ".txt"),
                          background       = paste0("./data/occurrences/back_", sp_names[i], ".txt"),
                          predictor        = "./data/climatic_vars/selected/soil/", # abiotic vars for soil 
                          biovar_current   = "./data/climatic_vars/selected/current/",
                          biovar_rcp26     = "./data/climatic_vars/selected/rcp26/",
                          biovar_rcp45     = "./data/climatic_vars/selected/rcp45/",
                          biovar_rcp60     = "./data/climatic_vars/selected/rcp60/",
                          biovar_rcp85     = "./data/climatic_vars/selected/rcp85/",
                          cross_validation = 20)
  
  ###.............. Saving predictions
  writeRaster(result[["output_current"]], paste0("./data/outputs/XP1", sp_names[i],"_current.bil"), format = "EHdr")
  writeRaster(result[["output_rcp26"]],   paste0("./data/outputs/XP1", sp_names[i],"_rcp26.bil"), format = "EHdr")
  writeRaster(result[["output_rcp45"]],   paste0("./data/outputs/XP1", sp_names[i],"_rcp45.bil"), format = "EHdr")
  writeRaster(result[["output_rcp60"]],   paste0("./data/outputs/XP1", sp_names[i],"_rcp60.bil"), format = "EHdr")
  writeRaster(result[["output_rcp85"]],   paste0("./data/outputs/XP1", sp_names[i],"_rcp85.bil"), format = "EHdr")
  
  
  ###.............. Saving evaluation data
  write.table(result[["TPR_c"]],       paste0("./data/outputs/XP1", sp_names[i], "_TPR_current.txt"), sep = "\t", row.names = F)
  write.table(result[["Threshold_c"]], paste0("./data/outputs/XP1", sp_names[i], "_t_current.txt"),   sep = "\t", row.names = F)
  write.table(result[["Pred_area_c"]], paste0("./data/outputs/XP1", sp_names[i], "_d_current.txt"),   sep = "\t", row.names = F)
}
beep(8)


### Creating predictors variables to be used througout XP2:XP7
#   .....................................................................................

###.............. Presence Absence (PA)
all_thr <- 


bioclim <- stack(output_current[[1]])[[which(thr[,"bioclim"]>=0.7)]]
gower   <- stack(output_current[[2]])[[which(thr[,"gower"]>=0.7)]]
maxent  <- stack(output_current[[3]])[[which(thr[,"maxent"]>=0.7)]]
SVM     <- stack(output_current[[4]])[[which(thr[,"svm"]>=0.7)]]

pa_sep  <- 
pa_stk  <- 
  
  
###.............. Suitabilities (SUIT)
bioclim_val <- values(output_current[[1]])
gower_val   <- values(output_current[[2]])
maxent_val  <- values(output_current[[3]])
SVM_val     <- values(output_current[[4]])

bioclim_pad <- decostand(bioclim_val, "standardize", 2)
gower_pad   <- decostand(gower_val,   "standardize", 2)
maxent_pad  <- decostand(maxent_val,  "standardize", 2)
SVM_pad     <- decostand(SVM_val,     "standardize", 2)

suit_sep    <- 
  suit_stk    <- 

  
###.............. Sanving predictor variabels for XP2:XP7

# PA_SEP, PA_STK, SUIT_SEP, SUIT_SKT.

writeRaster(result[["pa_sep"]],   "./data/outputs/predictors/pa_sep.bil",   format = "EHdr")
writeRaster(result[["pa_stk"]],   "./data/outputs/predictors/pa_stk.bil",   format = "EHdr")
writeRaster(result[["suit_sep"]], "./data/outputs/predictors/suit_sep.bil", format = "EHdr")
writeRaster(result[["suit_stk"]], "./data/outputs/predictors/suit_stk.bil", format = "EHdr")

# *************************************************************************************** 
## 07. XP2                                          ----
# Scenario	         Present + Future - 4 RCPs
# Biotic Predictor	 SEP/PA - Plants XP1.2
# Variables       	 abiotic + SEP/PA -  (12 vars = 5  + 7 )
# Input           	 bee
# Output	           input x 4 predictive models
# Result  	         1 Ensemble

occur_thinned <- read.csv("./data/occurrences/occur_thinned.csv", sep = ",")
sp <- occur_thinned[occur_thinned[, 1] == "Lithurgus_huberi", ]
sp <- gsub("C[1-9]","", sp$SPEC)
sp_name <- unique(sp)
sp_name
for (i in 1:length(sp_name))
{
  ###.............. Running the modedelling experiment
  result <- multiple_ENMs(occurrence       = paste0("./data/occurrences/var_",  sp_name[i], ".txt"),
                          background       = paste0("./data/occurrences/back_", sp_name[i], ".txt"),
                          predictor        = "./data/outputs/predictors/",
                          biovar_current   = "./data/climatic_vars/selected/current/",
                          biovar_rcp26     = "./data/climatic_vars/selected/rcp26/",
                          biovar_rcp45     = "./data/climatic_vars/selected/rcp45/",
                          biovar_rcp60     = "./data/climatic_vars/selected/rcp60/",
                          biovar_rcp85     = "./data/climatic_vars/selected/rcp85/",
                          cross_validation = 20)
  
  ###.............. Saving predictions
  writeRaster(result[["output_current"]], paste0("./data/outputs/XP2", sp_name[i],"_current.bil"), format = "EHdr")
  writeRaster(result[["output_rcp26"]],   paste0("./data/outputs/XP2", sp_name[i],"_rcp26.bil"), format = "EHdr")
  writeRaster(result[["output_rcp45"]],   paste0("./data/outputs/XP2", sp_name[i],"_rcp45.bil"), format = "EHdr")
  writeRaster(result[["output_rcp60"]],   paste0("./data/outputs/XP2", sp_name[i],"_rcp60.bil"), format = "EHdr")
  writeRaster(result[["output_rcp85"]],   paste0("./data/outputs/XP2", sp_name[i],"_rcp85.bil"), format = "EHdr")
  
  ###.............. Saving evaluation data
  write.table(result[["TPR_c"]],       paste0("./data/outputs/XP2", sp_name[i], "_TPR_current.txt"), sep = "\t", row.names = F)
  write.table(result[["Threshold_c"]], paste0("./data/outputs/XP2", sp_name[i], "_t_current.txt"),   sep = "\t", row.names = F)
  write.table(result[["Pred_area_c"]], paste0("./data/outputs/XP2", sp_name[i], "_d_current.txt"),   sep = "\t", row.names = F)
  
}
beep(8)

# ***************************************************************************************
## 08. XP3                                          ----
# Scenario	         Present + Future - 4 RCPs
# Biotic Predictor	 SEP/SUIT - Plants XP1.2
# Variables       	 abiotic + SEP/SUIT -  (12 vars = 5  + 7 )
# Input           	 bee
# Output	           input x 4 predictive models
# Result  	         1 Ensemble

occur_thinned <- read.csv("./data/occurrences/occur_thinned.csv", sep = ",")
sp <- occur_thinned[occur_thinned[, 1] == "Lithurgus_huberi", ]
sp <- gsub("C[1-9]","", sp$SPEC)
sp_name <- unique(sp)
sp_name
for (i in 1:length(sp_name))
{
  ###.............. Running the modedelling experiment
  result <- multiple_ENMs(occurrence       = paste0("./data/occurrences/var_",  sp_name[i], ".txt"),
                          background       = paste0("./data/occurrences/back_", sp_name[i], ".txt"),
                          predictor        = "./data/outputs/predictors/",
                          biovar_current   = "./data/climatic_vars/selected/current/",
                          biovar_rcp26     = "./data/climatic_vars/selected/rcp26/",
                          biovar_rcp45     = "./data/climatic_vars/selected/rcp45/",
                          biovar_rcp60     = "./data/climatic_vars/selected/rcp60/",
                          biovar_rcp85     = "./data/climatic_vars/selected/rcp85/",
                          cross_validation = 20)
  
  ###.............. Saving predictions
  writeRaster(result[["output_current"]], paste0("./data/outputs/XP3", sp_name[i],"_current.bil"), format = "EHdr")
  writeRaster(result[["output_rcp26"]],   paste0("./data/outputs/XP3", sp_name[i],"_rcp26.bil"), format = "EHdr")
  writeRaster(result[["output_rcp45"]],   paste0("./data/outputs/XP3", sp_name[i],"_rcp45.bil"), format = "EHdr")
  writeRaster(result[["output_rcp60"]],   paste0("./data/outputs/XP3", sp_name[i],"_rcp60.bil"), format = "EHdr")
  writeRaster(result[["output_rcp85"]],   paste0("./data/outputs/XP3", sp_name[i],"_rcp85.bil"), format = "EHdr")
  
  ###.............. Saving evaluation data
  write.table(result[["TPR_c"]],       paste0("./data/outputs/XP3", sp_name[i], "_TPR_current.txt"), sep = "\t", row.names = F)
  write.table(result[["Threshold_c"]], paste0("./data/outputs/XP3", sp_name[i], "_t_current.txt"),   sep = "\t", row.names = F)
  write.table(result[["Pred_area_c"]], paste0("./data/outputs/XP3", sp_name[i], "_d_current.txt"),   sep = "\t", row.names = F)
  
}
beep(8)

# ***************************************************************************************
## 09. XP4                                          ----
# Scenario	         Present + Future - 4 RCPs
# Biotic Predictor	 STK/PA - Plants XP1.2
# Variables       	 abiotic + STK/PA -  (6 vars = 5  + 1 )
# Input           	 bee
# Output	           input x 4 predictive models
# Result  	         1 Ensemble

occur_thinned <- read.csv("./data/occurrences/occur_thinned.csv", sep = ",")
sp <- occur_thinned[occur_thinned[, 1] == "Lithurgus_huberi", ]
sp <- gsub("C[1-9]","", sp$SPEC)
sp_name <- unique(sp)
sp_name
for (i in 1:length(sp_name))
{
  ###.............. Running the modedelling experiment
  result <- multiple_ENMs(occurrence       = paste0("./data/occurrences/var_",  sp_name[i], ".txt"),
                          background       = paste0("./data/occurrences/back_", sp_name[i], ".txt"),
                          predictor        = "./data/outputs/predictors/",
                          biovar_current   = "./data/climatic_vars/selected/current/",
                          biovar_rcp26     = "./data/climatic_vars/selected/rcp26/",
                          biovar_rcp45     = "./data/climatic_vars/selected/rcp45/",
                          biovar_rcp60     = "./data/climatic_vars/selected/rcp60/",
                          biovar_rcp85     = "./data/climatic_vars/selected/rcp85/",
                          cross_validation = 20)
  
  ###.............. Saving predictions
  writeRaster(result[["output_current"]], paste0("./data/outputs/XP4", sp_name[i],"_current.bil"), format = "EHdr")
  writeRaster(result[["output_rcp26"]],   paste0("./data/outputs/XP4", sp_name[i],"_rcp26.bil"), format = "EHdr")
  writeRaster(result[["output_rcp45"]],   paste0("./data/outputs/XP4", sp_name[i],"_rcp45.bil"), format = "EHdr")
  writeRaster(result[["output_rcp60"]],   paste0("./data/outputs/XP4", sp_name[i],"_rcp60.bil"), format = "EHdr")
  writeRaster(result[["output_rcp85"]],   paste0("./data/outputs/XP4", sp_name[i],"_rcp85.bil"), format = "EHdr")
  
  ###.............. Saving evaluation data
  write.table(result[["TPR_c"]],       paste0("./data/outputs/XP4", sp_name[i], "_TPR_current.txt"), sep = "\t", row.names = F)
  write.table(result[["Threshold_c"]], paste0("./data/outputs/XP4", sp_name[i], "_t_current.txt"),   sep = "\t", row.names = F)
  write.table(result[["Pred_area_c"]], paste0("./data/outputs/XP4", sp_name[i], "_d_current.txt"),   sep = "\t", row.names = F)
  
}
beep(8)
# ***************************************************************************************
## 10. XP5                                          ----
# Scenario	         Present + Future - 4 RCPs
# Biotic Predictor	 STK/SUIT - Plants XP1.2
# Variables       	 abiotic + STK/SUIT -  (6 vars = 5  + 1 )
# Input           	 bee
# Output	           input x 4 predictive models
# Result  	         1 Ensemble

occur_thinned <- read.csv("./data/occurrences/occur_thinned.csv", sep = ",")
sp <- occur_thinned[occur_thinned[, 1] == "Lithurgus_huberi", ]
sp <- gsub("C[1-9]","", sp$SPEC)
sp_name <- unique(sp)
sp_name
for (i in 1:length(sp_name))
{
  ###.............. Running the modedelling experiment
  result <- multiple_ENMs(occurrence       = paste0("./data/occurrences/var_",  sp_name[i], ".txt"),
                          background       = paste0("./data/occurrences/back_", sp_name[i], ".txt"),
                          predictor        = "./data/outputs/predictors/",
                          biovar_current   = "./data/climatic_vars/selected/current/",
                          biovar_rcp26     = "./data/climatic_vars/selected/rcp26/",
                          biovar_rcp45     = "./data/climatic_vars/selected/rcp45/",
                          biovar_rcp60     = "./data/climatic_vars/selected/rcp60/",
                          biovar_rcp85     = "./data/climatic_vars/selected/rcp85/",
                          cross_validation = 20)
  
  ###.............. Saving predictions
  writeRaster(result[["output_current"]], paste0("./data/outputs/XP5", sp_name[i],"_current.bil"), format = "EHdr")
  writeRaster(result[["output_rcp26"]],   paste0("./data/outputs/XP5", sp_name[i],"_rcp26.bil"), format = "EHdr")
  writeRaster(result[["output_rcp45"]],   paste0("./data/outputs/XP5", sp_name[i],"_rcp45.bil"), format = "EHdr")
  writeRaster(result[["output_rcp60"]],   paste0("./data/outputs/XP5", sp_name[i],"_rcp60.bil"), format = "EHdr")
  writeRaster(result[["output_rcp85"]],   paste0("./data/outputs/XP5", sp_name[i],"_rcp85.bil"), format = "EHdr")
  
  ###.............. Saving evaluation data
  write.table(result[["TPR_c"]],       paste0("./data/outputs/XP5", sp_name[i], "_TPR_current.txt"), sep = "\t", row.names = F)
  write.table(result[["Threshold_c"]], paste0("./data/outputs/XP5", sp_name[i], "_t_current.txt"),   sep = "\t", row.names = F)
  write.table(result[["Pred_area_c"]], paste0("./data/outputs/XP5", sp_name[i], "_d_current.txt"),   sep = "\t", row.names = F)
  
}
beep(8)
# ***************************************************************************************
## 11. XP6                                          ----
# Scenario	         Present + Future - 4 RCPs
# Biotic Predictor	 Resource PA - Plants XP1.2
# Variables       	 abiotic + Resource PA -  (6 vars = 5  + 1 )
# Input           	 bee
# Output	           input x 4 predictive models
# Result  	         1 Ensemble

occur_thinned <- read.csv("./data/occurrences/occur_thinned.csv", sep = ",")
sp <- occur_thinned[occur_thinned[, 1] == "Lithurgus_huberi", ]
sp <- gsub("C[1-9]","", sp$SPEC)
sp_name <- unique(sp)
sp_name
for (i in 1:length(sp_name))
{
  ###.............. Running the modedelling experiment
  result <- multiple_ENMs(occurrence       = paste0("./data/occurrences/var_",  sp_name[i], ".txt"),
                          background       = paste0("./data/occurrences/back_", sp_name[i], ".txt"),
                          predictor        = "./data/outputs/predictors/",
                          biovar_current   = "./data/climatic_vars/selected/current/",
                          biovar_rcp26     = "./data/climatic_vars/selected/rcp26/",
                          biovar_rcp45     = "./data/climatic_vars/selected/rcp45/",
                          biovar_rcp60     = "./data/climatic_vars/selected/rcp60/",
                          biovar_rcp85     = "./data/climatic_vars/selected/rcp85/",
                          cross_validation = 20)
  
  ###.............. Saving predictions
  writeRaster(result[["output_current"]], paste0("./data/outputs/XP6", sp_name[i],"_current.bil"), format = "EHdr")
  writeRaster(result[["output_rcp26"]],   paste0("./data/outputs/XP6", sp_name[i],"_rcp26.bil"), format = "EHdr")
  writeRaster(result[["output_rcp45"]],   paste0("./data/outputs/XP6", sp_name[i],"_rcp45.bil"), format = "EHdr")
  writeRaster(result[["output_rcp60"]],   paste0("./data/outputs/XP6", sp_name[i],"_rcp60.bil"), format = "EHdr")
  writeRaster(result[["output_rcp85"]],   paste0("./data/outputs/XP6", sp_name[i],"_rcp85.bil"), format = "EHdr")
  
  ###.............. Saving evaluation data
  write.table(result[["TPR_c"]],       paste0("./data/outputs/XP6", sp_name[i], "_TPR_current.txt"), sep = "\t", row.names = F)
  write.table(result[["Threshold_c"]], paste0("./data/outputs/XP6", sp_name[i], "_t_current.txt"),   sep = "\t", row.names = F)
  write.table(result[["Pred_area_c"]], paste0("./data/outputs/XP6", sp_name[i], "_d_current.txt"),   sep = "\t", row.names = F)
  
}
beep(8)

# ***************************************************************************************
## 12. XP7                                          ----
# Scenario	         Present + Future - 4 RCPs
# Biotic Predictor	 Resource SUIT - Plants XP1.2
# Variables       	 abiotic + resource SUIT -  (6 vars = 5  + 1 )
# Input           	 bee
# Output	           input x 4 predictive models
# Result  	         1 Ensemble

occur_thinned <- read.csv("./data/occurrences/occur_thinned.csv", sep = ",")
sp <- occur_thinned[occur_thinned[, 1] == "Lithurgus_huberi", ]
sp <- gsub("C[1-9]","", sp$SPEC)
sp_name <- unique(sp)
sp_name
for (i in 1:length(sp_name))
{
  ###.............. Running the modedelling experiment
  result <- multiple_ENMs(occurrence       = paste0("./data/occurrences/var_",  sp_name[i], ".txt"),
                          background       = paste0("./data/occurrences/back_", sp_name[i], ".txt"),
                          predictor        = "./data/outputs/predictors/",
                          biovar_current   = "./data/climatic_vars/selected/current/",
                          biovar_rcp26     = "./data/climatic_vars/selected/rcp26/",
                          biovar_rcp45     = "./data/climatic_vars/selected/rcp45/",
                          biovar_rcp60     = "./data/climatic_vars/selected/rcp60/",
                          biovar_rcp85     = "./data/climatic_vars/selected/rcp85/",
                          cross_validation = 20)
  
  ###.............. Saving predictions
  writeRaster(result[["output_current"]], paste0("./data/outputs/XP7", sp_name[i],"_current.bil"), format = "EHdr")
  writeRaster(result[["output_rcp26"]],   paste0("./data/outputs/XP7", sp_name[i],"_rcp26.bil"), format = "EHdr")
  writeRaster(result[["output_rcp45"]],   paste0("./data/outputs/XP7", sp_name[i],"_rcp45.bil"), format = "EHdr")
  writeRaster(result[["output_rcp60"]],   paste0("./data/outputs/XP7", sp_name[i],"_rcp60.bil"), format = "EHdr")
  writeRaster(result[["output_rcp85"]],   paste0("./data/outputs/XP7", sp_name[i],"_rcp85.bil"), format = "EHdr")
  
  ###.............. Saving evaluation data
  write.table(result[["TPR_c"]],       paste0("./data/outputs/XP7", sp_name[i], "_TPR_current.txt"), sep = "\t", row.names = F)
  write.table(result[["Threshold_c"]], paste0("./data/outputs/XP7", sp_name[i], "_t_current.txt"),   sep = "\t", row.names = F)
  write.table(result[["Pred_area_c"]], paste0("./data/outputs/XP7", sp_name[i], "_d_current.txt"),   sep = "\t", row.names = F)
  
}
beep(8)


# ***************************************************************************************
## 13. Selecting XP from XP2:XP7                    ----

# AVALIAÇAO DA REPRESENTATIVIDADE	
# Anova de Medidas Repetidas (mesmos subconjuntos de XP2 a XP7)	
# 1	
# preditor: adequabilidade ponderada XP2:XP7	
# resposta: huberi (mesma resposta dividida nos mesmos subgrupos %treino/%teste no 6 XPs)	
# 2	
# preditor: tamanho de range  XP2:XP7	
# resposta: huberi	

## 14. Preparing analysis factors                   ----

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

### Getting factors for each species
# sp1
sp_names[1]
d_sp1    <- all_d[[1]]
suit_sp1 <- all_pad[[1]]

# sp2
sp_names[2]
d_sp2  <- all_d[[2]]

# sp3
sp_names[3]
d_sp3  <- all_d[[3]]

# sp4
sp_names[4]
d_sp4  <- all_d[[4]]

# sp5
sp_names[5]
d_sp5  <- all_d[[5]]

# sp6
sp_names[6]
d_sp6  <- all_d[[6]]

# sp7
sp_names[7]
d_sp7  <- all_d[[7]]

# sp8
sp_names[8]
d_sp8  <- all_d[[8]]

# ***************************************************************************************
## 15. Ensemble                                     ----






# ***************************************************************************************
## 16. Uncertainty Evaluation                       ----

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



# **** List of improvements to the scritp           ----

# 1. Implement validation by the checkerboards method.
# 2. Implement occurrence filtering at the ambiental space and compare with spThin (geographical space)
# 3. Transform maps in frequencies instead of suitabilities.
# 4. Impelement multi cores for running several models simultaneously.
# 5. Reduce code by implementing subfunctions, lopps.
# 6. Rewrite the code using tidy



