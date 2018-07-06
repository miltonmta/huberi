# Selecting AOGCMs - Cluster analysis####

require(raster)
require(rgdal)
require(abind)
require(amap)
require(stats)


# This script has an index table. If you are in RStudio go to Code > Show Document Outline (shift + command/clrt + o)

# The directories with the variables utilized here can be viewed and downloaded from OneDrive:
browseURL(" https://1drv.ms/f/s!ApJZaitgpPr7gZtfS9n9mU9DDzXQMg")

# At worldclim.com at the resolution of 2.5min we selected  only the variables that apper simultaneously in all the Representative Concetration Pathways projetions (RCP26, RCP45, RCP60, RCP80). The codes of the 11 GCMs utilized are: bc, cc, gs, hd, he, ip, mi, mr, mc, mg, no.

## Wolrdclim GCM	code----
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
# model_names <- c("BCC-CSM1-1", "CCSM4", "GISS-EZ-R", "HadGEM2-AO", "HadGEM2-ES", "IPSL-CM5A-LR", "MIROC5", "MRI-CGCM3", "MIROC-ESM-CHEM", "MIROC-ESM", "NorESM1-M")# naming must be in the same reading order of the directory.


### Read the variables ----
# Tow options of functions for creating the RCP array

# "I've got the same cluster results running both functions over and over again (see uploaded plot "Cluster_RCP26_TinkerBell", "Cluster_RCP_TinkerBell", also tables of results). But this results are different form the ones I have found while running the long way, one model at the time ("Cluster_RCP26", "Cluster_RCP45")."

## Option 1 ----
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

rcp_26 <- ToothFairy( x = "./data/climatic_vars/26bi70/")
rcp_45 <- ToothFairy( x = "./data/climatic_vars/45bi70/")
rcp_60 <- ToothFairy( x = "./data/climatic_vars/60bi70/")
rcp_85 <- ToothFairy( x = "./data/climatic_vars/85bi70/")

## Option 2 ----
## Its ill-advised to grow complex objects like arrays in a loop. So here is another solution.

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
apn <- function(...) abind(..., along = 3)  # empty function for setting appending par to main function

# Models from RCP 26
x <- list.dirs("./data/climatic_vars/26bi70/", full.names = TRUE)[-1]
model_list <- lapply(x, tinker_bell)
rcp_26 <- do.call("apn", model_list)

# Models from RCP 45
x <- list.dirs("./data/climatic_vars/45bi70/", full.names = TRUE)[-1]
model_list <- lapply(x, tinker_bell)
rcp_45 <- do.call("apn", model_list)
rm(rcp_45_TinkerBell)

# Models from RCP 60
x <- list.dirs("./data/climatic_vars/60bi70/", full.names = TRUE)[-1]
model_list <- lapply(x, tinker_bell)
rcp_60 <- do.call("apn", model_list)

# Models from RCP 85
x <- list.dirs("./data/climatic_vars/85bi70/", full.names = TRUE)[-1]
model_list <- lapply(x, tinker_bell)
rcp_85 <- do.call("apn", model_list)


## plot variables

### Standard Deviation of the variables

# determining the Quartile Coefficient of Deviation (qcd)

### map sd

### identifying areas of high heterogeneity between models

### absolute change, using thresholds

### Correlation between predictions----
# library (amap)
model_names <- c("BCC-CSM1-1", "CCSM4", "GISS-EZ-R", "HadGEM2-AO", "HadGEM2-ES", "IPSL-CM5A-LR", "MIROC5", "MRI-CGCM3", "MIROC-ESM-CHEM", "MIROC-ESM", "NorESM1-M") # must be in the same order of the directories. 
hc <- list()
for (i in 1:19)
{
  raw_data <- t(rcp_85[ , i+2, ]) # get the variable data except the first two columms (lat, long)
  rownames (raw_data) <- model_names 
  cor_bio <- hcluster (raw_data, method = "correlation")
  # rect.hclust(raw_data, k=i, border = "gray") Erro: $ operator is invalid for atomic vectors
  plot (cor_bio)
  hc [[i]] <- cor_bio
  
}

# head (cor_bio) 
names (hc)<- c(paste ("BIO", c(1:19), sep=""))
par (las = 1)
for (i in 1:19)
{
  plot (hc[[i]])
  mtext (names(hc)[i], side = 1, line = 2)
}


### Euclidean clusters ----
# library (stats)

hc_2 <- list()
for (i in 1:19)
{
  raw_data <- t(rcp_85[ , i+2, ])
  rownames (raw_data) <- model_names 
  cor_bio <- hcluster (raw_data, method = "euclidean")
  hc_2[[i]] <- cor_bio
}

names (hc_2)<- c(paste ("BIO", c(1:19), sep=""))
par (las = 1)
for (i in 1:19)
{
  plot (hc_2[[i]]) 
  mtext (names(hc_2)[i], side= 1, line=2)

}

# dev.off()
plot (hc[[4]], 
      main = "Cluster Dendrogram\n(RCP 85)")


## Response grouping by k means

res_k_85<- NULL
for (i in 1:19){
  res<- cutree(hc[[i]], k=4) 
  res_k_85<- rbind (res_k_85, res)
}

# Write the RCP cluster and the results table
rownames (res_k_85)<- c(paste ("BIO", c(1:19), sep="")) 
hc_k_85<- hcluster (t(res_k_85), method="euclidean")
plot (hc_k_85, 
      main = "Cluster Dendrogram by K means\n(RCP 85 tinker bell)")
t(res_k_85)

## Response grouping by Height
res_h_85<- NULL
for (i in 1:19){
  res<- cutree(hc[[i]], h=0.2) 
  res_h_85<- rbind (res_h_85, res)
}

# Write the RCP cluster and the results table
rownames (res_h_85)<- c(paste ("BIO", c(1:19), sep="")) 
hc_h_85<- hcluster (t(res_h_85), method="euclidean")
plot (hc_h_85,
      main = "Cluster Dendrogram by Height\n(RCP  tinker bell)")
t(res_h_85)


# Tables of cluster results ####

#### rcp 26 ####

# > t(res_groups_k_26)
#                  BIO1 BIO2 BIO3 BIO4 BIO5 BIO6 BIO7 BIO8 BIO9 BIO10 BIO11
# BCC-CSM1-1        1    1    1    1    1    1    1    1    1     1     1
# CCSM4             1    2    1    1    2    2    2    2    1     2     2
# GISS-EZ-R         2    2    2    1    2    2    2    1    1     2     3
# HadGEM2-AO        2    2    2    2    1    2    3    1    2     3     3
# HadGEM2-ES        2    2    2    2    1    2    2    1    2     3     3
# IPSL-CM5A-LR      1    1    2    3    3    3    2    3    3     2     2
# MIROC-ESM-CHEM    3    3    3    4    4    4    4    4    4     4     4
# MIROC-ESM         4    3    4    2    2    4    4    2    4     4     4
# MIROC5            2    2    2    2    2    2    2    2    2     2     4
# MRI-CGCM3         1    4    1    2    1    1    2    1    1     2     2
# NorESM1-M         1    2    1    1    2    2    2    2    1     2     3
#                  BIO12 BIO13 BIO14 BIO15 BIO16 BIO17 BIO18 BIO19
# BCC-CSM1-1         1     1     1     1     1     1     1     1
# CCSM4              1     1     1     2     1     1     1     1
# GISS-EZ-R          1     1     1     2     1     1     2     2
# HadGEM2-AO         1     2     2     3     2     2     2     1
# HadGEM2-ES         1     2     2     3     2     2     1     1
# IPSL-CM5A-LR       2     1     3     1     2     3     1     3
# MIROC-ESM-CHEM     3     3     4     4     3     4     3     4
# MIROC-ESM          4     3     4     4     3     4     4     4
# MIROC5             1     4     1     2     3     1     1     2
# MRI-CGCM3          1     1     1     3     4     1     4     1
# NorESM1-M          1     4     1     2     1     1     1     1

# Result with the function ToothFairy
# > t(res_groups_k_26)
# BIO1 BIO2 BIO3 BIO4 BIO5 BIO6 BIO7 BIO8 BIO9 BIO10 BIO11
# BCC-CSM1-1        1    1    1    1    1    1    1    1    1     1     1
# CCSM4             1    2    1    1    2    2    2    2    1     2     2
# GISS-EZ-R         2    2    2    1    2    2    2    1    1     2     3
# HadGEM2-AO        2    2    2    2    1    2    3    1    2     3     3
# HadGEM2-ES        2    2    2    2    1    2    2    1    2     3     3
# IPSL-CM5A-LR      1    1    2    3    3    3    2    3    3     2     2
# MIROC5            2    2    2    2    2    2    2    2    2     2     4
# MRI-CGCM3         1    3    1    2    1    1    2    1    1     2     2
# MIROC-ESM-CHEM    3    4    3    4    4    4    4    4    4     4     4
# MIROC-ESM         4    4    4    2    2    4    4    2    4     4     4
# NorESM1-M         1    2    1    1    2    2    2    2    1     2     3
# BIO12 BIO13 BIO14 BIO15 BIO16 BIO17 BIO18 BIO19
# BCC-CSM1-1         1     1     1     1     1     1     1     1
# CCSM4              1     1     1     2     1     1     1     1
# GISS-EZ-R          1     1     1     2     1     1     2     2
# HadGEM2-AO         1     2     2     3     2     2     2     1
# HadGEM2-ES         1     2     2     3     2     2     1     1
# IPSL-CM5A-LR       2     1     3     1     2     3     1     3
# MIROC5             1     3     1     2     3     1     1     2
# MRI-CGCM3          1     1     1     3     4     1     3     1
# MIROC-ESM-CHEM     3     4     4     4     3     4     4     4
# MIROC-ESM          4     4     4     4     3     4     3     4
# NorESM1-M          1     3     1     2     1     1     1     1


# > t(res_groups_h_26)
#                  BIO1 BIO2 BIO3 BIO4 BIO5 BIO6 BIO7 BIO8 BIO9 BIO10 BIO11
# BCC-CSM1-1        1    1    1    1    1    1    1    1    1     1     1
# CCSM4             1    1    1    1    1    1    1    1    1     2     1
# GISS-EZ-R         1    1    1    1    1    1    1    1    1     2     1
# HadGEM2-AO        1    1    1    1    1    1    1    1    1     2     1
# HadGEM2-ES        1    1    1    1    1    1    1    1    1     2     1
# IPSL-CM5A-LR      1    1    1    1    1    1    1    1    1     2     1
# MIROC5            1    1    1    1    1    1    1    1    1     2     2
# MRI-CGCM3         1    1    1    1    1    1    1    1    1     2     2
# MIROC-ESM-CHEM    1    1    1    1    1    1    1    1    1     2     2
# MIROC-ESM         1    1    1    1    1    1    1    1    1     2     1
# NorESM1-M         1    1    1    1    1    1    1    1    1     2     1
#                 BIO12 BIO13 BIO14 BIO15 BIO16 BIO17 BIO18 BIO19
# BCC-CSM1-1         1     1     1     1     1     1     1     1
# CCSM4              1     1     1     1     1     1     1     1
# GISS-EZ-R          1     1     1     1     1     1     1     1
# HadGEM2-AO         1     1     1     1     1     1     1     1
# HadGEM2-ES         1     1     1     1     1     1     1     1
# IPSL-CM5A-LR       1     1     1     1     1     1     1     1
# MIROC5             1     1     1     1     1     1     1     1
# MRI-CGCM3          1     1     1     1     1     1     1     1
# MIROC-ESM-CHEM     1     1     1     1     1     1     1     1
# MIROC-ESM          1     1     1     1     1     1     1     1
# NorESM1-M          1     1     1     1     1     1     1     1


#### rcp 45 -----

# > t(res_groups_k_45)
#                  BIO1 BIO2 BIO3 BIO4 BIO5 BIO6 BIO7 BIO8 BIO9 BIO10 BIO11
# BCC-CSM1-1        1    1    1    1    1    1    1    1    1     1     1
# CCSM4             1    1    1    2    2    2    1    2    2     2     1
# GISS-EZ-R         2    2    2    1    1    2    1    1    2     1     2
# HadGEM2-AO        2    3    1    2    2    2    2    2    2     2     3
# HadGEM2-ES        2    3    1    2    2    2    2    2    2     2     3
# IPSL-CM5A-LR      1    1    1    3    3    1    3    3    1     3     2
# MIROC5            3    4    3    4    4    3    4    4    3     4     4
# MRI-CGCM3         4    4    4    4    4    3    4    4    3     4     4
# MIROC-ESM-CHEM    2    3    2    2    2    4    4    4    4     1     4
# MIROC-ESM         1    2    1    1    1    2    1    1    2     1     2
# NorESM1-M         1    1    1    2    2    2    1    2    2     2     1
#                  BIO12 BIO13 BIO14 BIO15 BIO16 BIO17 BIO18 BIO19
# BCC-CSM1-1         1     1     1     1     1     1     1     1
# CCSM4              1     1     1     1     1     1     1     2
# GISS-EZ-R          1     1     1     2     1     1     1     2
# HadGEM2-AO         1     2     2     1     2     2     2     3
# HadGEM2-ES         1     2     2     1     2     2     2     2
# IPSL-CM5A-LR       2     1     2     3     3     2     3     1
# MIROC5             3     3     3     4     1     3     3     4
# MRI-CGCM3          4     3     3     4     1     4     3     4
# MIROC-ESM-CHEM     1     4     4     1     1     1     1     2
# MIROC-ESM          1     1     1     2     4     1     4     2
# NorESM1-M          1     1     1     1     1     1     1     2

# > t(res_k_45) # feito com a função TinkerBell
# BIO1 BIO2 BIO3 BIO4 BIO5 BIO6 BIO7 BIO8 BIO9 BIO10 BIO11
# BCC-CSM1-1        1    1    1    1    1    1    1    1    1     1     1
# CCSM4             1    2    1    1    2    2    2    2    1     2     2
# GISS-EZ-R         2    2    2    1    2    2    2    1    1     2     3
# HadGEM2-AO        2    2    2    2    1    2    3    1    2     3     3
# HadGEM2-ES        2    2    2    2    1    2    2    1    2     3     3
# IPSL-CM5A-LR      1    1    2    3    3    3    2    3    3     2     2
# MIROC5            2    2    2    2    2    2    2    2    2     2     4
# MRI-CGCM3         1    3    1    2    1    1    2    1    1     2     2
# MIROC-ESM-CHEM    3    4    3    4    4    4    4    4    4     4     4
# MIROC-ESM         4    4    4    2    2    4    4    2    4     4     4
# NorESM1-M         1    2    1    1    2    2    2    2    1     2     3
# BIO12 BIO13 BIO14 BIO15 BIO16 BIO17 BIO18 BIO19
# BCC-CSM1-1         1     1     1     1     1     1     1     1
# CCSM4              1     1     1     2     1     1     1     1
# GISS-EZ-R          1     1     1     2     1     1     2     2
# HadGEM2-AO         1     2     2     3     2     2     2     1
# HadGEM2-ES         1     2     2     3     2     2     1     1
# IPSL-CM5A-LR       2     1     3     1     2     3     1     3
# MIROC5             1     3     1     2     3     1     1     2
# MRI-CGCM3          1     1     1     3     4     1     3     1
# MIROC-ESM-CHEM     3     4     4     4     3     4     4     4
# MIROC-ESM          4     4     4     4     3     4     3     4
# NorESM1-M          1     3     1     2     1     1     1     1
# > 

# > t(res_groups_h_45)
#                  BIO1 BIO2 BIO3 BIO4 BIO5 BIO6 BIO7 BIO8 BIO9 BIO10 BIO11
# BCC-CSM1-1        1    1    1    1    1    1    1    1    1     1     1
# CCSM4             1    1    1    1    1    1    1    1    1     1     1
# GISS-EZ-R         1    1    1    1    1    1    1    1    1     1     1
# HadGEM2-AO        1    1    1    1    1    1    1    1    1     1     1
# HadGEM2-ES        1    1    1    1    1    1    1    1    1     1     1
# IPSL-CM5A-LR      1    1    1    1    1    1    1    1    1     2     1
# MIROC5            1    1    1    1    1    1    1    1    1     3     2
# MRI-CGCM3         1    1    1    1    1    1    1    1    1     3     2
# MIROC-ESM-CHEM    1    1    1    1    1    1    1    1    1     1     2
# MIROC-ESM         1    1    1    1    1    1    1    1    1     1     1
# NorESM1-M         1    1    1    1    1    1    1    1    1     1     1
#                 BIO12 BIO13 BIO14 BIO15 BIO16 BIO17 BIO18 BIO19
# BCC-CSM1-1         1     1     1     1     1     1     1     1
# CCSM4              1     1     1     1     1     1     1     1
# GISS-EZ-R          1     1     1     1     1     1     1     1
# HadGEM2-AO         1     1     1     1     1     1     1     1
# HadGEM2-ES         1     1     1     1     1     1     1     1
# IPSL-CM5A-LR       1     1     1     1     1     1     1     1
# MIROC5             1     1     1     1     1     1     1     1
# MRI-CGCM3          1     1     1     1     1     1     1     1
# MIROC-ESM-CHEM     1     1     1     1     1     1     1     1
# MIROC-ESM          1     1     1     1     1     1     1     1
# NorESM1-M          1     1     1     1     1     1     1     1

#### rcp 60 ----

# > t(res_groups_k_60)
#                  BIO1 BIO2 BIO3 BIO4 BIO5 BIO6 BIO7 BIO8 BIO9 BIO10 BIO11
# BCC-CSM1-1        1    1    1    1    1    1    1    1    1     1     1
# CCSM4             1    2    1    1    1    1    2    2    1     2     1
# GISS-EZ-R         2    2    2    2    1    1    1    3    1     1     2
# HadGEM2-AO        3    3    3    1    2    1    1    1    2     2     1
# HadGEM2-ES        3    3    3    1    2    1    1    1    2     2     1
# IPSL-CM5A-LR      1    1    2    3    3    2    3    4    3     3     3
# MIROC5            4    4    4    4    4    3    4    3    4     4     4
# MRI-CGCM3         4    4    4    4    4    3    4    3    4     4     4
# MIROC-ESM-CHEM    3    2    2    4    1    4    4    3    4     2     4
# MIROC-ESM         3    2    2    4    1    1    1    3    1     2     2
# NorESM1-M         1    2    1    1    2    4    2    2    2     2     1
#                  BIO12 BIO13 BIO14 BIO15 BIO16 BIO17 BIO18 BIO19
# BCC-CSM1-1         1     1     1     1     1     1     1     1
# CCSM4              1     2     1     1     1     1     2     1
# GISS-EZ-R          1     1     1     1     2     1     3     1
# HadGEM2-AO         1     2     2     1     2     2     3     1
# HadGEM2-ES         1     2     2     1     2     2     2     1
# IPSL-CM5A-LR       2     1     3     2     3     3     4     2
# MIROC5             3     3     4     3     4     4     1     3
# MRI-CGCM3          4     3     4     3     4     4     1     4
# MIROC-ESM-CHEM     1     2     1     1     2     1     3     3
# MIROC-ESM          1     4     1     4     2     1     1     1
# NorESM1-M          1     2     1     1     1     1     2     1
# > 


# > t(res_groups_h_60)
#                  BIO1 BIO2 BIO3 BIO4 BIO5 BIO6 BIO7 BIO8 BIO9 BIO10 BIO11
# BCC-CSM1-1        1    1    1    1    1    1    1    1    1     1     1
# CCSM4             1    1    1    1    1    1    1    1    1     2     1
# GISS-EZ-R         1    1    1    1    1    1    1    1    1     1     2
# HadGEM2-AO        1    1    1    1    1    1    1    1    1     2     1
# HadGEM2-ES        1    1    1    1    1    1    1    1    1     2     1
# IPSL-CM5A-LR      1    1    1    1    1    1    1    1    1     3     2
# MIROC5            1    1    1    1    1    1    1    1    1     2     1
# MRI-CGCM3         1    1    1    1    1    1    1    1    1     2     1
# MIROC-ESM-CHEM    1    1    1    1    1    1    1    1    1     2     1
# MIROC-ESM         1    1    1    1    1    1    1    1    1     2     2
# NorESM1-M         1    1    1    1    1    1    1    1    1     2     1
#                  BIO12 BIO13 BIO14 BIO15 BIO16 BIO17 BIO18 BIO19
# BCC-CSM1-1         1     1     1     1     1     1     1     1
# CCSM4              1     1     1     1     1     1     1     1
# GISS-EZ-R          1     1     1     1     1     1     1     1
# HadGEM2-AO         1     1     1     1     1     1     1     1
# HadGEM2-ES         1     1     1     1     1     1     1     1
# IPSL-CM5A-LR       1     1     1     1     1     1     1     1
# MIROC5             1     1     1     1     1     1     1     1
# MRI-CGCM3          1     1     1     1     1     1     1     1
# MIROC-ESM-CHEM     1     1     1     1     1     1     1     1
# MIROC-ESM          1     1     1     1     1     1     1     1
# NorESM1-M          1     1     1     1     1     1     1     1
# 

#### rcp 85 ----

# > t(res_groups_k_85)
#                  BIO1 BIO2 BIO3 BIO4 BIO5 BIO6 BIO7 BIO8 BIO9 BIO10 BIO11
# BCC-CSM1-1        1    1    1    1    1    1    1    1    1     1     1
# CCSM4             1    1    1    2    1    1    2    1    1     1     1
# GISS-EZ-R         2    1    2    3    1    2    1    2    2     2     2
# HadGEM2-AO        3    1    2    2    2    1    1    1    2     3     3
# HadGEM2-ES        3    1    2    2    2    1    1    1    2     3     3
# IPSL-CM5A-LR      1    2    3    4    3    3    3    3    3     4     4
# MIROC5            4    3    4    1    1    4    4    4    4     1     1
# MRI-CGCM3         4    3    4    1    1    4    4    4    4     1     1
# MIROC-ESM-CHEM    1    1    1    1    1    4    2    1    4     1     1
# MIROC-ESM         1    4    1    1    1    1    1    2    2     1     2
# NorESM1-M         1    1    1    2    4    1    2    1    1     1     1
#                  BIO12 BIO13 BIO14 BIO15 BIO16 BIO17 BIO18 BIO19
# BCC-CSM1-1         1     1     1     1     1     1     1     1
# CCSM4              1     2     1     1     1     1     2     2
# GISS-EZ-R          2     1     2     1     2     1     2     2
# HadGEM2-AO         1     2     1     1     2     2     3     2
# HadGEM2-ES         1     2     1     1     2     2     3     2
# IPSL-CM5A-LR       3     1     2     2     3     3     1     3
# MIROC5             4     3     3     3     4     4     4     4
# MRI-CGCM3          4     3     3     3     4     4     4     4
# MIROC-ESM-CHEM     1     2     1     1     1     1     2     2
# MIROC-ESM          1     4     4     4     1     1     1     1
# NorESM1-M          1     2     1     1     1     1     2     2



# > t(res_groups_h_85)
#                  BIO1 BIO2 BIO3 BIO4 BIO5 BIO6 BIO7 BIO8 BIO9 BIO10 BIO11
# BCC-CSM1-1        1    1    1    1    1    1    1    1    1     1     1
# CCSM4             1    1    1    1    1    1    1    1    1     1     1
# GISS-EZ-R         1    1    1    1    1    1    1    2    1     2     2
# HadGEM2-AO        1    1    1    1    2    1    1    1    1     3     1
# HadGEM2-ES        1    1    1    1    2    1    1    1    1     3     1
# IPSL-CM5A-LR      1    1    1    1    3    1    1    2    1     4     2
# MIROC5            1    1    1    1    1    1    1    2    1     1     1
# MRI-CGCM3         1    1    1    1    1    1    1    2    1     1     1
# MIROC-ESM-CHEM    1    1    1    1    1    1    1    1    1     1     1
# MIROC-ESM         1    1    1    1    1    1    1    2    1     1     2
# NorESM1-M         1    1    1    1    4    1    1    1    1     1     1
#                  BIO12 BIO13 BIO14 BIO15 BIO16 BIO17 BIO18 BIO19
# BCC-CSM1-1         1     1     1     1     1     1     1     1
# CCSM4              1     1     1     1     1     1     1     1
# GISS-EZ-R          1     1     1     1     1     1     1     1
# HadGEM2-AO         1     1     1     1     1     1     1     1
# HadGEM2-ES         1     1     1     1     1     1     1     1
# IPSL-CM5A-LR       1     1     1     1     1     1     1     1
# MIROC5             1     1     1     1     1     1     1     1
# MRI-CGCM3          1     1     1     1     1     1     1     1
# MIROC-ESM-CHEM     1     1     1     1     1     1     1     1
# MIROC-ESM          1     1     1     1     1     1     1     1
# NorESM1-M          1     1     1     1     1     1     1     1


### ??? Construct clusters with the GCMs---- 
# Is't this what I just did above?


# # read correlation between models to 
# library (tree)
# library (rpart)

## correlation: Variables

## correlation models


rm(list = ls())

# A. RCP 26 ----
# read the variables manually----

bc_26 <- stack(list.files("./data/climatic_vars/26bi70/bc26bi70",  pattern = ".tif$", full.names = TRUE))

cc_26 <- stack(list.files("./data/climatic_vars/26bi70/cc26bi70",  pattern = ".tif$", full.names = TRUE))

gs_26 <- stack(list.files("./data/climatic_vars/26bi70/gs26bi70",  pattern = ".tif$", full.names = TRUE))

hd_26 <- stack(list.files("./data/climatic_vars/26bi70/hd26bi70",  pattern = ".tif$", full.names = TRUE))

he_26 <- stack(list.files("./data/climatic_vars/26bi70/he26bi70",  pattern = ".tif$", full.names = TRUE))

ip_26 <- stack(list.files("./data/climatic_vars/26bi70/ip26bi70",  pattern = ".tif$", full.names = TRUE))

mi_26 <- stack(list.files("./data/climatic_vars/26bi70/mi26bi70",  pattern = ".tif$", full.names = TRUE))

mr_26 <- stack(list.files("./data/climatic_vars/26bi70/mr26bi70",  pattern = ".tif$", full.names = TRUE))

mc_26 <- stack(list.files("./data/climatic_vars/26bi70/mc26bi70",  pattern = ".tif$", full.names = TRUE))

mg_26 <- stack(list.files("./data/climatic_vars/26bi70/mg26bi70",  pattern = ".tif$", full.names = TRUE))

no_26 <- stack(list.files("./data/climatic_vars/26bi70/no26bi70",  pattern = ".tif$", full.names = TRUE))


### cut raster to studied area (xmin, xmax, ymin, ymax) ----
e <- extent(-122, -18, -56, 14) 

bc_26_e <- crop (bc_26, e)
cc_26_e <- crop (cc_26, e)
gs_26_e <- crop (gs_26, e)
hd_26_e <- crop (hd_26, e)
he_26_e <- crop (he_26, e)
ip_26_e <- crop (ip_26, e)
mi_26_e <- crop (mi_26, e)
mr_26_e <- crop (mr_26, e)
mc_26_e <- crop (mc_26, e)
mg_26_e <- crop (mg_26, e)
no_26_e <- crop (no_26, e)

# rm(bc_26, cc_26, gs_26, hd_26, he_26, ip_26, mi_26, mr_26, mc_26, mg_26, no_26)
# 
# 
# rm(bc_26_e, cc_26_e, gs_26_e, hd_26_e, he_26_e, ip_26_e, mi_26_e, mr_26_e, mc_26_e, mg_26_e, no_26_e)
# 


## Masking with the shapefile
# shape by Löwenberg-Neto, P. (2014) Neotropical region: a shapefile of Morrone's (2014) biogeographical regionalisation. Zootaxa, 3802(2): 300-300. 
# browseURL("http://purl.org/biochartis/neo2014shp")

# another shape file 
# brouseURL("https://www.naturalearthdata.com/downloads/110m-physical-vectors/")

# shape_neot <- shapefile("./data/shape/Lowenberg_Neto_2014_shapefile/Lowenberg_Neto_2014.shp")
# rcp_26_neot <- mask(crop(rcp_26_e, shape_neot), shape_neot) 
# 
# dim(rcp_26)
# str(rcp_26)
# dim(rcp_26_neot)
# str(rcp_26_neot)
# print(raster(rcp_26_e))

### extracting values from raster ----

bc_26_val <- getValues(bc_26_e)
coord_bc_26 <- xyFromCell( bc_26_e, 1:ncell(bc_26_e)) # retrieve the variables coordinates.
bc_26_final <- cbind(coord_bc_26, bc_26_val) # merge coordinates with the extrated values data.
bc_26_final <- na.omit(bc_26_final)

cc_26_val <- values(cc_26_e)
coord_cc_26 <- xyFromCell( cc_26_e, 1:ncell(cc_26_e)) 
cc_26_final <- cbind(coord_cc_26, cc_26_val) 
cc_26_final <- na.omit(cc_26_final)

gs_26_val <- values(gs_26_e)
coord_gs_26 <- xyFromCell( gs_26_e, 1:ncell(gs_26_e))
gs_26_final <- cbind(coord_gs_26, gs_26_val) 
gs_26_final <- na.omit(gs_26_final)

hd_26_val <- values(hd_26_e)
coord_hd_26 <- xyFromCell( hd_26_e, 1:ncell(hd_26_e))
hd_26_final <- cbind(coord_hd_26, hd_26_val) 
hd_26_final <- na.omit(hd_26_final)

he_26_val <- values(he_26_e)
coord_he_26 <- xyFromCell( he_26_e, 1:ncell(he_26_e))
he_26_final <- cbind(coord_he_26, he_26_val) 
he_26_final <- na.omit(he_26_final)

ip_26_val <- values(ip_26_e)
coord_ip_26 <- xyFromCell( ip_26_e, 1:ncell(ip_26_e))
ip_26_final <- cbind(coord_ip_26, ip_26_val) 
ip_26_final <- na.omit(ip_26_final)

mi_26_val <- values(mi_26_e)
coord_mi_26 <- xyFromCell( mi_26_e, 1:ncell(mi_26_e))
mi_26_final <- cbind(coord_mi_26, mi_26_val) 
mi_26_final <- na.omit(mi_26_final)

mr_26_val <- values(mr_26_e)
coord_mr_26 <- xyFromCell( mr_26_e, 1:ncell(mr_26_e))
mr_26_final <- cbind(coord_mr_26, mr_26_val) 
mr_26_final <- na.omit(mr_26_final)

mc_26_val <- values(mc_26_e)
coord_mc_26 <- xyFromCell( mc_26_e, 1:ncell(mc_26_e))
mc_26_final <- cbind(coord_mc_26, mc_26_val) 
mc_26_final <- na.omit(mc_26_final)

mg_26_val <- values(mg_26_e)
coord_mg_26 <- xyFromCell( mg_26_e, 1:ncell(mg_26_e))
mg_26_final <- cbind(coord_mg_26, mg_26_val) 
mg_26_final <- na.omit(mg_26_final)

no_26_val <- values(no_26_e)
coord_no_26 <- xyFromCell( no_26_e, 1:ncell(no_26_e))
no_26_final <- cbind(coord_no_26, no_26_val) 
no_26_final <- na.omit(no_26_final)

## making array with abind 

rcp_26_no_function <- abind(bc_26_final, cc_26_final, gs_26_final, hd_26_final, he_26_final, ip_26_final, mi_26_final, mr_26_final, mc_26_final, mg_26_final, no_26_final, along = 3)

# B. RCP 45  ----
### read the variables manually ----


bc_45 <- stack(list.files("./data/climatic_vars/45bi70/bc45bi70",  pattern = ".tif$", full.names = TRUE))

cc_45 <- stack(list.files("./data/climatic_vars/45bi70/cc45bi70",  pattern = ".tif$", full.names = TRUE))

gs_45 <- stack(list.files("./data/climatic_vars/45bi70/gs45bi70",  pattern = ".tif$", full.names = TRUE))

hd_45 <- stack(list.files("./data/climatic_vars/45bi70/hd45bi70",  pattern = ".tif$", full.names = TRUE))

he_45 <- stack(list.files("./data/climatic_vars/45bi70/he45bi70",  pattern = ".tif$", full.names = TRUE))

ip_45 <- stack(list.files("./data/climatic_vars/45bi70/ip45bi70",  pattern = ".tif$", full.names = TRUE))

mi_45 <- stack(list.files("./data/climatic_vars/45bi70/mi45bi70",  pattern = ".tif$", full.names = TRUE))

mr_45 <- stack(list.files("./data/climatic_vars/45bi70/mr45bi70",  pattern = ".tif$", full.names = TRUE))

mc_45 <- stack(list.files("./data/climatic_vars/45bi70/mc45bi70",  pattern = ".tif$", full.names = TRUE))

mg_45 <- stack(list.files("./data/climatic_vars/45bi70/mg45bi70",  pattern = ".tif$", full.names = TRUE))

no_45 <- stack(list.files("./data/climatic_vars/45bi70/no45bi70",  pattern = ".tif$", full.names = TRUE))


### cut raster to studied area (xmin, xmax, ymin, ymax) ----
e <- extent(-122, -18, -56, 14) 

bc_45_e <- crop (bc_45, e)
cc_45_e <- crop (cc_45, e)
gs_45_e <- crop (gs_45, e)
hd_45_e <- crop (hd_45, e)
he_45_e <- crop (he_45, e)
ip_45_e <- crop (ip_45, e)
mi_45_e <- crop (mi_45, e)
mr_45_e <- crop (mr_45, e)
mc_45_e <- crop (mc_45, e)
mg_45_e <- crop (mg_45, e)
no_45_e <- crop (no_45, e)


### extracting values from raster ----

bc_45_val <- getValues(bc_45_e)
coord_bc_45 <- xyFromCell( bc_45_e, 1:ncell(bc_45_e)) 
bc_45_final <- cbind(coord_bc_45, bc_45_val) 
bc_45_final <- na.omit(bc_45_final)

cc_45_val <- values(cc_45_e)
coord_cc_45 <- xyFromCell( cc_45_e, 1:ncell(cc_45_e)) 
cc_45_final <- cbind(coord_cc_45, cc_45_val) 
cc_45_final <- na.omit(cc_45_final)

gs_45_val <- values(gs_45_e)
coord_gs_45 <- xyFromCell( gs_45_e, 1:ncell(gs_45_e))
gs_45_final <- cbind(coord_gs_45, gs_45_val) 
gs_45_final <- na.omit(gs_45_final)

hd_45_val <- values(hd_45_e)
coord_hd_45 <- xyFromCell( hd_45_e, 1:ncell(hd_45_e))
hd_45_final <- cbind(coord_hd_45, hd_45_val) 
hd_45_final <- na.omit(hd_45_final)

he_45_val <- values(he_45_e)
coord_he_45 <- xyFromCell( he_45_e, 1:ncell(he_45_e))
he_45_final <- cbind(coord_he_45, he_45_val) 
he_45_final <- na.omit(he_45_final)

ip_45_val <- values(ip_45_e)
coord_ip_45 <- xyFromCell( ip_45_e, 1:ncell(ip_45_e))
ip_45_final <- cbind(coord_ip_45, ip_45_val) 
ip_45_final <- na.omit(ip_45_final)

mi_45_val <- values(mi_45_e)
coord_mi_45 <- xyFromCell( mi_45_e, 1:ncell(mi_45_e))
mi_45_final <- cbind(coord_mi_45, mi_45_val) 
mi_45_final <- na.omit(mi_45_final)

mr_45_val <- values(mr_45_e)
coord_mr_45 <- xyFromCell( mr_45_e, 1:ncell(mr_45_e))
mr_45_final <- cbind(coord_mr_45, mr_45_val) 
mr_45_final <- na.omit(mr_45_final)

mc_45_val <- values(mc_45_e)
coord_mc_45 <- xyFromCell( mc_45_e, 1:ncell(mc_45_e))
mc_45_final <- cbind(coord_mc_45, mc_45_val) 
mc_45_final <- na.omit(mc_45_final)

mg_45_val <- values(mg_45_e)
coord_mg_45 <- xyFromCell( mg_45_e, 1:ncell(mg_45_e))
mg_45_final <- cbind(coord_mg_45, mg_45_val) 
mg_45_final <- na.omit(mg_45_final)

no_45_val <- values(no_45_e)
coord_no_45 <- xyFromCell( no_45_e, 1:ncell(no_45_e))
no_45_final <- cbind(coord_no_45, no_45_val) 
no_45_final <- na.omit(no_45_final)

## making array with abind 

rcp_45 <- abind(bc_45_final, cc_45_final, gs_45_final, hd_45_final, he_45_final, ip_45_final, mi_45_final, mr_45_final, mc_45_final, mg_45_final, no_45_final, along = 3)




# C. RCP 60 ----
### read the variables manually----


bc_60 <- stack(list.files("./data/climatic_vars/60bi70/bc60bi70",  pattern = ".tif$", full.names = TRUE))

cc_60 <- stack(list.files("./data/climatic_vars/60bi70/cc60bi70",  pattern = ".tif$", full.names = TRUE))

gs_60 <- stack(list.files("./data/climatic_vars/60bi70/gs60bi70",  pattern = ".tif$", full.names = TRUE))

hd_60 <- stack(list.files("./data/climatic_vars/60bi70/hd60bi70",  pattern = ".tif$", full.names = TRUE))

he_60 <- stack(list.files("./data/climatic_vars/60bi70/he60bi70",  pattern = ".tif$", full.names = TRUE))

ip_60 <- stack(list.files("./data/climatic_vars/60bi70/ip60bi70",  pattern = ".tif$", full.names = TRUE))

mi_60 <- stack(list.files("./data/climatic_vars/60bi70/mi60bi70",  pattern = ".tif$", full.names = TRUE))

mr_60 <- stack(list.files("./data/climatic_vars/60bi70/mr60bi70",  pattern = ".tif$", full.names = TRUE))

mc_60 <- stack(list.files("./data/climatic_vars/60bi70/mc60bi70",  pattern = ".tif$", full.names = TRUE))

mg_60 <- stack(list.files("./data/climatic_vars/60bi70/mg60bi70",  pattern = ".tif$", full.names = TRUE))

no_60 <- stack(list.files("./data/climatic_vars/60bi70/no60bi70",  pattern = ".tif$", full.names = TRUE))


### cut raster to studied area (xmin, xmax, ymin, ymax) ----
e <- extent(-122, -18, -56, 14) 

bc_60_e <- crop (bc_60, e)
cc_60_e <- crop (cc_60, e)
gs_60_e <- crop (gs_60, e)
hd_60_e <- crop (hd_60, e)
he_60_e <- crop (he_60, e)
ip_60_e <- crop (ip_60, e)
mi_60_e <- crop (mi_60, e)
mr_60_e <- crop (mr_60, e)
mc_60_e <- crop (mc_60, e)
mg_60_e <- crop (mg_60, e)
no_60_e <- crop (no_60, e)


### extracting values from raster ----

bc_60_val <- getValues(bc_60_e)
coord_bc_60 <- xyFromCell( bc_60_e, 1:ncell(bc_60_e)) 
bc_60_final <- cbind(coord_bc_60, bc_60_val) 
bc_60_final <- na.omit(bc_60_final)

cc_60_val <- values(cc_60_e)
coord_cc_60 <- xyFromCell( cc_60_e, 1:ncell(cc_60_e)) 
cc_60_final <- cbind(coord_cc_60, cc_60_val) 
cc_60_final <- na.omit(cc_60_final)

gs_60_val <- values(gs_60_e)
coord_gs_60 <- xyFromCell( gs_60_e, 1:ncell(gs_60_e))
gs_60_final <- cbind(coord_gs_60, gs_60_val) 
gs_60_final <- na.omit(gs_60_final)

hd_60_val <- values(hd_60_e)
coord_hd_60 <- xyFromCell( hd_60_e, 1:ncell(hd_60_e))
hd_60_final <- cbind(coord_hd_60, hd_60_val) 
hd_60_final <- na.omit(hd_60_final)

he_60_val <- values(he_60_e)
coord_he_60 <- xyFromCell( he_60_e, 1:ncell(he_60_e))
he_60_final <- cbind(coord_he_60, he_60_val) 
he_60_final <- na.omit(he_60_final)

ip_60_val <- values(ip_60_e)
coord_ip_60 <- xyFromCell( ip_60_e, 1:ncell(ip_60_e))
ip_60_final <- cbind(coord_ip_60, ip_60_val) 
ip_60_final <- na.omit(ip_60_final)

mi_60_val <- values(mi_60_e)
coord_mi_60 <- xyFromCell( mi_60_e, 1:ncell(mi_60_e))
mi_60_final <- cbind(coord_mi_60, mi_60_val) 
mi_60_final <- na.omit(mi_60_final)

mr_60_val <- values(mr_60_e)
coord_mr_60 <- xyFromCell( mr_60_e, 1:ncell(mr_60_e))
mr_60_final <- cbind(coord_mr_60, mr_60_val) 
mr_60_final <- na.omit(mr_60_final)

mc_60_val <- values(mc_60_e)
coord_mc_60 <- xyFromCell( mc_60_e, 1:ncell(mc_60_e))
mc_60_final <- cbind(coord_mc_60, mc_60_val) 
mc_60_final <- na.omit(mc_60_final)

mg_60_val <- values(mg_60_e)
coord_mg_60 <- xyFromCell( mg_60_e, 1:ncell(mg_60_e))
mg_60_final <- cbind(coord_mg_60, mg_60_val) 
mg_60_final <- na.omit(mg_60_final)

no_60_val <- values(no_60_e)
coord_no_60 <- xyFromCell( no_60_e, 1:ncell(no_60_e))
no_60_final <- cbind(coord_no_60, no_60_val) 
no_60_final <- na.omit(no_60_final)

## making array with abind 

rcp_60_hand_work <- abind(bc_60_final, cc_60_final, gs_60_final, hd_60_final, he_60_final, ip_60_final, mi_60_final, mr_60_final, mc_60_final, mg_60_final, no_60_final, along = 3)






# D. RCP 85 ----
### read the variables manually ----



bc_85 <- stack(list.files("./data/climatic_vars/85bi70/bc85bi70",  pattern = ".tif$", full.names = TRUE))

cc_85 <- stack(list.files("./data/climatic_vars/85bi70/cc85bi70",  pattern = ".tif$", full.names = TRUE))

gs_85 <- stack(list.files("./data/climatic_vars/85bi70/gs85bi70",  pattern = ".tif$", full.names = TRUE))

hd_85 <- stack(list.files("./data/climatic_vars/85bi70/hd85bi70",  pattern = ".tif$", full.names = TRUE))

he_85 <- stack(list.files("./data/climatic_vars/85bi70/he85bi70",  pattern = ".tif$", full.names = TRUE))

ip_85 <- stack(list.files("./data/climatic_vars/85bi70/ip85bi70",  pattern = ".tif$", full.names = TRUE))

mi_85 <- stack(list.files("./data/climatic_vars/85bi70/mi85bi70",  pattern = ".tif$", full.names = TRUE))

mr_85 <- stack(list.files("./data/climatic_vars/85bi70/mr85bi70",  pattern = ".tif$", full.names = TRUE))

mc_85 <- stack(list.files("./data/climatic_vars/85bi70/mc85bi70",  pattern = ".tif$", full.names = TRUE))

mg_85 <- stack(list.files("./data/climatic_vars/85bi70/mg85bi70",  pattern = ".tif$", full.names = TRUE))

no_85 <- stack(list.files("./data/climatic_vars/85bi70/no85bi70",  pattern = ".tif$", full.names = TRUE))


### cut raster to studied area (xmin, xmax, ymin, ymax) ----
e <- extent(-122, -18, -56, 14) 

bc_85_e <- crop (bc_85, e)
cc_85_e <- crop (cc_85, e)
gs_85_e <- crop (gs_85, e)
hd_85_e <- crop (hd_85, e)
he_85_e <- crop (he_85, e)
ip_85_e <- crop (ip_85, e)
mi_85_e <- crop (mi_85, e)
mr_85_e <- crop (mr_85, e)
mc_85_e <- crop (mc_85, e)
mg_85_e <- crop (mg_85, e)
no_85_e <- crop (no_85, e)


### extracting values from raster ----

bc_85_val <- getValues(bc_85_e)
coord_bc_85 <- xyFromCell( bc_85_e, 1:ncell(bc_85_e)) 
bc_85_final <- cbind(coord_bc_85, bc_85_val) 
bc_85_final <- na.omit(bc_85_final)

cc_85_val <- values(cc_85_e)
coord_cc_85 <- xyFromCell( cc_85_e, 1:ncell(cc_85_e)) 
cc_85_final <- cbind(coord_cc_85, cc_85_val) 
cc_85_final <- na.omit(cc_85_final)

gs_85_val <- values(gs_85_e)
coord_gs_85 <- xyFromCell( gs_85_e, 1:ncell(gs_85_e))
gs_85_final <- cbind(coord_gs_85, gs_85_val) 
gs_85_final <- na.omit(gs_85_final)

hd_85_val <- values(hd_85_e)
coord_hd_85 <- xyFromCell( hd_85_e, 1:ncell(hd_85_e))
hd_85_final <- cbind(coord_hd_85, hd_85_val) 
hd_85_final <- na.omit(hd_85_final)

he_85_val <- values(he_85_e)
coord_he_85 <- xyFromCell( he_85_e, 1:ncell(he_85_e))
he_85_final <- cbind(coord_he_85, he_85_val) 
he_85_final <- na.omit(he_85_final)

ip_85_val <- values(ip_85_e)
coord_ip_85 <- xyFromCell( ip_85_e, 1:ncell(ip_85_e))
ip_85_final <- cbind(coord_ip_85, ip_85_val) 
ip_85_final <- na.omit(ip_85_final)

mi_85_val <- values(mi_85_e)
coord_mi_85 <- xyFromCell( mi_85_e, 1:ncell(mi_85_e))
mi_85_final <- cbind(coord_mi_85, mi_85_val) 
mi_85_final <- na.omit(mi_85_final)

mr_85_val <- values(mr_85_e)
coord_mr_85 <- xyFromCell( mr_85_e, 1:ncell(mr_85_e))
mr_85_final <- cbind(coord_mr_85, mr_85_val) 
mr_85_final <- na.omit(mr_85_final)

mc_85_val <- values(mc_85_e)
coord_mc_85 <- xyFromCell( mc_85_e, 1:ncell(mc_85_e))
mc_85_final <- cbind(coord_mc_85, mc_85_val) 
mc_85_final <- na.omit(mc_85_final)

mg_85_val <- values(mg_85_e)
coord_mg_85 <- xyFromCell( mg_85_e, 1:ncell(mg_85_e))
mg_85_final <- cbind(coord_mg_85, mg_85_val) 
mg_85_final <- na.omit(mg_85_final)

no_85_val <- values(no_85_e)
coord_no_85 <- xyFromCell( no_85_e, 1:ncell(no_85_e))
no_85_final <- cbind(coord_no_85, no_85_val) 
no_85_final <- na.omit(no_85_final)

## making array with abind 

rcp_85 <- abind(bc_85_final, cc_85_final, gs_85_final, hd_85_final, he_85_final, ip_85_final, mi_85_final, mr_85_final, mc_85_final, mg_85_final, no_85_final, along = 3)




