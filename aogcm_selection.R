# Selecting AOGCMs - Cluster analysis####

require(raster)
require(rgdal)
require(abind)
require(amap)
require(stats)


# This script has an index table. If you are in RStudio go to Code > Show Document Outline (shift + command / clrt + o)

# The directories with the data utilized and the ones outputted here can be downloaded from the following OneDrive repositorium:
browseURL("https://1drv.ms/f/s!ApJZaitgpPr7gZtfS9n9mU9DDzXQMg")

# The models projected for 2070 (average for 2061-2080) were obtain at "worldclim.com" by the spatial resolution of 2.5min (0.04º or ≈ 4.4km). We selected the variables that appear simultaneously in all the representative concentration pathways scnarios (RCP26, RCP45, RCP60, RCP80). The codes of the 11 GCMs utilized are: bc, cc, gs, hd, he, ip, mi, mr, mc, mg, no.
browseURL("http://www.worldclim.org/cmip5_2.5m")

## Wolrdclim GCM	code----

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
# model_names <- c("BCC-CSM1-1", "CCSM4", "GISS-EZ-R", "HadGEM2-AO", "HadGEM2-ES", "IPSL-CM5A-LR", "MIROC5", "MRI-CGCM3", "MIROC-ESM-CHEM", "MIROC-ESM", "NorESM1-M")# naming must be in the same reading order of the origin directory.

# model_names <- c("BC", "CC", "GS", "HD", "HE", "IP", "MC", "MG", "MI", "MR", "NO") # must be in the same order of the directories.

### Read the variables ----

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

rcp_26 <- tooth_fairy( x = "./data/climatic_vars/26bi70/")
rcp_45 <- tooth_fairy( x = "./data/climatic_vars/45bi70/")
rcp_60 <- tooth_fairy( x = "./data/climatic_vars/60bi70/")
rcp_85 <- tooth_fairy( x = "./data/climatic_vars/85bi70/")



## plot variables

### Standard Deviation of the variables-----

# determining the Quartile Coefficient of Deviation (qcd)

### Map sd -----

### Identifying areas of high heterogeneity between models -----

### Absolute change, using thresholds ----

### Correlation between predictions----
# library (amap)
model_names <- c("BCC-CSM1-1", "CCSM4", "GISS-EZ-R", "HadGEM2-AO", "HadGEM2-ES", "IPSL-CM5A-LR", "MIROC5", "MRI-CGCM3", "MIROC-ESM-CHEM", "MIROC-ESM", "NorESM1-M") # must be in the same order of the directories.

# model_names <- c("BC", "CC", "GS", "HD", "HE", "IP", "MC", "MG", "MI", "MR", "NO") # must be in the same order of the directories.

# model_names <- c("ALTER-CCSM4", "ALTER-NorESM1-M", "BCC-CSM1-1", "GISS-EZ-R", "HadGEM2-AO", "HadGEM2-ES", "IPSL-CM5A-LR", "MIROC5", "MRI-CGCM3", "MIROC-ESM-CHEM", "MIROC-ESM") # must be in the same order of th

hc <- list()
for (i in 1:19)
{
  raw_data <- t(rcp_26[ , i+2, ]) # get the variable data except the first two columms (lat, long)
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
  raw_data <- t(rcp_26[ , i+2, ])
  rownames (raw_data) <- model_names 
  cor_bio <- hcluster (raw_data, method = "euclidean")
  hc_2[[i]] <- cor_bio
}

names (hc_2) <- c(paste ("BIO", c(1:19), sep = ""))
par (las = 1)
for (i in 1:19)
{
  plot (hc_2[[i]]) 
  mtext (names(hc_2)[i], side = 1, line = 2)

}

dev.off()
plot (hc[[4]],
      hang = -1,
      main = "Cluster Dendrogram\n(RCP 26)")


## Response grouping by k means

res_k_26 <- NULL
for (i in 1:19){
  res <- cutree(hc[[i]], k = 4) 
  res_k_26 <- rbind (res_k_26, res)
}

# Write the RCP cluster and the results table
rownames (res_k_26) <- c(paste ("BIO", c(1:19), sep = "")) 
hc_k_26 <- hcluster (t(res_k_26), method = "euclidean")
<<<<<<< HEAD
hcd <- as.dendrogram((hc_k_26))
nodePar <- list(lab.cex = 0.9, pch = c(NA,19), cex = 0.7, col = "blue")
par (oma = c(0, 0, 0, 4))
plot (hcd,
=======
dev.off()

nodePar <- list(lab.cex = 0.9, pch = c(NA,19), cex = 0.7, col = "blue")
par(oma = c(0, 0, 0, 4))
plot (as.dendrogram(hc_k_26),
>>>>>>> 2a2296c5b970e0b19669f81d1362e2a49b46a345
      horiz   = TRUE,
      nodePar = nodePar,
      edgePar = list(col = 1:1, lwd = 2:1),
      xlim    = c(14, 0),
      xlab    = "Height",
      main    = "Cluster Dendrogram by K means\n(RCP 26)")
t(res_k_26)

## Response grouping by Height
res_h_26 <- NULL
for (i in 1:19){
  res <- cutree(hc[[i]], h = 0.2) 
  res_h_26 <- rbind (res_h_26, res)
}

# Write the RCP cluster and the results table
rownames (res_h_26) <- c(paste ("BIO", c(1:19), sep = "")) 
hc_h_26 <- hcluster (t(res_h_26), method = "euclidean")
plot (hc_h_26,
      hang = -1,
      # cex  = 0.6,
      main = "Cluster Dendrogram by Height\n(RCP 26)")
t(res_h_26)


# Tables of cluster results ####

#### rcp 26 ####

# > t(res_groups_k_26)
#                 BIO1 BIO2 BIO3 BIO4 BIO5 BIO6 BIO7 BIO8 BIO9 BIO10 BIO11
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
#                 BIO12 BIO13 BIO14 BIO15 BIO16 BIO17 BIO18 BIO19
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


#### rcp 45 -----

# > t(res_k_45)
#                 BIO1 BIO2 BIO3 BIO4 BIO5 BIO6 BIO7 BIO8 BIO9 BIO10 BIO11
# BCC-CSM1-1        1    1    1    1    1    1    1    1    1     1     1
# CCSM4             1    1    1    2    2    2    1    2    2     2     1
# GISS-EZ-R         2    2    2    1    1    2    1    1    2     1     2
# HadGEM2-AO        2    3    1    2    2    2    2    2    2     2     3
# HadGEM2-ES        2    3    1    2    2    2    2    2    2     2     3
# IPSL-CM5A-LR      1    1    1    3    3    1    3    3    1     3     2
# MIROC5            2    3    2    2    2    3    4    4    3     1     4
# MRI-CGCM3         1    2    1    1    1    2    1    1    2     1     2
# MIROC-ESM-CHEM    3    4    3    4    4    4    4    4    4     4     4
# MIROC-ESM         4    4    4    4    4    4    4    4    4     4     4
# NorESM1-M         1    1    1    2    2    2    1    2    2     2     1
#                 BIO12 BIO13 BIO14 BIO15 BIO16 BIO17 BIO18 BIO19
# BCC-CSM1-1         1     1     1     1     1     1     1     1
# CCSM4              1     1     1     1     1     1     1     2
# GISS-EZ-R          1     1     1     2     1     1     1     2
# HadGEM2-AO         1     2     2     1     2     2     2     3
# HadGEM2-ES         1     2     2     1     2     2     2     2
# IPSL-CM5A-LR       2     1     2     3     3     2     3     1
# MIROC5             1     3     3     1     1     1     1     2
# MRI-CGCM3          1     1     1     2     4     1     4     2
# MIROC-ESM-CHEM     3     4     4     4     1     3     3     4
# MIROC-ESM          4     4     4     4     1     4     3     4
# NorESM1-M          1     1     1     1     1     1     1     2
# > 
  
#### rcp 60 ----

# > t(res_k_60)
#                 BIO1 BIO2 BIO3 BIO4 BIO5 BIO6 BIO7 BIO8 BIO9 BIO10 BIO11
# BCC-CSM1-1        1    1    1    1    1    1    1    1    1     1     1
# CCSM4             1    2    1    1    1    1    2    2    1     2     1
# GISS-EZ-R         2    2    2    2    1    1    1    3    1     1     2
# HadGEM2-AO        3    3    3    1    2    1    1    1    2     2     1
# HadGEM2-ES        3    3    3    1    2    1    1    1    2     2     1
# IPSL-CM5A-LR      1    1    2    3    3    2    3    4    3     3     3
# MIROC5            3    2    2    4    1    3    4    3    4     2     4
# MRI-CGCM3         3    2    2    4    1    1    1    3    1     2     2
# MIROC-ESM-CHEM    4    4    4    4    4    4    4    3    4     4     4
# MIROC-ESM         4    4    4    4    4    4    4    3    4     4     4
# NorESM1-M         1    2    1    1    2    3    2    2    2     2     1
#                 BIO12 BIO13 BIO14 BIO15 BIO16 BIO17 BIO18 BIO19
# BCC-CSM1-1         1     1     1     1     1     1     1     1
# CCSM4              1     2     1     1     1     1     2     1
# GISS-EZ-R          1     1     1     1     2     1     3     1
# HadGEM2-AO         1     2     2     1     2     2     3     1
# HadGEM2-ES         1     2     2     1     2     2     2     1
# IPSL-CM5A-LR       2     1     3     2     3     3     4     2
# MIROC5             1     2     1     1     2     1     3     3
# MRI-CGCM3          1     3     1     3     2     1     1     1
# MIROC-ESM-CHEM     3     4     4     4     4     4     1     3
# MIROC-ESM          4     4     4     4     4     4     1     4
# NorESM1-M          1     2     1     1     1     1     2     1
# > 

#### rcp 85 ----

# > t(res_k_85)
#                 BIO1 BIO2 BIO3 BIO4 BIO5 BIO6 BIO7 BIO8 BIO9 BIO10 BIO11
# BCC-CSM1-1        1    1    1    1    1    1    1    1    1     1     1
# CCSM4             1    1    1    2    1    1    2    1    1     1     1
# GISS-EZ-R         2    1    2    3    1    2    1    2    2     2     2
# HadGEM2-AO        3    1    2    2    2    1    1    1    2     3     3
# HadGEM2-ES        3    1    2    2    2    1    1    1    2     3     3
# IPSL-CM5A-LR      1    2    3    4    3    3    3    3    3     4     4
# MIROC5            1    1    1    1    1    4    2    1    4     1     1
# MRI-CGCM3         1    3    1    1    1    1    1    2    2     1     2
# MIROC-ESM-CHEM    4    4    4    1    1    4    4    4    4     1     1
# MIROC-ESM         4    4    4    1    1    4    4    4    4     1     1
# NorESM1-M         1    1    1    2    4    1    2    1    1     1     1
#                 BIO12 BIO13 BIO14 BIO15 BIO16 BIO17 BIO18 BIO19
# BCC-CSM1-1         1     1     1     1     1     1     1     1
# CCSM4              1     2     1     1     1     1     2     2
# GISS-EZ-R          2     1     2     1     2     1     2     2
# HadGEM2-AO         1     2     1     1     2     2     3     2
# HadGEM2-ES         1     2     1     1     2     2     3     2
# IPSL-CM5A-LR       3     1     2     2     3     3     1     3
# MIROC5             1     2     1     1     1     1     2     2
# MRI-CGCM3          1     3     3     3     1     1     1     1
# MIROC-ESM-CHEM     4     4     4     4     4     4     4     4
# MIROC-ESM          4     4     4     4     4     4     4     4
# NorESM1-M          1     2     1     1     1     1     2     2
# > 

### ??? Construct clusters with the GCMs---- 
# Is't this what I just did above?


# # read correlation between models to 
# library (tree)
# library (rpart)

## correlation: Variables

## correlation models
