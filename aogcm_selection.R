
# Selecting AOGCMs - Cluster analysis####

require(raster)
require(rgdal)
require(abind)
require(amap)
require(stats)

# This script has an index table. If you are in RStudio go to Code > Show Document Outline (shift + command/clrt + o)

# Global Circulation Models (GCMs) selection through cluster analysis to reduce bias and improve uncertainty analysis.
# At worldclim.com at the resolution of 2.5min we selected the only the variables that apper in all the Representative Concetration Pathways projetions (RCP26, RCP45, RCP60, RCP80). The codes of the 11 GCMs utilized are: bc, cc, gs, hd, he, ip, mi, mr, mc, mg, no.

#####
## [optional] renaming the variables from all GCMs at all four RCPs (bio1:bio19):

# myPath <- './data/climatic_vars/60bi70/no60bi70'
# fileList <- dir(path = myPath, pattern = '*.tif')  # list of file names, not including their paths
# sapply(X = fileList, FUN = function(x) {
#   file.rename(paste0(myPath, x),     # paste0() the path and old name
#               paste0(myPath, 'bio', substring(x, first = 9))) })     # paste0() the path and new name
# substring('smaug', first = 2) returns 'aug' (starting at number 2, all following characters are returned) 

# A. RCP 26 ----
### read the variables ----

##importing all 19 variables from each one of the 11 GCMs of RCP 26

get()


bc_26 <- stack(list.files("./data/climatic_vars/26bi70/bc26bi70",  pattern = ".tif$", full.names = TRUE))

cc_26 <- stack(list.files("./data/climatic_vars/26bi70/cc26bi70",  pattern = ".tif$", full.names = TRUE))

gs_26 <- stack(list.files("./data/climatic_vars/26bi70/gs26bi70",  pattern = ".tif$", full.names = TRUE))

hd_26 <- stack(list.files("./data/climatic_vars/26bi70/hd26bi70",  pattern = ".tif$", full.names = TRUE))

he_26 <- stack(list.files("./data/climatic_vars/26bi70/he26bi70",  pattern = ".tif$", full.names = TRUE))

ip_26 <- stack(list.files("./data/climatic_vars/26bi70/ip26bi70",  pattern = ".tif$", full.names = TRUE))

mr_26 <- stack(list.files("./data/climatic_vars/26bi70/mr26bi70",  pattern = ".tif$", full.names = TRUE))

mi_26 <- stack(list.files("./data/climatic_vars/26bi70/mi26bi70",  pattern = ".tif$", full.names = TRUE))

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

# rm(bc_26_e, cc_26_e, gs_26_e, hd_26_e, he_26_e, ip_26_e, mi_26_e, mr_26_e, mc_26_e, mg_26_e, no_26_e)

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

rcp_26 <- abind(bc_26_final, cc_26_final, gs_26_final, hd_26_final, he_26_final, ip_26_final, mi_26_final, mr_26_final, mc_26_final, mg_26_final, no_26_final, along = 3)

print(rcp_26)
## plot variables

### Standard Deviation of the variables----

# determining the Quartile Coefficient of Deviation (qcd)

### map sd -----

### identifying areas of high heterogeneity between models---- 

### absolute change, using thresholds----

### correlation between predictions and Hierarchical cluster analysis----

## Correlation between predictions
# library (amap)

head (bc_26_final)
# > head (bc_26_final)
# x        y bio1 bio10 bio11 bio12 bio13 bio14 bio15 bio16 bio17
# [1,] -91.39583 13.97917  284   294   273  1455   408     1   101   837    11
# [2,] -91.35417 13.97917  284   294   273  1472   409     1   101   844    11
# [3,] -91.31250 13.97917  284   294   274  1485   410     1   100   849    12
# [4,] -91.27083 13.97917  284   294   274  1502   410     2   100   855    11
# [5,] -91.22917 13.97917  285   295   274  1510   410     2    99   856    11
# [6,] -91.18750 13.97917  285   295   275  1511   408     2    99   855    11
# bio18 bio19 bio2 bio3 bio4 bio5 bio6 bio7 bio8 bio9
# [1,]   364    11  112   72  817  357  202  155  284  273
# [2,]   368    11  112   72  812  357  203  154  284  273
# [3,]   369    12  111   72  806  357  203  154  284  274
# [4,]   374    11  111   72  804  356  203  153  284  274
# [5,]   374    11  110   72  799  357  205  152  284  274
# [6,]   374    11  110   72  791  357  205  151  285  275
print(raster (bc_26_final))
# > print(raster (bc_26_final))
# class       : RasterLayer 
# dimensions  : 912559, 21, 19163739  (nrow, ncol, ncell)
# resolution  : 0.04761905, 1.09582e-06  (x, y)
# extent      : 0, 1, 0, 1  (xmin, xmax, ymin, ymax)
# coord. ref. : NA 
# data source : in memory
# names       : layer 
# values      : -210, 11012  (min, max)

# Here i'm attempting to  creat a cluster for only one aogcm. Them to this for all the 11. I'm not sure how to incorporate all the aogcm model names here. Since I imported them one by one, I do not have any object from witch I could extract it from.

hc <- list()
for (i in 1:19)
{
  row_data <- bc_26_final[ nrow = i, ncol = -2 ] # get only the varible data but not the lat and long (first tow col)
  cor_bio <- hcluster (row_data, method = "correlation")
  plot (cor_bio)
  hc [[i]] <- cor_bio
  
}

head (cor_bio)
names (hc)<- c(paste ("BIO", c(1:19), sep=""))
par (las = 1)
for (i in 1:19)
{
  plot (hc[[i]])
  mtext (names(hc)[i], side = 1, line = 2)
}

## Euclidean clusters
# library (stats)

hc_2 <- list()
for (i in 1:19)
{
  row_data <- bc_26_final[ nrow = i, ncol = -2 ]
  clust_bio <- hcluster (row_data, method = "euclidean")
  hc2[[i]] <- clust_bio
}

names (hc2)<- c(paste ("BIO", c(1:19), sep=""))
par (las = 1)
for (i in 1:19)
{
  plot (hc2[[i]]) 
  mtext (names(hc2)[i], side= 1, line=2)
}

dev.off()
plot (hc[[4]])





# B. RCP 45 ----
### read the variables ----



bc_45 <- stack(list.files("./data/climatic_vars/45bi70/bc45bi70",  pattern = ".tif$", full.names = TRUE))

cc_45 <- stack(list.files("./data/climatic_vars/45bi70/cc45bi70",  pattern = ".tif$", full.names = TRUE))

gs_45 <- stack(list.files("./data/climatic_vars/45bi70/gs45bi70",  pattern = ".tif$", full.names = TRUE))

hd_45 <- stack(list.files("./data/climatic_vars/45bi70/hd45bi70",  pattern = ".tif$", full.names = TRUE))

he_45 <- stack(list.files("./data/climatic_vars/45bi70/he45bi70",  pattern = ".tif$", full.names = TRUE))

ip_45 <- stack(list.files("./data/climatic_vars/45bi70/ip45bi70",  pattern = ".tif$", full.names = TRUE))

mr_45 <- stack(list.files("./data/climatic_vars/45bi70/mr45bi70",  pattern = ".tif$", full.names = TRUE))

mi_45 <- stack(list.files("./data/climatic_vars/45bi70/mi45bi70",  pattern = ".tif$", full.names = TRUE))

mc_45 <- stack(list.files("./data/climatic_vars/45bi70/mc45bi70",  pattern = ".tif$", full.names = TRUE))

mg_45 <- stack(list.files("./data/climatic_vars/45bi70/mg45bi70",  pattern = ".tif$", full.names = TRUE))

no_45 <- stack(list.files("./data/climatic_vars/45bi70/no45bi70",  pattern = ".tif$", full.names = TRUE))


rcp_45 <- stack(bc_45, cc_45, gs_45, hd_45, he_45, ip_45, mi_45, mr_45, mc_45, mg_45, no_45)

### cut raster to studied area (xmin, xmax, ymin, ymax) ----
print(raster(rcp_45))

# e <- extent(-122, -18, -56, 14) 
rcp_26_e <- crop(rcp_45, e)
print(raster(rcp_45_e))
plot(rcp_45_e[[1]])
map(add=T)

?array

#extraindo valores do raster
ccsm.0k.val <- values(ccsm.0k.ASr)
ccsm.0k.val[1:5,]
nrow(ccsm.0k.val)

coord.AS <- xyFromCell(ccsm.0k.ASr, 1:ncell(ccsm.0k.ASr))
coord.AS[1:5,]
nrow(coord.AS)

ccsm.0k.ASm <- cbind(coord.AS, ccsm.0k.val)
ccsm.0k.ASm[1:5,]
nrow(ccsm.0k.ASm)
ccsm.0k.ASm <- na.omit(ccsm.0k.ASm)
nrow(ccsm.0k.ASm)

## Salving the stack of the GCMs at RCP45 cutted to studied area
writeRaste("./data/climatic_vars/45bi70/", rcp_45_e, "stack_rcp45_extent" ,format = "raster")

### Standard Deviation of the variables----

### identifying areas of high heterogeneity between models---- 

### absolute change, using thresholds----

### correlation between predictions and Hierarchical cluster analysis----



# C. RCP 60 ----
### read the variables ----


bc_60 <- stack(list.files("./data/climatic_vars/60bi70/bc60bi70",  pattern = ".tif$", full.names = TRUE))

cc_60 <- stack(list.files("./data/climatic_vars/60bi70/cc60bi70",  pattern = ".tif$", full.names = TRUE))

gs_60 <- stack(list.files("./data/climatic_vars/60bi70/gs60bi70",  pattern = ".tif$", full.names = TRUE))

hd_60 <- stack(list.files("./data/climatic_vars/60bi70/hd60bi70",  pattern = ".tif$", full.names = TRUE))

he_60 <- stack(list.files("./data/climatic_vars/60bi70/he60bi70",  pattern = ".tif$", full.names = TRUE))

ip_60 <- stack(list.files("./data/climatic_vars/60bi70/ip60bi70",  pattern = ".tif$", full.names = TRUE))

mr_60 <- stack(list.files("./data/climatic_vars/60bi70/mr60bi70",  pattern = ".tif$", full.names = TRUE))

mi_60 <- stack(list.files("./data/climatic_vars/60bi70/mi60bi70",  pattern = ".tif$", full.names = TRUE))

mc_60 <- stack(list.files("./data/climatic_vars/60bi70/mc60bi70",  pattern = ".tif$", full.names = TRUE))

mg_60 <- stack(list.files("./data/climatic_vars/60bi70/mg60bi70",  pattern = ".tif$", full.names = TRUE))

no_60 <- stack(list.files("./data/climatic_vars/60bi70/no60bi70",  pattern = ".tif$", full.names = TRUE))

rcp_60 <- stack(bc_60, cc_60, gs_60, hd_60, he_60, ip_60, mi_60, mr_60, mc_60, mg_60, no_60)


### crop raster to studied area (xmin, xmax, ymin, ymax) ----
print(raster(rcp_60))

# e <- extent(-122, -18, -56, 14) 
rcp_60_e <- crop(rcp_60, e)
print(raster(rcp_60_e))
plot(rcp_60_e[[1]])
map(add=T)

## Salving the stack of the GCMs at RCP60 cutted to studied area
writeRaste("./data/climatic_vars/60bi70/", rcp_60_e, "stack_rcp60_extent" ,format = "raster")

### Standard Deviation of the variables----

### identifying areas of high heterogeneity between models---- 

### absolute change, using thresholds----

### correlation between predictions----

### Hierarchical cluster analysis----

# D. RCP 85 ----
### read the variables ----

# importing all 19 variables from each one of the 11 chossen  worldclim GCMs of RCP 85


bc_85 <- stack(list.files("./data/climatic_vars/85bi70/bc85bi70",  pattern = ".asc$", full.names = TRUE))

cc_85 <- stack(list.files("./data/climatic_vars/85bi70/cc85bi70",  pattern = ".asc$", full.names = TRUE))

gs_85 <- stack(list.files("./data/climatic_vars/85bi70/gs85bi70",  pattern = ".asc$", full.names = TRUE))

hd_85 <- stack(list.files("./data/climatic_vars/85bi70/hd85bi70",  pattern = ".asc$", full.names = TRUE))

he_85 <- stack(list.files("./data/climatic_vars/85bi70/he85bi70",  pattern = ".asc$", full.names = TRUE))

ip_85 <- stack(list.files("./data/climatic_vars/85bi70/ip85bi70",  pattern = ".asc$", full.names = TRUE))

mr_85 <- stack(list.files("./data/climatic_vars/85bi70/mr85bi70",  pattern = ".asc$", full.names = TRUE))

mi_85 <- stack(list.files("./data/climatic_vars/85bi70/mi85bi70",  pattern = ".asc$", full.names = TRUE))

mc_85 <- stack(list.files("./data/climatic_vars/85bi70/mc85bi70",  pattern = ".asc$", full.names = TRUE))

mg_85 <- stack(list.files("./data/climatic_vars/85bi70/mg85bi70",  pattern = ".asc$", full.names = TRUE))

no_85 <- stack(list.files("./data/climatic_vars/85bi70/no85bi70",  pattern = ".asc$", full.names = TRUE))

rcp_85 <- stack(bc_85, cc_85, gs_85, hd_85, he_85, ip_85, mi_85, mr_85, mc_85, mg_85, no_85)

print(raster(cc_26))
print(raster(bc_85))
print(raster(cc_45))
bc_85_e <- crop(bc_85, e)


### crop raster to studied area (xmin, xmax, ymin, ymax) ----
e <- extent(-122, -18, -56, 14)
rcp_85_e <- crop(rcp_85, e)
plot(rcp_26_e[[1]])
map(add=T)
print(raster(rcp_85_e))

## Salving the stack of the GCMs at RCP85 cutted to studied area
writeRaste("./data/climatic_vars/85bi70/", rcp_85_e, "stack_rcp85_extent" ,format = "raster")

### Standard Deviation of the variables----

### identifying areas of high heterogeneity between models---- 

### absolute change, using thresholds----

### correlation between predictions----

###Hierarchical cluster analysis----
