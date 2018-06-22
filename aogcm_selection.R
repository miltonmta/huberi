library (raster)
library (rgdal)

# Selecting AOGCMs - Cluster analysis####
#####
# Global Circulation Models (GCMs) selection through cluster analysis to reduce bias and improve uncertainty analysis.
# At worldclim.com at the resolution of 2.5min we selected the only the variables that apper in all the Representative Concetration Pathways projetions (RCP26, RCP45, RCP60, RCP80). The codes of the 11 GCMs utilized are: bc, cc, gs, hd, he, ip, mi, mr, mc, mg, no.


# [optional] if you prefer to rename all the variables so you have a more direct understanding when comparing then with each other, see below:

myPath <- './data/climatic_vars/60bi70/no60bi70'
fileList <- dir(path = myPath, pattern = '*.tif')  # list of file names, not including their paths
sapply(X = fileList, FUN = function(x) {
  file.rename(paste0(myPath, x),     # paste0() the path and old name
              paste0(myPath, 'bio', substring(x, first = 9))) })     # paste0() the path and new name

# substring('hello world', first = 7) returns 'orld' (starting at number 7, all following characters are returned) it should indicate the maximum charatcter in the folder you want to change. For instance if in some point of the folder you had "hello 1world" if you leave "7" all of your files would be called "world". Changing it for "8" would rename of them to "orld".

# A. RCP 26 ----
### read the variables ----

#importing all 19 variables from each one of the 11 GCMs of RCP 26


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


rcp_26 <- stack(bc_26, cc_26, gs_26, hd_26, he_26, ip_26, mi_26, mr_26, mc_26, mg_26, no_26)

# > print(raster(rcp_26))
# class       : RasterLayer 
# dimensions  : 3600, 8640, 31104000  (nrow, ncol, ncell)
# resolution  : 0.04166667, 0.04166667  (x, y)
# extent      : -180, 180, -60, 90  (xmin, xmax, ymin, ymax)
# coord. ref. : +proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs 


# crop raster - Neotropical (xmin, xmax, ymin, ymax) ----
# shape by Löwenberg-Neto, P. (2014) Neotropical region: a shapefile of Morrone's (2014) biogeographical regionalisation. Zootaxa, 3802(2): 300-300. 
# browseURL("http://purl.org/biochartis/neo2014shp")

print(raster(rcp_26))

e <- extent(-85,-30,-60,15)
rcp_26_e <- crop(rcp_26, e)
print(raster(rcp_26_e))
plot(rcp_26_e[[1]])
map(add=T)

shape_neo <- shapefile("./data/shape/Lowenberg_Neto_2014_shapefile/Lowenberg_Neto_2014.shp")
rcp_26_neot <- mask(crop(rcp_26_e, shape_neo), shape_neo) 

dim(rcp_26)
str(rcp_26)
dim(rcp_26_neo)
str(rcp_26_neo)
print(rcp_26)


#sd of the variables----



#map sd----

#identifying areas of high heterogeneity between models using the quartile coeff----

# absolute change, using thresholds----

# correlation between predictions and Hierarchical cluster analysis----


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

# crop raster - Neotropical (xmin, xmax, ymin, ymax) ----

#sd of the variables----

#map sd----

#identifying areas of high heterogeneity between models using the quartile coeff----

# absolute change, using thresholds----

# correlation between predictions and Hierarchical cluster analysis----


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

# crop raster - Neotropical (xmin, xmax, ymin, ymax) ----

#sd of the variables----

#map sd----

#identifying areas of high heterogeneity between models using the quartile coeff----

# absolute change, using thresholds----

# correlation between predictions and Hierarchical cluster analysis----

# D. RCP 85 ----
### read the variables ----

#importing all 19 variables from each one of the 11 chossen  worldclim GCMs of RCP 85


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
print(raster(cc_85))
print(raster(cc_45))



# crop raster - Neotropical (xmin, xmax, ymin, ymax) ----
e <- extent(-85,-30,-60,15)
rcp_85_e <- crop(rcp_85, e)
plot(rcp_26_e[[1]])
map(add=T)
print(raster(rcp_85_e))



#sd of the variables----

#map sd----

#identifying areas of high heterogeneity between models using the quartile coeff----

# absolute change, using thresholds----

# correlation between predictions and Hierarchical cluster analysis----
