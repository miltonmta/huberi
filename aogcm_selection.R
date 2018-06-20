# Selecting AOGCMs - Cluster analysis####
#####
# Global Circulation Models (GCMs) selection through cluster analysis to reduce bias and improve uncertainty analysis.
# At worldclim.com at the resolution of 2.5min we selected the only the variables that apper in all the Representative Concetration Pathways projetions (RCP26, RCP45, RCP60, RCP80). The codes of the GCMs utilized are: bc, cc, gs, hd, he, ip, mi, mr, mc, mg, no.

library (raster)
library (rgdal)

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

## plot variables
dim(rcp_26)
str(rcp_26)

colores<- colorRampPalette (c("darkblue", "blue", "lightblue", 
                              "white", "salmon", "red"))
par (mar=c(0,0,0,0))
plot (rcp_26, axes=F, box=F, col=colores(100), legend=F)



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

#sd of the variables----

#map sd----

#identifying areas of high heterogeneity between models using the quartile coeff----

# absolute change, using thresholds----

# correlation between predictions and Hierarchical cluster analysis----

# D. RCP 85 ----
### read the variables ----

#importing all 19 variables from each one of the 11 chossen  worldclim GCMs of RCP 85


bc_85 <- stack(list.files("./data/climatic_vars/85bi70/bcc-csm1-1(bc)",  pattern = ".asc$", full.names = TRUE))

cc_85 <- stack(list.files("./data/climatic_vars/85bi70/ccsm4(cc)",  pattern = ".asc$", full.names = TRUE))

gs_85 <- stack(list.files("./data/climatic_vars/85bi70/giss-e2-r(gs)",  pattern = ".asc$", full.names = TRUE))

hd_85 <- stack(list.files("./data/climatic_vars/85bi70/hadgem2-ao(hd)",  pattern = ".asc$", full.names = TRUE))

he_85 <- stack(list.files("./data/climatic_vars/85bi70/hadgem2-es(he)",  pattern = ".asc$", full.names = TRUE))

ip_85 <- stack(list.files("./data/climatic_vars/85bi70/ipsl-cm5a-lr(ip)",  pattern = ".asc$", full.names = TRUE))

mr_85 <- stack(list.files("./data/climatic_vars/85bi70/miroc-esm(mr)",  pattern = ".asc$", full.names = TRUE))

mi_85 <- stack(list.files("./data/climatic_vars/85bi70/miroc-esm-chem(mi)",  pattern = ".asc$", full.names = TRUE))

mc_85 <- stack(list.files("./data/climatic_vars/85bi70/miroc5(mc)",  pattern = ".asc$", full.names = TRUE))

mg_85 <- stack(list.files("./data/climatic_vars/85bi70/mri-cgcm3(mg)",  pattern = ".asc$", full.names = TRUE))

no_85 <- stack(list.files("./data/climatic_vars/85bi70/noresm1-m(no)",  pattern = ".asc$", full.names = TRUE))

rcp_85 <- stack(bc_85, cc_85, gs_85, hd_85, he_85, ip_85, mi_85, mr_85, mc_85, mg_85, no_85)


#sd of the variables----

#map sd----

#identifying areas of high heterogeneity between models using the quartile coeff----

# absolute change, using thresholds----

# correlation between predictions and Hierarchical cluster analysis----
