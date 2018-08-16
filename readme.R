
# --- General Info                                  ----

### This script has an index table. If you are in RStudio go to Code > Show Document Outline (shift + command / clrt + o)

### Worldclim variables for current conditions (~1960-1990), by the spatial resolution of 2.5min. 
# browseURL("http://www.worldclim.org/current")

### The models projected for 2070 (average for 2061-2080) were obtain at "worldclim.com" by the spatial resolution of 2.5min (0.04º or ≈ 4.4km). We selected the variables that appear simultaneously in all the representative concentration pathways scnarios (RCP26, RCP45, RCP60, RCP80). The codes of the 11 GCMs utilized are: bc, cc, gs, hd, he, ip, mi, mr, mc, mg, no.
# browseURL("http://www.worldclim.org/cmip5_2.5m")

###.............. wolrdclim GCM codes
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

###.............. selected AOGCMs
# based on the cluster analysis we imported only the selected variables at each RCP scenario.

# RCP 26: CCSM4(CC), HADGEM2-AO(HD),   IPSL-CMSA-LR(IP), MIROC-ESM(MR) 
# RCP 45: CCSM4(CC), HADGEM2-AO(HD),   IPSL-CMSA-LR(IP), MIROC-ESM(MR) 
# RCP 60: CCSM4(CC), IPSL-CMSA-LR(IP), MRI-CGCM3(MG),    MIROC-ESM(MR)
# RCP 85: CCSM4(CC), IPSL-CMSA-LR(IP), MIROC5(MC),       MIROC-ESM(MR)

### for maintanining comparability... 

# RCP 26: CCSM4(CC),                   IPSL-CMSA-LR(IP), MIROC-ESM(MR) 
# RCP 45: CCSM4(CC),                   IPSL-CMSA-LR(IP), MIROC-ESM(MR) 
# RCP 60: CCSM4(CC), IPSL-CMSA-LR(IP),                   MIROC-ESM(MR)
# RCP 85: CCSM4(CC), IPSL-CMSA-LR(IP),                   MIROC-ESM(MR)

###.............. bioclimatic variables descriptions

#BIO01 = Annual Mean Temperature
#BIO02 = Mean Diurnal Range (Mean of monthly (max temp - min temp))
#BIO03 = Isothermality (BIO2/BIO7) (* 100)
#BIO04 = Temperature Seasonality (standard deviation *100)
#BIO05 = Max Temperature of Warmest Month
#BIO06 = Min Temperature of Coldest Month
#BIO07 = Temperature Annual Range (BIO5-BIO6)
#BIO08 = Mean Temperature of Wettest Quarter
#BIO09 = Mean Temperature of Driest Quarter
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


### --- Maxent/rJava troubleshoot                     ----

### MaxEnt is available as a standalone Java program. Dismo has a function 'maxent' that communicates with this program. To use it you must first download the program from:
# browseURL("http://www.cs.princeton.edu/~schapire/maxent/")
### Put the file 'maxent.jar' in the 'java' folder of the 'dismo' package. That is the folder returned by system.file("java", package="dismo"). You need MaxEnt version 3.3.3b or higher.

### While loagind the package rJava, if you have the error image not found, follow the steps:
### 1. Check if JDK is installed in you machine. If not:
# browseURL("http://www.oracle.com/technetwork/java/javase/downloads/jdk8-downloads-2133151.html")

### 2. Checking if the jar file is present. 
jar <- paste(system.file(package="dismo"), "/java/maxent.jar/", sep='')
jar
file.exists(jar)

### 3. Use dyn.load for ataching the libjvm.dylib file before running the `rJava` Package.
# dyn.load('/Library/Java/JavaVirtualMachines/jdk1.8.0_131.jdk/Contents/Home/jre/lib/server/libjvm.dylib')

### 4. If you are in MacOS, a permanent solution is to link libjvm.dylib to /usr/local.lib at terminal:
# sudo ln -f -s $(/usr/libexec/java_home)/jre/lib/server/libjvm.dylib /usr/local/lib
# browseURL ("https://stackoverflow.com/questions/30738974/rjava-load-error-in-rstudio-r-after-upgrading-to-osx-yosemite")