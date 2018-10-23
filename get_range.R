get.range <- function(biovar_current,
                      sp_name,
                      pout,
                      thr,
                      cross_val,
                      AOGCMs)
{
  current <- stack(list.files(biovar_current,  pattern = ".grd$", full.names = TRUE))
  coords <- na.omit(cbind(xyFromCell(current, 1:ncell(current)), values(current)))[,1:2]
  thr <- read.table(paste0(thr, sp_name, "t_current.txt"), sep = "\t", h = T)
  n_cells <- nrow(na.omit(values(current)))
  
  #... reading predictions for present
  bioclim_c   <- stack(paste0(pout, sp_name,"_bioclim_c.tif"))
  SVM_c       <- stack(paste0(pout, sp_name,"_SVM_c.tif"))
  maxent_c    <- stack(paste0(pout, sp_name,"_maxent_c.tif"))
  
  #... creating objects to save results
  bioclim_pres <- maxent_pres <- SVM_pres <- NULL
  bioclim_26_r <- maxent_26_r <- SVM_26_r <- NULL
  bioclim_45_r <- maxent_45_r <- SVM_45_r <- NULL
  bioclim_60_r <- maxent_60_r <- SVM_60_r <- NULL
  bioclim_85_r <- maxent_85_r <- SVM_85_r <- NULL

  for (i in 1:cross_val)
  {
    #... reading predictions for futures
    bioclim_26  <- stack(paste0(pout, sp_name, i, "_bioclim_rcp26.tif"))
    SVM_26      <- stack(paste0(pout, sp_name, i, "_SVM_rcp26.tif"))
    maxent_26   <- stack(paste0(pout, sp_name, i, "_maxent_rcp26.tif"))
    
    bioclim_45  <- stack(paste0(pout, sp_name, i, "_bioclim_rcp45.tif"))
    SVM_45      <- stack(paste0(pout, sp_name, i, "_SVM_rcp45.tif"))
    maxent_45   <- stack(paste0(pout, sp_name, i, "_maxent_rcp45.tif"))
    
    bioclim_60  <- stack(paste0(pout, sp_name, i, "_bioclim_rcp60.tif"))
    SVM_60      <- stack(paste0(pout, sp_name, i, "_SVM_rcp60.tif"))
    maxent_60   <- stack(paste0(pout, sp_name, i, "_maxent_rcp60.tif"))
    
    bioclim_85  <- stack(paste0(pout, sp_name, i, "_bioclim_rcp85.tif"))
    SVM_85      <- stack(paste0(pout, sp_name, i, "_SVM_rcp85.tif"))
    maxent_85   <- stack(paste0(pout, sp_name, i, "_maxent_rcp85.tif"))
    
    #... calculating ranges
    #pres
    
    val_bioclim <- na.omit(values(bioclim_c[[i]]))
    PAb <- ifelse(val_bioclim >= thr$bioclim[i], 1, 0)
    # bioclim_pres <- c(bioclim_pres, sum(PAb)/n_cells) # range prop number of cells of the studied grid.
    bioclim_pres <- c(bioclim_pres, sum(PAb)/1000) # range = number of cells with presence (sum 1s)
    
    val_maxent <- na.omit(values(maxent_c[[i]]))
    PAm <- ifelse(val_maxent >= thr$maxent[i], 1, 0)
    # maxent_pres <- c(maxent_pres, sum(PAm)/n_cells)
    maxent_pres <- c(maxent_pres, sum(PAm)/1000)
    
    val_SVM <- na.omit(values(SVM_c[[i]]))
    PAs <- ifelse(val_SVM >= thr$SVM[i], 1, 0)
    # SVM_pres <- c(SVM_pres, sum(PAs)/n_cells)
    SVM_pres <- c(SVM_pres, sum(PAs)/1000)
    
    #fut
    bioclim_fut26 <- maxent_fut26 <- SVM_fut26 <- NULL
    bioclim_fut45 <- maxent_fut45 <- SVM_fut45 <- NULL
    bioclim_fut60 <- maxent_fut60 <- SVM_fut60 <- NULL
    bioclim_fut85 <- maxent_fut85 <- SVM_fut85 <- NULL
    
    for (j in 1:AOGCMs)
    {
      #RCP26
      val_bioclim <- na.omit(values(bioclim_26[[j]]))
      PAb <- ifelse(val_bioclim >= thr$bioclim[i], 1, 0)
      bioclim_fut26 <- cbind(bioclim_fut26, sum(PAb)/1000) 
      # bioclim_fut26 <- cbind(bioclim_fut26, sum(PAb)/n_cells) 
      
      val_maxent <- na.omit(values(maxent_26[[j]]))
      PAm <- ifelse(val_maxent >= thr$maxent[i], 1, 0)
      maxent_fut26 <- cbind(maxent_fut26, sum(PAm)/1000)
      # maxent_fut26 <- cbind(maxent_fut26, sum(PAm)/n_cells)
      
      val_SVM <- na.omit(values(SVM_26[[j]]))
      PAs <- ifelse(val_SVM >= thr$SVM[i], 1, 0)
      SVM_fut26 <- cbind(SVM_fut26, sum(PAs)/1000)
      # SVM_fut26 <- cbind(SVM_fut26, sum(PAs)/n_cells)
      
      #RCP45
      val_bioclim <- na.omit(values(bioclim_45[[j]]))
      PAb <- ifelse(val_bioclim >= thr$bioclim[i], 1, 0)
      bioclim_fut45 <- cbind(bioclim_fut45, sum(PAb)/1000)
      # bioclim_fut45 <- cbind(bioclim_fut45, sum(PAb)/n_cells)
      
      val_maxent <- na.omit(values(maxent_45[[j]]))
      PAm <- ifelse(val_maxent >= thr$maxent[i], 1, 0)
      maxent_fut45 <- cbind(maxent_fut45, sum(PAm)/1000)
      # maxent_fut45 <- cbind(maxent_fut45, sum(PAm)/n_cells)
      
      val_SVM <- na.omit(values(SVM_45[[j]]))
      PAs <- ifelse(val_SVM >= thr$SVM[i], 1, 0)
      SVM_fut45 <- cbind(SVM_fut45, sum(PAs)/1000)
      # SVM_fut45 <- cbind(SVM_fut45, sum(PAs)/n_cells)
      
      #RCP60
      val_bioclim <- na.omit(values(bioclim_60[[j]]))
      PAb <- ifelse(val_bioclim >= thr$bioclim[i], 1, 0)
      bioclim_fut60 <- cbind(bioclim_fut60, sum(PAb)/1000)
      # bioclim_fut60 <- cbind(bioclim_fut60, sum(PAb)/n_cells)
      
      val_maxent <- na.omit(values(maxent_60[[j]]))
      PAm <- ifelse(val_maxent >= thr$maxent[i], 1, 0)
      maxent_fut60 <- cbind(maxent_fut60, sum(PAm)/1000)
      # maxent_fut60 <- cbind(maxent_fut60, sum(PAm)/n_cells)
      
      val_SVM <- na.omit(values(SVM_60[[j]]))
      PAs <- ifelse(val_SVM >= thr$SVM[i], 1, 0)
      SVM_fut60 <- cbind(SVM_fut60, sum(PAs)/1000)
      # SVM_fut60 <- cbind(SVM_fut60, sum(PAs)/n_cells)
      
      #RCP85
      val_bioclim <- na.omit(values(bioclim_85[[j]]))
      PAb <- ifelse(val_bioclim >= thr$bioclim[i], 1, 0)
      bioclim_fut85 <- cbind(bioclim_fut85, sum(PAb)/1000)
      # bioclim_fut85 <- cbind(bioclim_fut85, sum(PAb)/n_cells)
      
      val_maxent <- na.omit(values(maxent_85[[j]]))
      PAm <- ifelse(val_maxent >= thr$maxent[i], 1, 0)
      maxent_fut85 <- cbind(maxent_fut85, sum(PAm)/1000)
      # maxent_fut85 <- cbind(maxent_fut85, sum(PAm)/n_cells)
      
      val_SVM <- na.omit(values(SVM_85[[j]]))
      PAs <- ifelse(val_SVM >= thr$SVM[i], 1, 0)
      SVM_fut85 <- cbind(SVM_fut85, sum(PAs)/1000)
      # SVM_fut85 <- cbind(SVM_fut85, sum(PAs)/n_cells)
    }
    
    bioclim_26_r <- rbind(bioclim_26_r, bioclim_fut26)
    bioclim_45_r <- rbind(bioclim_45_r, bioclim_fut45)
    bioclim_60_r <- rbind(bioclim_60_r, bioclim_fut60)
    bioclim_85_r <- rbind(bioclim_85_r, bioclim_fut85)
    
    maxent_26_r <- rbind(maxent_26_r, maxent_fut26)
    maxent_45_r <- rbind(maxent_45_r, maxent_fut45)
    maxent_60_r <- rbind(maxent_60_r, maxent_fut60)
    maxent_85_r <- rbind(maxent_85_r, maxent_fut85)
    
    SVM_26_r <- rbind(SVM_26_r, SVM_fut26)
    SVM_45_r <- rbind(SVM_45_r, SVM_fut45)
    SVM_60_r <- rbind(SVM_60_r, SVM_fut60)
    SVM_85_r <- rbind(SVM_85_r, SVM_fut85)
  }
  pres  <- data.frame(bioclim = bioclim_pres, maxent = maxent_pres, SVM = SVM_pres)
  fut26 <- data.frame(bioclim = bioclim_26_r, maxent = maxent_26_r, SVM = SVM_26_r)
  fut45 <- data.frame(bioclim = bioclim_45_r, maxent = maxent_45_r, SVM = SVM_45_r)
  fut60 <- data.frame(bioclim = bioclim_60_r, maxent = maxent_60_r, SVM = SVM_60_r)
  fut85 <- data.frame(bioclim = bioclim_85_r, maxent = maxent_85_r, SVM = SVM_85_r)
  
  return(list("pres" = pres,
              "fut26" = fut26,
              "fut45" = fut45,
              "fut60" = fut60,
              "fut85" = fut85))
}
