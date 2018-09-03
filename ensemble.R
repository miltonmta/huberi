ensemble <- function(dir)
{
  ### Reading predictions             -----
  ###......................................
  bioclim_c     <- stack(list.files(dir, pattern = "bioclim_c.tif$",     full.names = TRUE))
  gower_c       <- stack(list.files(dir, pattern = "gower_c.tif$",       full.names = TRUE))
  maxent_c      <- stack(list.files(dir, pattern = "maxent_c.tif$",      full.names = TRUE))
  SVM_c         <- stack(list.files(dir, pattern = "SVM_c.tif$",         full.names = TRUE))
  
  bioclim_rcp26 <- stack(list.files(dir, pattern = "bioclim_rcp26.tif$", full.names = TRUE))
  gower_rcp26   <- stack(list.files(dir, pattern = "gower_rcp26.tif$",   full.names = TRUE))
  maxent_rcp26  <- stack(list.files(dir, pattern = "maxent_rcp26.tif$",  full.names = TRUE))
  SVM_rcp26     <- stack(list.files(dir, pattern = "SVM_rcp26.tif$",     full.names = TRUE))
  
  bioclim_rcp45 <- stack(list.files(dir, pattern = "bioclim_rcp45.tif$", full.names = TRUE))
  gower_rcp45   <- stack(list.files(dir, pattern = "gower_rcp45.tif$",   full.names = TRUE))
  maxent_rcp45  <- stack(list.files(dir, pattern = "maxent_rcp45.tif$",  full.names = TRUE))
  SVM_rcp45     <- stack(list.files(dir, pattern = "SVM_rcp45.tif$",     full.names = TRUE))
  
  bioclim_rcp60 <- stack(list.files(dir, pattern = "bioclim_rcp60.tif$", full.names = TRUE))
  gower_rcp60   <- stack(list.files(dir, pattern = "gower_rcp60.tif$",   full.names = TRUE))
  maxent_rcp60  <- stack(list.files(dir, pattern = "maxent_rcp60.tif$",  full.names = TRUE))
  SVM_rcp60     <- stack(list.files(dir, pattern = "SVM_rcp60.tif$",     full.names = TRUE))
  
  bioclim_rcp85 <- stack(list.files(dir, pattern = "bioclim_rcp85.tif$", full.names = TRUE))
  gower_rcp85   <- stack(list.files(dir, pattern = "gower_rcp85.tif$",   full.names = TRUE))
  maxent_rcp85  <- stack(list.files(dir, pattern = "maxent_rcp85.tif$",  full.names = TRUE))
  SVM_rcp85     <- stack(list.files(dir, pattern = "SVM_rcp85.tif$",    full.names = TRUE))
  
  ### Partial outputs                 -----
  ###......................................
  
  bioclim_Pout_c     <- na.omit(values(bioclim_c))
  bioclim_Pout_rcp26 <- na.omit(values(bioclim_rcp26))
  bioclim_Pout_rcp45 <- na.omit(values(bioclim_rcp45))
  bioclim_Pout_rcp60 <- na.omit(values(bioclim_rcp60))
  bioclim_Pout_rcp85 <- na.omit(values(bioclim_rcp85))
  
  gower_Pout_c       <- na.omit(values(gower_c))
  gower_Pout_rcp26   <- na.omit(values(gower_rcp26))
  gower_Pout_rcp45   <- na.omit(values(gower_rcp45))
  gower_Pout_rcp60   <- na.omit(values(gower_rcp60))
  gower_Pout_rcp85   <- na.omit(values(gower_rcp85))
  
  maxent_Pout_c      <- na.omit(values(maxent_c))
  maxent_Pout_rcp26  <- na.omit(values(maxent_rcp26))
  maxent_Pout_rcp45  <- na.omit(values(maxent_rcp45))
  maxent_Pout_rcp60  <- na.omit(values(maxent_rcp60))
  maxent_Pout_rcp85  <- na.omit(values(maxent_rcp85))
  
  SVM_Pout_c         <- na.omit(values(SVM_c))
  SVM_Pout_rcp45     <- na.omit(values(SVM_rcp45))
  SVM_Pout_rcp26     <- na.omit(values(SVM_rcp26))
  SVM_Pout_rcp60     <- na.omit(values(SVM_rcp60))
  SVM_Pout_rcp85     <- na.omit(values(SVM_rcp85))
  
  ### Saving the predictive methods   -----
  ###......................................
  
  output_current <- cbind(Bioclim = bioclim_Pout_c,
                          Gower   = gower_Pout_c,
                          Maxent  = maxent_Pout_c,
                          SVM     = SVM_Pout_c )
  
  output_rcp26   <- cbind(Bioclim = bioclim_Pout_rcp26,
                          Gower   = gower_Pout_rcp26,
                          Maxent  = maxent_Pout_rcp26,
                          SVM     = SVM_Pout_rcp26 )
  
  output_rcp45   <- cbind(Bioclim = bioclim_Pout_rcp45,
                          Gower   = gower_Pout_rcp45,
                          Maxent  = maxent_Pout_rcp45,
                          SVM     = SVM_Pout_rcp45 )
  
  output_rcp60   <- cbind(Bioclim = bioclim_Pout_rcp60,
                          Gower   = gower_Pout_rcp60,
                          Maxent  = maxent_Pout_rcp60,
                          SVM     = SVM_Pout_rcp60 )
  
  output_rcp85   <- cbind(Bioclim = bioclim_Pout_rcp85,
                          Gower   = gower_Pout_rcp85,
                          Maxent  = maxent_Pout_rcp85,
                          SVM     = SVM_Pout_rcp85 )
  
  
  
  ### Standardizing predictions       ----
  ###.....................................
  id_col_fut <- rep(1:ncol(output_current), each = 3)
  id_time    <- c(rep("c", nrow(output_current)), rep(c("rcp26", "rcp45", "rcp60", "rcp85"),
                                                      each = nrow(output_current) * length(AOGCMs)))
  
  pad_c <- pad_rcp26 <- pad_rcp45 <- pad_rcp60 <- pad_rcp85 <- NULL
  for(p in 1:ncol(output_current))
  {
    suit     <- cbind(output_current[, p],
                      output_rcp26[, which(id_col_fut == p)],
                      output_rcp45[, which(id_col_fut == p)],
                      output_rcp60[, which(id_col_fut == p)],
                      output_rcp85[, which(id_col_fut == p)])
    suit     <- as.numeric(suit)
    suit_pad <- decostand(x = suit, method = "range") # requires vegan
    
    pad_c     <- cbind(pad_c, suit_pad[which(id_time == "c"), 1])
    
    pad_rcp26 <- cbind(pad_rcp26, matrix(data = suit_pad[which(id_time == "rcp26"), 1],
                                         nrow = nrow(output_current), ncol = length(AOGCMs)))
    
    pad_rcp45 <- cbind(pad_rcp45, matrix(data = suit_pad[which(id_time == "rcp45"), 1],
                                         nrow = nrow(output_current), ncol = length(AOGCMs)))
    
    pad_rcp60 <- cbind(pad_rcp60, matrix(data = suit_pad[which(id_time == "rcp60"), 1],
                                         nrow = nrow(output_current), ncol = length(AOGCMs)))
    
    pad_rcp85 <- cbind(pad_rcp85, matrix(data = suit_pad[which(id_time == "rcp85"), 1],
                                         nrow = nrow(output_current), ncol = length(AOGCMs)))
  }
  
  
  ### Calculating Ensembles           ----
  ###.....................................
  dStat_c <- c(bioclim_d, gower_d, maxent_d, SVM_d)
  dStat_fut <- rep(dStat_c, each= length(AOGCMs))
  
  ensemble_c     <- apply(output_current, 1, function(x) sum(x * dStat_c)   / sum(dStat_c))
  ensemble_rcp26 <- apply(output_rcp26,   1, function(x) sum(x * dStat_fut) / sum(dStat_fut))
  ensemble_rcp45 <- apply(output_rcp45,   1, function(x) sum(x * dStat_fut) / sum(dStat_fut))
  ensemble_rcp60 <- apply(output_rcp60,   1, function(x) sum(x * dStat_fut) / sum(dStat_fut))
  ensemble_rcp85 <- apply(output_rcp85,   1, function(x) sum(x * dStat_fut) / sum(dStat_fut))
  
  ### Inserting coords to outputs
  coords <- na.omit(cbind(xyFromCell(current, 1:ncell(current)), values(current)))[,1:2]
  
  output_current <- cbind(coords, output_current)
  output_rcp26   <- cbind(coords, output_rcp26)
  output_rcp45   <- cbind(coords, output_rcp45)
  output_rcp60   <- cbind(coords, output_rcp60)
  output_rcp85   <- cbind(coords, output_rcp85)
  
  ### Saving Ensembles
  FULLensemble       <- data.frame(coords,
                                   Ensemble_present = ensemble_c,
                                   Ensemble_rcp26   = ensemble_rcp26,
                                   Ensemble_rcp45   = ensemble_rcp45,
                                   Ensemble_rcp60   = ensemble_rcp60,
                                   Ensemble_rcp85   = ensemble_rcp85)
  save(FULLensemble,   file = paste0(Pout, i, "_ENSEMBLE.Rdata"))
  
  ### Saving Evaluation data          ----
  ##.....................................
  
  models_e <- data.frame(bioclim = bioclim_e, gower = gower_e, maxent = maxent_e, SVM = SVM_e)
  models_t <- data.frame(bioclim = bioclim_t, gower = gower_t, maxent = maxent_t, SVM = SVM_t)
  models_d <- data.frame(bioclim = bioclim_d, gower = gower_d, maxent = maxent_d, SVM = SVM_d)
  
  
  ### Returning Function data         ----
  ###.....................................
  
  return(list(c("output_current"  = output_current, 
                "output_rcp26"    = output_rcp26,
                "output_rcp45"    = output_rcp45,
                "output_rcp60"    = output_rcp60,
                "output_rcp85"    = output_rcp85,
                "Ensemble"        = FULLensemble)))
}
  