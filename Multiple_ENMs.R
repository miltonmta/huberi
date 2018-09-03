multiple_ENMs <- function(occurrence, 
                          background, 
                          biovar_current,
                          biovar_rcp26,
                          biovar_rcp45,
                          biovar_rcp60,
                          biovar_rcp85,
                          newvar_current,
                          newvar_rcp26,
                          newvar_rcp45,
                          newvar_rcp60,
                          newvar_rcp85,
                          trainning,
                          testing,
                          AOGCMs,
                          Pout,
                          cross_validation)
{
  ### Loading data                    ----

  ###.............. Reading the selected climatic variables
  
  if (is.numeric(newvar_current)){
    current <- stack(list.files(biovar_current,  pattern = ".grd$", full.names = TRUE))
  }else{
    current <- addLayer(stack(list.files(biovar_current,  pattern = ".grd$", full.names = TRUE)), 
                        stack(list.files(newvar_current,  pattern = ".grd$", full.names = TRUE)))
  }
  
  ###.............. Reading ocurrence and background
  occur <- read.table(occurrence, sep = ";", h = T)
  back  <- read.table(background, sep = ";", h = T)
  
  ###.............. Creating objects for saving the results
  
  bioclim_c <- gower_c  <- maxent_c <- SVM_c <- stack()
  
  bioclim_e <- gower_e <- maxent_e <- SVM_e <- NULL # True presence rate (TPR)
  bioclim_t <- gower_t <- maxent_t <- SVM_t <- NULL # The highest threshold at which there is no omission
  bioclim_d <- gower_d <- maxent_d <- SVM_d <- NULL # Area predicted as presence
  

  ### Cross validation                ----
  n_cells <- nrow(na.omit(values(current)))
  for (i in 1:cross_validation)
  {
    ### OPEN "i" ----
    
    ###.............. Loading trainning-testing subsets
    
    if (is.character(trainning)){
     trainning <- read.table(paste0(trainning, i, ".txt"), sep = ";")
     testing   <- read.table(paste0(testing,   i, ".txt"), sep = ";")
    }else{
      sample_occur <- sample(1:nrow(occur), round(0.75 * nrow(occur), 0))
      trainning <- prepareData(x = current, 
                               p = occur[sample_occur,  1:2], 
                               b = back [sample_occur,  1:2]) 
      testing   <- prepareData(x = current, 
                               p = occur[-sample_occur, 1:2], 
                               b = back [-sample_occur, 1:2])
    }
    
    # ***************************************************************************************
    ### Bioclim -----------------------------------------------------------------------------
    
    ## Adjusting models
    bioclim_model <- bioclim(trainning[trainning[, "pb"] == 1, -1])
    
    ## Predicting
    bioclim_c <- stack(bioclim_c, predict(object = bioclim_model, x = current))
    writeRaster(bioclim_c, paste0(Pout, "_bioclim_c.tif"), format = "GTiff", overwrite = TRUE)
    
    ## Evaluating models
    thr <- quantile(extract(bioclim_c[[i]], occur[,1:2]), 0.05) 
    TPR <- evaluate(p = testing[testing[, "pb"] == 1, -1], a = testing[testing[, "pb"] == 0, -1], 
                    model = bioclim_model, tr = thr)@TPR
    Pi <- sum(values(bioclim_c[[i]] >= thr), na.rm = T) / n_cells # predicted area proportion.
    
    bioclim_e <- c(bioclim_e, TPR)
    bioclim_t <- c(bioclim_t, thr) 
    bioclim_d <- c(bioclim_d, TPR * (1 - Pi))
    rm(thr, TPR, Pi)
    
    ### Gower -------------------------------------------------------------------------------

    # ## Adjusting models
    # gower_model <- domain(trainning[trainning[,"pb"] == 1, -1])
    # 
    # ## Predicting
    # gower_c <- stack(gower_c, predict(object = gower_model, x = current))
    # writeRaster(gower_c, paste0(Pout, "_gower_c.tif"), format = "GTiff", overwrite = TRUE)
    # 
    # ## Evaluating models
    # thr <- quantile(extract(gower_c[[i]], occur[,1:2]), 0.05)
    # TPR <- evaluate(p = testing[testing[, "pb"] == 1, -1], a = testing[testing[, "pb"] == 0, -1],
    #                 model = gower_model, tr = thr)@TPR
    # Pi <- sum(values(gower_c[[i]] >= thr), na.rm = T) / n_cells
    # 
    # gower_e <- c(gower_e, TPR)
    # gower_t <- c(gower_t, thr)
    # gower_d <- c(gower_d, TPR * (1 - Pi))
    # rm(thr, TPR, Pi)
    
    ### Maxent -------------------------------------------------------------------------------
    
    ## Adjusting models
    Sys.setenv(NOAWT = TRUE)
    maxent_model <- maxent(x = trainning[, -1], p = trainning[, 1])
    
    ## Predicting
    maxent_c <- stack(maxent_c, predict(object = maxent_model, x = current))
    writeRaster(maxent_c, paste0(Pout, "_maxent_c.tif"), format = "GTiff", overwrite = TRUE)    
    
    ## Evaluating models
    thr <- quantile(extract(maxent_c[[i]], occur[,1:2]), 0.05) 
    TPR <- evaluate(p = testing[testing[, "pb"] == 1, -1], a = testing[testing[, "pb"] == 0, -1], 
                    model = maxent_model, tr = thr)@TPR
    Pi <- sum(values(maxent_c[[i]] >= thr), na.rm = T) / n_cells
    
    maxent_e <- c(maxent_e, TPR) 
    maxent_t <- c(maxent_t, thr)
    maxent_d <- c(maxent_d, TPR * (1 - Pi))
    rm(thr, TPR, Pi)
    
    ### SVM ----------------------------------------------------------------------------------
    
    ## Adjusting models
    SVM_model <- ksvm(pb ~ ., data = trainning)
    
    ## Predicting
    SVM_c <- stack(SVM_c, predict(model = SVM_model, object = current))
    writeRaster(SVM_c, paste0(Pout, "_SVM_c.tif"), format = "GTiff", overwrite = TRUE)    
    
    ## Evaluating models
    thr <- quantile(extract(SVM_c[[i]], occur[,1:2]), 0.05) 
    TPR <- evaluate(p = testing[testing[, "pb"] == 1, -1], a = testing[testing[, "pb"] == 0, -1], 
                    model = SVM_model, tr = thr)@TPR
    Pi <- sum(values(SVM_c[[i]] >= thr), na.rm = T) / n_cells
    
    SVM_e <- c(SVM_e, TPR)
    SVM_t <- c(SVM_t, thr)
    SVM_d <- c(SVM_d, TPR * (1 - Pi))
    rm(thr, TPR, Pi)
    
    # ***************************************************************************************
 
    ###.............. Making predictions for the RCPs
 
    # Creanting objects for saving results at each loop
    bioclim_rcp26 <- gower_rcp26 <- maxent_rcp26 <- SVM_rcp26 <- stack() 
    bioclim_rcp45 <- gower_rcp45 <- maxent_rcp45 <- SVM_rcp45 <- stack()
    bioclim_rcp60 <- gower_rcp60 <- maxent_rcp60 <- SVM_rcp60 <- stack()
    bioclim_rcp85 <- gower_rcp85 <- maxent_rcp85 <- SVM_rcp85 <- stack()
    
    # mdid <- paste0('.',1:3,'.grd$')
    # #### OPEN "j" ----
    # for (j in AOGCMs)
    # {
    #  
    #   ### Reading the selected AOGCMs climatic models
    #   if(is.numeric(newvar_rcp26)){
    #     # ..............
    #     mdls <- lapply(mdid, function(x){stack(list.files(biovar_rcp26, 
    #                                                       pattern = x, full.names = TRUE))})
    #     rcp26  <- mdls[[j]] 
    #     names(rcp26) <- names(current)
    #     rm(mdls)
    #     
    #     # ..............
    #     mdls <- lapply(mdid, function(x){stack(list.files(biovar_rcp45, 
    #                                                       pattern = x, full.names = TRUE))})
    #     rcp45  <- mdls[[j]]
    #     names(rcp45) <- names(current)
    #     rm(mdls)
    #     
    #     # ..............
    #     mdls  <- lapply(mdid, function(x){stack(list.files(biovar_rcp60, 
    #                                                        pattern = x, full.names = TRUE))})
    #     rcp60  <- mdls[[j]]
    #     names(rcp60) <- names(current)
    #     rm(mdls)
    #     
    #     # ..............
    #     mdls    <- lapply(mdid, function(x){stack(list.files(biovar_rcp85, 
    #                                                        pattern = x, full.names = TRUE))})
    #     rcp85  <- mdls[[j]]
    #     names(rcp85) <- names(current)
    #     rm(mdls)
    #     
    #   }else{
    #     # ..............
    #     mdls_bio <- lapply(mdid, function(x){stack(list.files(biovar_rcp26, 
    #                                                           pattern = x, full.names = TRUE))})
    #     mdls_new <- lapply(mdid, function(x){stack(list.files(newvar_rcp26, 
    #                                                           pattern = x, full.names = TRUE))})
    #     mdls     <- c(mdls_bio, mdls_new)
    #     rcp26  <- addLayer(mdls[[j]], mdls[[j+3]])
    #     names(rcp26) <- names(current)
    #     rm(mdls_bio, mdls_new, mdls)
    #     
    #     # ..............
    #     mdls_bio <- lapply(mdid, function(x){stack(list.files(biovar_rcp45, 
    #                                                           pattern = x, full.names = TRUE))})
    #     mdls_new <- lapply(mdid, function(x){stack(list.files(newvar_rcp45, 
    #                                                           pattern = x, full.names = TRUE))})
    #     mdls     <- c(mdls_bio, mdls_new)
    #     rcp45  <- addLayer(mdls[[j]], mdls[[j+3]])
    #     names(rcp45) <- names(current)
    #     rm(mdls_bio, mdls_new, mdls)
    #     
    #     # ..............
    #     mdls_bio <- lapply(mdid, function(x){stack(list.files(biovar_rcp60, 
    #                                                           pattern = x, full.names = TRUE))})
    #     mdls_new <- lapply(mdid, function(x){stack(list.files(newvar_rcp60, 
    #                                                           pattern = x, full.names = TRUE))})
    #     mdls     <- c(mdls_bio, mdls_new)
    #     rcp60  <- addLayer(mdls[[j]], mdls[[j+3]])
    #     names(rcp60) <- names(current)
    #     rm(mdls_bio, mdls_new, mdls)
    #     
    #     # ..............
    #     mdls_bio <- lapply(mdid, function(x){stack(list.files(biovar_rcp85, 
    #                                                           pattern = x, full.names = TRUE))})
    #     mdls_new <- lapply(mdid, function(x){stack(list.files(newvar_rcp85, 
    #                                                           pattern = x, full.names = TRUE))})
    #     mdls     <- c(mdls_bio, mdls_new)
    #     rcp85  <- addLayer(mdls[[j]], mdls[[j+3]])
    #     names(rcp85) <- names(current)
    #     rm(mdls_bio, mdls_new, mdls)
    #     
    #   }
    #   
    #   ### Predicting
    #   #..........
    #   bioclim_rcp26 <- stack(bioclim_rcp26, predict(object = bioclim_model, x = rcp26))
    #   writeRaster(bioclim_rcp26, paste0(Pout, i, "_bioclim_rcp26.tif"), format = "GTiff", overwrite = TRUE)
    #   
    #   bioclim_rcp45 <- stack(bioclim_rcp45, predict(object = bioclim_model, x = rcp45))
    #   writeRaster(bioclim_rcp45, paste0(Pout, i, "_bioclim_rcp45.tif"), format = "GTiff", overwrite = TRUE)
    #   
    #   bioclim_rcp60 <- stack(bioclim_rcp60, predict(object = bioclim_model, x = rcp60))
    #   writeRaster(bioclim_rcp60, paste0(Pout, i, "_bioclim_rcp60.tif"), format = "GTiff", overwrite = TRUE)
    #   
    #   bioclim_rcp85 <- stack(bioclim_rcp85, predict(object = bioclim_model, x = rcp85))
    #   writeRaster(bioclim_rcp85, paste0(Pout, i, "_bioclim_rcp85.tif"), format = "GTiff", overwrite = TRUE)
    #   
    #   #..........
    #   gower_rcp26   <- stack(gower_rcp26,   predict(object = gower_model, x = rcp26))
    #   writeRaster(gower_rcp26, paste0(Pout, i, "_gower_rcp26.tif"), format = "GTiff", overwrite = TRUE)
    # 
    #   gower_rcp45   <- stack(gower_rcp45,   predict(object = gower_model, x = rcp45))
    #   writeRaster(gower_rcp45, paste0(Pout, i, "_gower_rcp45.tif"), format = "GTiff", overwrite = TRUE)
    # 
    #   gower_rcp60   <- stack(gower_rcp60,   predict(object = gower_model, x = rcp60))
    #   writeRaster(gower_rcp60, paste0(Pout, i, "_gower_rcp60.tif"), format = "GTiff", overwrite = TRUE)
    # 
    #   gower_rcp85   <- stack(gower_rcp85,   predict(object = gower_model, x = rcp85))
    #   writeRaster(gower_rcp85, paste0(Pout, i, "_gower_rcp85.tif"), format = "GTiff", overwrite = TRUE)
    #   
    #   #..........
    #   maxent_rcp26  <- stack(maxent_rcp26,  predict(object = maxent_model, x = rcp26))
    #   writeRaster(maxent_rcp26, paste0(Pout, i, "_maxent_rcp26.tif"), format = "GTiff", overwrite = TRUE)
    #   
    #   maxent_rcp45  <- stack(maxent_rcp45,  predict(object = maxent_model, x = rcp45))
    #   writeRaster(maxent_rcp45, paste0(Pout, i, "_maxent_rcp45.tif"), format = "GTiff", overwrite = TRUE)
    #   
    #   maxent_rcp60  <- stack(maxent_rcp60,  predict(object = maxent_model, x = rcp60))
    #   writeRaster(maxent_rcp60, paste0(Pout, i, "_maxent_rcp60.tif"), format = "GTiff", overwrite = TRUE)
    #   
    #   maxent_rcp85  <- stack(maxent_rcp85,  predict(object = maxent_model, x = rcp85))
    #   writeRaster(maxent_rcp85, paste0(Pout, i, "_maxent_rcp85.tif"), format = "GTiff", overwrite = TRUE)
    #   
    #   #...........
    #   SVM_rcp26     <- stack(SVM_rcp26,     predict(model = SVM_model, object = rcp26))
    #   writeRaster(SVM_rcp26, paste0(Pout, i, "_SVM_rcp26.tif"), format = "GTiff", overwrite = TRUE)
    #   
    #   SVM_rcp45     <- stack(SVM_rcp45,     predict(model = SVM_model, object = rcp45))
    #   writeRaster(SVM_rcp45, paste0(Pout, i, "_SVM_rcp45.tif"), format = "GTiff", overwrite = TRUE)
    #   
    #   SVM_rcp60     <- stack(SVM_rcp60,     predict(model = SVM_model, object = rcp60))
    #   writeRaster(SVM_rcp60, paste0(Pout, i, "_SVM_rcp60.tif"), format = "GTiff", overwrite = TRUE)
    #   
    #   SVM_rcp85     <- stack(SVM_rcp85,     predict(model = SVM_model, object = rcp85))
    #   writeRaster(SVM_rcp85, paste0(Pout, i, "_SVM_rcp85.tif"), format = "GTiff", overwrite = TRUE)
    #   
    #   # ...........
    #   rm(rcp26, rcp45, rcp60, rcp85)
    #   # CLOSE "i" ----
    #   # AOGCMs
    # }
    # rm(bioclim_rcp26, bioclim_rcp45, bioclim_rcp60, bioclim_rcp85,
    #    #gower_rcp26,   gower_rcp45,   gower_rcp60,   gower_rcp85,
    #    maxent_rcp26,  maxent_rcp45,  maxent_rcp60,  maxent_rcp85,
    #    SVM_rcp26,     SVM_rcp45,     SVM_rcp60,     SVM_rcp85)
    # # CLOSE "j"   ---- 
    # # cross-validation
  }
  
  ### Saving Evaluation data          ----
  ##.....................................
  
  models_e <- data.frame(bioclim = bioclim_e,  maxent = maxent_e, SVM = SVM_e)
  models_t <- data.frame(bioclim = bioclim_t,  maxent = maxent_t, SVM = SVM_t)
  models_d <- data.frame(bioclim = bioclim_d,  maxent = maxent_d, SVM = SVM_d)
  
  
  ### Returning Function data         ----
  ###.....................................
  
  return(list(c("TPR_c"           = models_e, 
                "Threshold_c"     = models_t, 
                "Pred_area_c"     = models_d)))
  
} # CLOSE "Multiple_ENMs"