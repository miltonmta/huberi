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
  
  e <- extent(-109.4583, -29.29167, -56, 14)
  ###.............. Reading the selected climatic variables
  
  if (is.numeric(newvar_current)){
    current <- crop(stack(list.files(biovar_current,  pattern = ".grd$", full.names = TRUE)), e)
  }else{
    current <- addLayer(crop(stack(list.files(biovar_current,  pattern = ".grd$", full.names = TRUE)), e), 
                             stack(list.files(newvar_current,  pattern = "_c",    full.names = TRUE)))
  }
  
  ###.............. Species data
  occur <- read.table(occurrence, sep = ";", h = T)
  back  <- read.table(background, sep = ";", h = T)
  coord.sp <- occur[,c(1,2)]
  
  ###.............. Creating objects to save the results
  
  bioclim_c <- gower_c <- maxent_c <- SVM_c <- Enfa_c <- stack()
  
  bioclim_e <- gower_e <- maxent_e <- SVM_e <- Enfa_e <- NULL # True presence rate (TPR)
  bioclim_t <- gower_t <- maxent_t <- SVM_t <- Enfa_t <- NULL # The highest threshold at which there is no omission
  bioclim_d <- gower_d <- maxent_d <- SVM_d <- Enfa_d <- NULL # Area predicted as presence
  bioclim_t_sens <- gower_t_sens <- maxent_t_sens <- SVM_t_sens <- Enfa_t_sens <- NULL
  bioclim_t_auc <- gower_t_auc <- maxent_t_auc <- SVM_t_auc <- Enfa_t_auc <- NULL
  
  ### Cross validation                ----
  n_cells <- nrow(na.omit(values(current)))
  for (i in 1:cross_validation)
  {
    ### OPEN "i" ----
    
    ###.............. Loading trainning-testing subsets
    
    if (is.character(trainning)){ # For variation control, the PRESENT was modeled usig a fix set of subsets.
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
    # D Statistics and thr no-omission
    thr <- quantile(raster::extract(bioclim_c[[i]], occur[,1:2]), 0.05) 
    TPR <- evaluate(p = testing[testing[, "pb"] == 1, -1], a = testing[testing[, "pb"] == 0, -1], 
                    model = bioclim_model, tr = thr)@TPR
    Pi <- sum(values(bioclim_c[[i]] >= thr), na.rm = T) / n_cells # predicted area proportion.
    
    bioclim_e <- c(bioclim_e, TPR)
    bioclim_t <- c(bioclim_t, thr) 
    bioclim_d <- c(bioclim_d, TPR * (1 - Pi))
    rm(thr, TPR, Pi)
    
    # AUC and thr spec_sens
    bioclim_eval <- evaluate(p = teste[teste[,"pb"] == 1, -1], a = teste[teste[, "pb"] == 0, -1], model = bioclim_model)
    bioclim_t_sens <- c(bioclim_t_sens, threshold(bioclim.eval, "spec_sens"))
    bioclim_auc <- c(bioclim_auc, bioclim_eval@auc)
    
    ### Gower -------------------------------------------------------------------------------
    
    ## Adjusting models
    gower_model <- domain(trainning[trainning[,"pb"] == 1, -1]) 
    
    ## Predicting
    gower_c <- stack(gower_c, predict(object = gower_model, x = current))
    
    ## Evaluating models
    # D Statistics and thr no-omission
    thr <- quantile(raster::extract(gower_c[[i]], occur[,1:2]), 0.05) 
    TPR <- evaluate(p = testing[testing[, "pb"] == 1, -1], a = testing[testing[, "pb"] == 0, -1], 
                    model = gower_model, tr = thr)@TPR
    Pi <- sum(values(gower_c[[i]] >= thr), na.rm = T) / n_cells
    
    gower_e <- c(gower_e, TPR)
    gower_t <- c(gower_t, thr)
    gower_d <- c(gower_d, TPR * (1 - Pi))
    rm(thr, TPR, Pi)
    
    # AUC and thr spec_sens
    gower_eval <- evaluate(p = teste[teste[,"pb"] == 1, -1], a = teste[teste[, "pb"] == 0, -1], model = gower_model)
    gower_t_sens <- c(gower_t_sens, threshold(gower.eval, "spec_sens"))
    gower_auc <- c(gower_auc, gower_eval@auc)
    
    ### Maxent -------------------------------------------------------------------------------
    
    ## Adjusting models
    Sys.setenv(NOAWT = TRUE)
    maxent_model <- maxent(x = trainning[, -1], p = trainning[, 1])
    
    ## Predicting
    maxent_c <- stack(maxent_c, predict(object = maxent_model, x = current))
    writeRaster(maxent_c, paste0(Pout, "_maxent_c.tif"), format = "GTiff", overwrite = TRUE)    
    
    ## Evaluating models
    # D Statistics and thr no-omission
    thr <- quantile(raster::extract(maxent_c[[i]], occur[,1:2]), 0.05) 
    TPR <- evaluate(p = testing[testing[, "pb"] == 1, -1], a = testing[testing[, "pb"] == 0, -1], 
                    model = maxent_model, tr = thr)@TPR
    Pi <- sum(values(maxent_c[[i]] >= thr), na.rm = T) / n_cells
    
    maxent_e <- c(maxent_e, TPR) 
    maxent_t <- c(maxent_t, thr)
    maxent_d <- c(maxent_d, TPR * (1 - Pi))
    rm(thr, TPR, Pi)
    
    # AUC and thr spec_sens
    maxent_eval <- evaluate(p = teste[teste[,"pb"] == 1, -1], a = teste[teste[, "pb"] == 0, -1], model = maxent_model)
    maxent_t_sens <- c(maxent_t_sens, threshold(maxent.eval, "spec_sens"))
    maxent_auc <- c(maxent_auc, maxent_eval@auc)
    
    ### SVM ----------------------------------------------------------------------------------
    
    ## Adjusting models
    SVM_model <- ksvm(pb ~ ., data = trainning)
    
    ## Predicting
    SVM_c <- stack(SVM_c, predict(model = SVM_model, object = current))
    writeRaster(SVM_c, paste0(Pout, "_SVM_c.tif"), format = "GTiff", overwrite = TRUE)    
    
    ## Evaluating models
    # D Statistics and thr no-omission
    thr <- quantile(raster::extract(SVM_c[[i]], occur[,1:2]), 0.05) 
    TPR <- evaluate(p = testing[testing[, "pb"] == 1, -1], a = testing[testing[, "pb"] == 0, -1], 
                    model = SVM_model, tr = thr)@TPR
    Pi <- sum(values(SVM_c[[i]] >= thr), na.rm = T) / n_cells
    
    SVM_e <- c(SVM_e, TPR)
    SVM_t <- c(SVM_t, thr)
    SVM_d <- c(SVM_d, TPR * (1 - Pi))
    rm(thr, TPR, Pi)
    
    # AUC and thr spec_sens
    SVM_eval <- evaluate(p = teste[teste[,"pb"] == 1, -1], a = teste[teste[, "pb"] == 0, -1], model = SVM_model)
    SVM_t_sens <- c(SVM_t_sens, threshold(SVM.eval, "spec_sens"))
    SVM_auc <- c(SVM_auc, SVM_eval@auc)
    
    ### Enfa ----------------------------------------------------------------------------------
    current.values <- values(current) 
    current.spdf <- na.omit(data.frame(xyFromCell(current, 1:ncell(current)), current.values))
    gridded(current.spdf) <- ~ x + y
    
    pr.cell <- extract(current, coord.sp, cellnumber=T) #"coord.sp" deve ter as coordenadas da espécie relacionadas ao subconjunto 'treino'
    pr <- data.frame(pr= rep(0, ncell(current)), current.values)
    pr[pr.cell[,"cells"],] <- 1
    pr <- na.omit(pr)
    pr <- pr[,1]	#o objeto "pr" contêm um vetor de 0/1s indicando as células onde a sp está presente (valores 1s)
    
    ## Standardizing variables
    media.current <- apply(slot(current.spdf, "data"), 2, mean)
    sd.current <- apply(slot(current.spdf, "data"), 2, sd) 
    current.scale<- sweep(slot(current.spdf, "data"),2, media.current) 
    current.scale<- as.matrix(current.scale) %*% diag(1/sd.current)
    colnames(current.scale) <- names(current.spdf)
    
    ## Adjusting models
    # library(adehabitatHR)
    Enfa_model <- madifa(dudi.pca(current.scale, center=F, scale=F, scannf=F), pr, scannf=F)
    
    ## Predicting
    # source("./predict_enfa.R")
    Enfa_c <- Predict.enfa(object.enfa= Enfa_model, baseline.climate= current.spdf, new.climate= current.spdf)
    rm(current.values, pr.cell, pr)
    
    ## Evaluating models
    # D Statistics and thr no-omission
    thr <- quantile(raster::extract(Enfa_c[[i]], occur[,1:2]), 0.05) 
    TPR <- evaluate(p = testing[testing[, "pb"] == 1, -1], a = testing[testing[, "pb"] == 0, -1], 
                    model = Enfa_model, tr = thr)@TPR
    Pi <- sum(values(Enfa_c[[i]] >= thr), na.rm = T) / n_cells
    
    Enfa_e <- c(Enfa_e, TPR)
    Enfa_t <- c(Enfa_t, thr)
    Enfa_d <- c(Enfa_d, TPR * (1 - Pi))
    rm(thr, TPR, Pi)
    
    # AUC and thr spec_sens
    Enfa_eval <- evaluate(p = teste[teste[,"pb"] == 1, -1], a = teste[teste[, "pb"] == 0, -1], model = Enfa_model)
    Enfa_t_sens <- c(Enfa_t_sens, threshold(Enfa.eval, "spec_sens"))
    Enfa_auc <- c(Enfa_auc, Enfa_eval@auc)
    
    # ***************************************************************************************
 
    ###.............. Making predictions for the RCPs
 
    # Creanting objects for saving results at each loop
    bioclim_rcp26 <- gower_rcp26 <- maxent_rcp26 <- SVM_rcp26 <- Enfa_rcp26 <-  stack() 
    bioclim_rcp45 <- gower_rcp45 <- maxent_rcp45 <- SVM_rcp45 <- Enfa_rcp45 <-  stack()
    bioclim_rcp60 <- gower_rcp60 <- maxent_rcp60 <- SVM_rcp60 <- Enfa_rcp60 <-  stack()
    bioclim_rcp85 <- gower_rcp85 <- maxent_rcp85 <- SVM_rcp85 <- Enfa_rcp85 <-  stack()
    
    mdid   <- paste0('.',1:3,'.grd$')
    mdid2  <- paste0('.',1:3,'.tif$')
    #### OPEN "j" ----
    
    for (j in AOGCMs)
    {

      ### Reading the selected AOGCMs climatic models
      if(is.numeric(newvar_rcp26)){
        # ..............
        mdls <- lapply(mdid, function(x){stack(list.files(biovar_rcp26,
                                                          pattern = x, full.names = TRUE))})
        rcp26  <- crop(mdls[[j]], e)
        names(rcp26) <- names(current)
        rm(mdls)
        gc()

        # ..............
        mdls <- lapply(mdid, function(x){stack(list.files(biovar_rcp45,
                                                          pattern = x, full.names = TRUE))})
        rcp45  <- crop(mdls[[j]], e)
        names(rcp45) <- names(current)
        rm(mdls)
        gc()

        # ..............
        mdls  <- lapply(mdid, function(x){stack(list.files(biovar_rcp60,
                                                           pattern = x, full.names = TRUE))})
        rcp60  <- crop(mdls[[j]], e)
        names(rcp60) <- names(current)
        rm(mdls)
        gc()

        # ..............
        mdls    <- lapply(mdid, function(x){stack(list.files(biovar_rcp85,
                                                           pattern = x, full.names = TRUE))})
        rcp85  <- crop(mdls[[j]], e)
        names(rcp85) <- names(current)
        rm(mdls)
        gc()

      }else{
        # ..............
        mdls_bio <- lapply(mdid,  function(x){stack(list.files(biovar_rcp26,
                                                               pattern = x, full.names = TRUE))})
        mdls_new <- lapply(mdid2, function(x){stack(list.files(newvar_rcp26, 
                                                               pattern = x, full.names = TRUE))})
                                                              
        rcp26    <- addLayer(crop(mdls_bio[[j]], e), mdls_new[[j]])
        names(rcp26) <- names(current)
        rm(mdls_bio, mdls_new)
        gc()

        # ..............
        mdls_bio <- lapply(mdid,  function(x){stack(list.files(biovar_rcp45,
                                                               pattern = x, full.names = TRUE))})
        mdls_new <- lapply(mdid2, function(x){stack(list.files(newvar_rcp45, 
                                                               pattern = x, full.names = TRUE))})
                                                            
        rcp45    <- addLayer(crop(mdls_bio[[j]], e), mdls_new[[j]])
        names(rcp45) <- names(current)
        rm(mdls_bio, mdls_new)
        gc()

        # ..............
        mdls_bio <- lapply(mdid, function(x){stack(list.files(biovar_rcp60,
                                                              pattern = x, full.names = TRUE))})
        mdls_new <- lapply(mdid2, function(x){stack(list.files(newvar_rcp60, 
                                                              pattern = x, full.names = TRUE))})
        
        rcp60  <- addLayer(crop(mdls_bio[[j]], e), mdls_new[[j]])
        names(rcp60) <- names(current)
        rm(mdls_bio, mdls_new)
        gc()

        # ..............
        mdls_bio <- lapply(mdid,  function(x){stack(list.files(biovar_rcp85,
                                                               pattern = x, full.names = TRUE))})
        mdls_new <- lapply(mdid2, function(x){stack(list.files(newvar_rcp85, 
                                                               pattern = x, full.names = TRUE))})
        
        rcp85  <- addLayer(crop(mdls_bio[[j]], e), mdls_new[[j]])
        names(rcp85) <- names(current)
        rm(mdls_bio, mdls_new)
        gc()}

      ### Predicting
      #..........
      bioclim_rcp26 <- stack(bioclim_rcp26, predict(object = bioclim_model, x = rcp26))
      writeRaster(bioclim_rcp26, paste0(Pout, i, "_bioclim_rcp26.tif"), format = "GTiff", overwrite = TRUE)

      bioclim_rcp45 <- stack(bioclim_rcp45, predict(object = bioclim_model, x = rcp45))
      writeRaster(bioclim_rcp45, paste0(Pout, i, "_bioclim_rcp45.tif"), format = "GTiff", overwrite = TRUE)

      bioclim_rcp60 <- stack(bioclim_rcp60, predict(object = bioclim_model, x = rcp60))
      writeRaster(bioclim_rcp60, paste0(Pout, i, "_bioclim_rcp60.tif"), format = "GTiff", overwrite = TRUE)

      bioclim_rcp85 <- stack(bioclim_rcp85, predict(object = bioclim_model, x = rcp85))
      writeRaster(bioclim_rcp85, paste0(Pout, i, "_bioclim_rcp85.tif"), format = "GTiff", overwrite = TRUE)

      #..........
      gower_rcp26   <- stack(gower_rcp26,   predict(object = gower_model, x = rcp26))
      writeRaster(gower_rcp26, paste0(Pout, i, "_gower_rcp26.tif"), format = "GTiff", overwrite = TRUE)
      
      gower_rcp45   <- stack(gower_rcp45,   predict(object = gower_model, x = rcp45))
      writeRaster(gower_rcp45, paste0(Pout, i, "_gower_rcp26.tif"), format = "GTiff", overwrite = TRUE)
      
      gower_rcp60   <- stack(gower_rcp60,   predict(object = gower_model, x = rcp60))
      writeRaster(gower_rcp60, paste0(Pout, i, "_gower_rcp26.tif"), format = "GTiff", overwrite = TRUE)
      
      gower_rcp85   <- stack(gower_rcp85,   predict(object = gower_model, x = rcp85))
      writeRaster(gower_rcp85, paste0(Pout, i, "_gower_rcp26.tif"), format = "GTiff", overwrite = TRUE)
      
      #..........
      maxent_rcp26  <- stack(maxent_rcp26,  predict(object = maxent_model, x = rcp26))
      writeRaster(maxent_rcp26, paste0(Pout, i, "_maxent_rcp26.tif"), format = "GTiff", overwrite = TRUE)

      maxent_rcp45  <- stack(maxent_rcp45,  predict(object = maxent_model, x = rcp45))
      writeRaster(maxent_rcp45, paste0(Pout, i, "_maxent_rcp45.tif"), format = "GTiff", overwrite = TRUE)

      maxent_rcp60  <- stack(maxent_rcp60,  predict(object = maxent_model, x = rcp60))
      writeRaster(maxent_rcp60, paste0(Pout, i, "_maxent_rcp60.tif"), format = "GTiff", overwrite = TRUE)

      maxent_rcp85  <- stack(maxent_rcp85,  predict(object = maxent_model, x = rcp85))
      writeRaster(maxent_rcp85, paste0(Pout, i, "_maxent_rcp85.tif"), format = "GTiff", overwrite = TRUE)

      #...........
      SVM_rcp26     <- stack(SVM_rcp26,     predict(model = SVM_model, object = rcp26))
      writeRaster(SVM_rcp26, paste0(Pout, i, "_SVM_rcp26.tif"), format = "GTiff", overwrite = TRUE)

      SVM_rcp45     <- stack(SVM_rcp45,     predict(model = SVM_model, object = rcp45))
      writeRaster(SVM_rcp45, paste0(Pout, i, "_SVM_rcp45.tif"), format = "GTiff", overwrite = TRUE)

      SVM_rcp60     <- stack(SVM_rcp60,     predict(model = SVM_model, object = rcp60))
      writeRaster(SVM_rcp60, paste0(Pout, i, "_SVM_rcp60.tif"), format = "GTiff", overwrite = TRUE)

      SVM_rcp85     <- stack(SVM_rcp85,     predict(model = SVM_model, object = rcp85))
      writeRaster(SVM_rcp85, paste0(Pout, i, "_SVM_rcp85.tif"), format = "GTiff", overwrite = TRUE)

      #...........
      Enfa_rcp26 <- stack(Enfa_rcp26, Predict.enfa(object.enfa = Enfa_model, baseline.climate = current.spdf, new.climate = rcp26))
      writeRaster(Enfa_rcp26, paste0(Pout, i, "_Enfa_rcp26.tif"), format = "GTiff", overwrite = TRUE)
      
      Enfa_rcp45 <- stack(Enfa_rcp45, Predict.enfa(object.enfa = Enfa_model, baseline.climate = current.spdf, new.climate = rcp45))
      writeRaster(Enfa_rcp45, paste0(Pout, i, "_Enfa_rcp45.tif"), format = "GTiff", overwrite = TRUE)
      
      Enfa_rcp60 <- stack(Enfa_rcp60, Predict.enfa(object.enfa = Enfa_model, baseline.climate = current.spdf, new.climate = rcp60))
      writeRaster(Enfa_rcp60, paste0(Pout, i, "_Enfa_rcp60.tif"), format = "GTiff", overwrite = TRUE)
      
      Enfa_rcp85 <- stack(Enfa_rcp85, Predict.enfa(object.enfa = Enfa_model, baseline.climate = current.spdf, new.climate = rcp85))
      writeRaster(Enfa_rcp85, paste0(Pout, i, "_Enfa_rcp85.tif"), format = "GTiff", overwrite = TRUE)
      
     
      # ...........
      rm(rcp26, rcp45, rcp60, rcp85)
      gc()
      # CLOSE "i" ----
      # AOGCMs
    }
    rm(bioclim_rcp26, bioclim_rcp45, bioclim_rcp60, bioclim_rcp85,
       gower_rcp26,   gower_rcp45,   gower_rcp60,   gower_rcp85,
       maxent_rcp26,  maxent_rcp45,  maxent_rcp60,  maxent_rcp85,
       SVM_rcp26,     SVM_rcp45,     SVM_rcp60,     SVM_rcp85,
       Enfa_rcp26,    Enfa_rcp45,    Enfa_rcp60,    Enfa_rcp85)
    gc()
    # CLOSE "j"   ----
    # cross-validation
  }
  
  ### Saving Evaluation data          ----
  ##.....................................
  
  models_e <- data.frame(bioclim = bioclim_e,  gower = gower_e, maxent = maxent_e, SVM = SVM_e)
  models_t <- data.frame(bioclim = bioclim_t,  gower = gower_t, maxent = maxent_t, SVM = SVM_t)
  models_d <- data.frame(bioclim = bioclim_d,  gower = gower_d, maxent = maxent_d, SVM = SVM_d)
  
  models_t_sens <- data.frame(bioclim = bioclim_t_sens,  gower = gower_t_sens, maxent = maxent_t_sens, SVM = SVM_t_sens, ENFA = Enfa_t_sens)
  models_auc <- data.frame(bioclim = bioclim_auc,  gower = gower_auc, maxent = maxent_auc, SVM = SVM_auc, ENFA = Enfa_auc)

  
  return(list("TPR_c"           = models_e, 
              "Threshold_c"     = models_t, 
              "Pred_area_c"     = models_d,
              "AUC"             = models_auc,
              "Threshold_sens"  = models_t_sens))
  
} # CLOSE "Multiple_ENMs"