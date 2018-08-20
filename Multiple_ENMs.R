multiple_ENMs <- function(occurrence, 
                          background, 
                          biovar_current,
                          biovar_rcp26,
                          biovar_rcp45,
                          biovar_rcp60,
                          biovar_rcp85 ,
                          cross_validation)
{
  
  ###.............. Reading the selected bioclimatic variables
  
  ## Reading the selected current model
  current <- stack(list.files(biovar_current,  pattern = ".grd$", full.names = TRUE))
  
  
  ## Reading the selected aogcm models
  mdid <- paste0('.',1:3,'.grd$') # models identification. 
  model_names <- c("CCSM4", "IPSL-CM5A-LR", "MIROC-ESM") # for changing mdid to the proper models names
  
  mdls <- lapply(mdid, function(x){stack(list.files(biovar_rcp26, pattern = x, full.names = TRUE))})
  names(mdls) <- model_names
  rcp26  <- mdls[[i]] # now only the ith model will be extracted from the mdls.
  
  mdls <- lapply(mdid, function(x){stack(list.files(biovar_rcp45, pattern = x, full.names = TRUE))})
  names(mdls) <- model_names
  rcp45  <- mdls[[i]]
  
  mdls  <- lapply(mdid, function(x){stack(list.files(biovar_rcp60, pattern = x, full.names = TRUE))})
  names(mdls) <- model_names
  rcp60  <- mdls[[i]]
  
  mdls  <- lapply(mdid, function(x){stack(list.files(biovar_rcp85, pattern = x, full.names = TRUE))})
  names(mdls) <- model_names
  rcp85  <- mdls[[i]]
  
  ###.............. Reading ocurrence and background
  occur <- read.table(occurrence, sep = ";", h = T)
  back  <- read.table(background, sep = ";", h = T)
  
  ###.............. Creating objects for saving the results
  ## Predictive models
  bioclim_c <- gower_c  <- maxent_c <- SVM_c <- stack()
  bioclim_rcp26 <- gower_rcp26 <- maxent_rcp26 <- SVM_rcp26 <- stack() 
  bioclim_rcp45 <- gower_rcp45 <- maxent_rcp45 <- SVM_rcp45 <- stack()
  bioclim_rcp60 <- gower_rcp60 <- maxent_rcp60 <- SVM_rcp60 <- stack()
  bioclim_rcp85 <- gower_rcp85 <- maxent_rcp85 <- SVM_rcp85 <- stack()
  
  ## Evaluation
  # present
  bioclim_e <- gower_e <- maxent_e <- SVM_e <- NULL # True presence rate (TPR)
  bioclim_t <- gower_t <- maxent_t <- SVM_t <- NULL # The highest threshold at which there is no omission
  bioclim_d <- gower_d <- maxent_d <- SVM_d <- NULL # Area predicted as presence
  
  # rcp26
  bioclim_e_rcp26 <- gower_e_rcp26 <- maxent_e_rcp26 <- SVM_e_rcp26 <- NULL 
  bioclim_t_rcp26 <- gower_t_rcp26 <- maxent_t_rcp26 <- SVM_t_rcp26 <- NULL 
  bioclim_d_rcp26 <- gower_d_rcp26 <- maxent_d_rcp26 <- SVM_d_rcp26 <- NULL 
  
  # rcp45
  bioclim_e_rcp45 <- gower_e_rcp45 <- maxent_e_rcp45 <- SVM_e_rcp45 <- NULL 
  bioclim_t_rcp45 <- gower_t_rcp45 <- maxent_t_rcp45 <- SVM_t_rcp45 <- NULL 
  bioclim_d_rcp45 <- gower_d_rcp45 <- maxent_d_rcp45 <- SVM_d_rcp45 <- NULL 
  
  # rcp60
  bioclim_e_rcp60 <- gower_e_rcp60 <- maxent_e_rcp60 <- SVM_e_rcp60 <- NULL 
  bioclim_t_rcp60 <- gower_t_rcp60 <- maxent_t_rcp60 <- SVM_t_rcp60 <- NULL 
  bioclim_d_rcp60 <- gower_d_rcp60 <- maxent_d_rcp60 <- SVM_d_rcp60 <- NULL 
  
  # rcp85
  bioclim_e_rcp85 <- gower_e_rcp85 <- maxent_e_rcp85 <- SVM_e_rcp85 <- NULL 
  bioclim_t_rcp85 <- gower_t_rcp85 <- maxent_t_rcp85 <- SVM_t_rcp85 <- NULL 
  bioclim_d_rcp85 <- gower_d_rcp85 <- maxent_d_rcp85 <- SVM_d_rcp85 <- NULL 
  
  ## Partial results
  bioclim_Pout_c <- gower_Pout_c <- maxent_Pout_c <- SVM_Pout_c <- NULL
  bioclim_Pout_rcp26 <- gower_Pout_rcp26 <- maxent_Pout_rcp26 <- SVM_Pout_rcp26 <- NULL
  bioclim_Pout_rcp45 <- gower_Pout_rcp45 <- maxent_Pout_rcp45 <- SVM_Pout_rcp45 <- NULL
  bioclim_Pout_rcp60 <- gower_Pout_rcp60 <- maxent_Pout_rcp60 <- SVM_Pout_rcp60 <- NULL
  bioclim_Pout_rcp85 <- gower_Pout_rcp85 <- maxent_Pout_rcp85 <- SVM_Pout_rcp85 <- NULL
  
  ## Predictions results
  output_current <- output_rcp26 <- output_rcp45 <- output_rcp60 <- output_rcp85 <- NULL
  
  
  for (j in 1:cross_validation)
  {
    ### OPEN "j" ----
    
    ###.............. Creating trainning-testing subsets
    sample_occur <- sample(1:nrow(occur), round(0.75 * nrow(occur), 0))
    trainning <- prepareData(x = current, p = occur[sample_occur,  1:2], b = back[sample_occur,  1:2]) 
    testing   <- prepareData(x = current, p = occur[-sample_occur, 1:2], b = back[-sample_occur, 1:2])
    
    
    # ***************************************************************************************
    ### Bioclim -----------------------------------------------------------------------------
    
    ## Adjusting models
    bioclim_model <- bioclim(trainning[trainning[, "pb"] == 1, -1])
    
    ## Predicting
    bioclim_c <- stack(bioclim_c, predict(object = bioclim_model, x = current))
    
    ## Evaluating models
    bioclim_eval <- evaluate(p = testing[testing[, "pb"] == 1, -1], a = testing[testing[, "pb"] == 0, -1], model = bioclim_model)
    
    ## Saving evaluation metrics
    thr <- threshold(bioclim_eval, "no_omission") 
    TPR <- bioclim_eval@TPR[which( bioclim_eval@t == thr)]
    n_cells <- nrow(na.omit(values(current))) 
    pi <- sum(values(bioclim_c >= thr), na.rm = T) / n_cells # predicted area proportion.
    
    bioclim_e <- c(bioclim_e, TPR)
    bioclim_t <- c(bioclim_t, thr) 
    bioclim_d <- TPR * (1 - pi)
    
    
    ### Gower -------------------------------------------------------------------------------
    
    ## Adjusting models
    gower_model <- domain(trainning[trainning[,"pb"] == 1, -1]) 
    
    ## Predicting
    gower_c <- stack(gower_c, predict(object = gower_model, x = current))
    
    ## Evaluating models
    gower_eval <- evaluate(p = testing[testing[, "pb"] == 1, -1], a = testing[testing[, "pb"] == 0, -1], model = gower_model)
    
    ## Saving evaluation metrics
    thr <- threshold(gower_eval, "no_omission")
    TPR <- gower_eval@TPR[which( gower_eval@t == thr)]
    n_cells <- nrow(na.omit(values(current))) 
    pi <- sum(values(gower_c >= thr), na.rm = T) / n_cells
    
    gower_e <- c(gower_e, TPR)
    gower_t <- c(gower_t, thr)
    gower_d <- TPR * (1 - pi)
    
    
    ### Maxent -------------------------------------------------------------------------------
    
    ## Adjusting models
    Sys.setenv(NOAWT = TRUE)
    maxent_model <- maxent(x = trainning[, -1], p = trainning[, 1])
    
    ## Predicting
    maxent_c <- stack(maxent_c, predict(object = maxent_model, x = current))
    
    ## Evaluating models
    maxent_eval <- evaluate(p = testing[testing[, "pb"] == 1, -1], a = testing[testing[, "pb"] == 0, -1], model = maxent_model, tr = 0.5414856)
    
    ## Saving evaluation metrics
    thr <- threshold(maxent_eval, "no_omission")
    TPR <- maxent_eval@TPR[which( maxent_eval@t == thr)]
    n_cells <- nrow(na.omit(values(current))) 
    pi <- sum(values(maxent_c >= thr), na.rm = T) / n_cells
    
    maxent_e <- c(maxent_e, TPR) 
    maxent_t <- c(maxent_t, thr)
    maxent_d <- TPR * (1 - pi)
    
    
    ### SVM ----------------------------------------------------------------------------------
    
    ## Adjusting models
    SVM_model <- ksvm(pb ~ ., data = trainning)
    
    ## Predicting
    SVM_c <- stack(SVM_c, predict(model = SVM_model, object = current)) 
    
    ## Evaluating models
    SVM_eval <- evaluate(p = testing[testing[, "pb"] == 1, -1], a = testing[testing[, "pb"] == 0, -1], model = SVM_model, tr = 0.5566005)
    
    ## Saving evaluation metrics
    thr <- threshold(SVM_eval, "no_omission")
    TPR <- SVM_eval@TPR[which( SVM_eval@t == thr)]
    n_cells <- nrow(na.omit(values(current))) 
    pi <- sum(values(SVM_c >= thr), na.rm = T) / n_cells
    
    SVM_e <- c(SVM_e, TPR)
    SVM_t <- c(SVM_t, thr)
    SVM_d <- TPR * (1 - pi)
    
    
    # ***************************************************************************************
    ###.............. Saving partial outputs for current model 
    
    bioclim_Pout_c <- cbind(bioclim_Pout_c, values(bioclim_c))
    gower_Pout_c   <- cbind(gower_Pout_c,   values(gower_c))
    maxent_Pout_c  <- cbind(maxent_Pout_c,  values(maxent_c))
    SVM_Pout_c     <- cbind(SVM_Pout_c,     values(SVM_c))
    
    
    ###.............. Making predictions for the RCPs
    
    AOGCMs <- c("CCSM4", "IPSL-CM5A-LR", "MIROC-ESM")
    for (i in AOGCMs)
    {
      # OPEN "i" ----
      
      ### Predicting
      bioclim_rcp26 <- stack(bioclim_rcp26, predict(object = bioclim_model, x = rcp26))
      bioclim_rcp45 <- stack(bioclim_rcp45, predict(object = bioclim_model, x = rcp45))
      bioclim_rcp60 <- stack(bioclim_rcp60, predict(object = bioclim_model, x = rcp60))
      bioclim_rcp85 <- stack(bioclim_rcp85, predict(object = bioclim_model, x = rcp85))
      
      gower_rcp26   <- stack(gower_rcp26,   predict(object = gower_model, x = rcp26))
      gower_rcp45   <- stack(gower_rcp45,   predict(object = gower_model, x = rcp45))
      gower_rcp60   <- stack(gower_rcp60,   predict(object = gower_model, x = rcp60))
      gower_rcp85   <- stack(gower_rcp85,   predict(object = gower_model, x = rcp85))
      
      maxent_rcp26  <- stack(maxent_rcp26,  predict(object = maxent_model, x = rcp26))
      maxent_rcp45  <- stack(maxent_rcp45,  predict(object = maxent_model, x = rcp45))
      maxent_rcp60  <- stack(maxent_rcp60,  predict(object = maxent_model, x = rcp60))
      maxent_rcp85  <- stack(maxent_rcp85,  predict(object = maxent_model, x = rcp85))
      
      SVM_rcp26     <- stack(SVM_rcp26,     predict(model = SVM_model, object = rcp26))
      SVM_rcp45     <- stack(SVM_rcp45,     predict(model = SVM_model, object = rcp45))
      SVM_rcp60     <- stack(SVM_rcp60,     predict(model = SVM_model, object = rcp60))
      SVM_rcp85     <- stack(SVM_rcp85,     predict(model = SVM_model, object = rcp85))
      
      ### Salving evaluation metrics
      ###.............. bioclim
      
      thr <- threshold(bioclim_eval, "no_omission")
      TPR <- bioclim_eval@TPR[which( bioclim_eval@t == thr)]
      
      ### rcp26
      n_cells <- nrow(na.omit(values(rcp26))) 
      pi <- sum(values(bioclim_rcp26 >= thr), na.rm = T) / n_cells
      bioclim_e_rcp26 <- c(bioclim_e, TPR)
      bioclim_t_rcp26 <- c(bioclim_t, thr)
      bioclim_d_rcp26 <- TPR * (1 - pi)
      
      ### rcp45
      n_cells <- nrow(na.omit(values(rcp45))) 
      pi <- sum(values(bioclim_rcp45 >= thr), na.rm = T) / n_cells
      bioclim_e_rcp45 <- c(bioclim_e, TPR)
      bioclim_t_rcp45 <- c(bioclim_t, thr)
      bioclim_d_rcp45 <- TPR * (1 - pi)
      
      ### rcp60
      n_cells <- nrow(na.omit(values(rcp60))) 
      pi <- sum(values(bioclim_rcp60 >= thr), na.rm = T) / n_cells
      bioclim_e_rcp60 <- c(bioclim_e, TPR)
      bioclim_t_rcp60 <- c(bioclim_t, thr)
      bioclim_d_rcp60 <- TPR * (1 - pi)
      
      ### rcp85
      n_cells <- nrow(na.omit(values(rcp85))) 
      pi <- sum(values(bioclim_rcp85 >= thr), na.rm = T) / n_cells
      bioclim_e_rcp85 <- c(bioclim_e, TPR)
      bioclim_t_rcp85 <- c(bioclim_t, thr)
      bioclim_d_rcp85 <- TPR * (1 - pi)
      
      ###.............. gower
      
      thr <- threshold(gower_eval, "no_omission")
      TPR <- gower_eval@TPR[which( gower_eval@t == thr)]
      
      ### rcp26
      n_cells <- nrow(na.omit(values(rcp26))) 
      pi <- sum(values(gower_rcp26 >= thr), na.rm = T) / n_cells
      gower_e_rcp26 <- c(gower_e, TPR)
      gower_t_rcp26 <- c(gower_t, thr)
      gower_d_rcp26 <- TPR * (1 - pi)
      
      ### rcp45
      n_cells <- nrow(na.omit(values(rcp45))) 
      pi <- sum(values(gower_rcp45 >= thr), na.rm = T) / n_cells
      gower_e_rcp45 <- c(gower_e, TPR)
      gower_t_rcp45 <- c(gower_t, thr)
      gower_d_rcp45 <- TPR * (1 - pi)
      
      ### rcp60
      n_cells <- nrow(na.omit(values(rcp60))) 
      pi <- sum(values(gower_rcp60 >= thr), na.rm = T) / n_cells
      gower_e_rcp60 <- c(gower_e, TPR)
      gower_t_rcp60 <- c(gower_t, thr)
      gower_d_rcp60 <- TPR * (1 - pi)
      
      ### rcp85
      n_cells <- nrow(na.omit(values(rcp85))) 
      pi <- sum(values(gower_rcp85 >= thr), na.rm = T) / n_cells
      gower_e_rcp85 <- c(gower_e, TPR)
      gower_t_rcp85 <- c(gower_t, thr)
      gower_d_rcp85 <- TPR * (1 - pi)
      
      ###.............. maxent
      
      thr <- threshold(maxent_eval, "no_omission")
      TPR <- maxent_eval@TPR[which( maxent_eval@t == thr)]
      
      ### rcp26
      n_cells <- nrow(na.omit(values(rcp26))) 
      pi <- sum(values(maxent_rcp26 >= thr), na.rm = T) / n_cells
      maxent_e_rcp26 <- c(maxent_e, TPR)
      maxent_t_rcp26 <- c(maxent_t, thr)
      maxent_d_rcp26 <- TPR * (1 - pi)
      
      ### rcp45
      n_cells <- nrow(na.omit(values(rcp45))) 
      pi <- sum(values(maxent_rcp45 >= thr), na.rm = T) / n_cells
      maxent_e_rcp45 <- c(maxent_e, TPR)
      maxent_t_rcp45 <- c(maxent_t, thr)
      maxent_d_rcp45 <- TPR * (1 - pi)
      
      ### rcp60
      n_cells <- nrow(na.omit(values(rcp60))) 
      pi <- sum(values(maxent_rcp60 >= thr), na.rm = T) / n_cells
      maxent_e_rcp60 <- c(maxent_e, TPR)
      maxent_t_rcp60 <- c(maxent_t, thr)
      maxent_d_rcp60 <- TPR * (1 - pi)
      
      ### rcp85
      n_cells <- nrow(na.omit(values(rcp85))) 
      pi <- sum(values(maxent_rcp85 >= thr), na.rm = T) / n_cells
      maxent_e_rcp85 <- c(maxent_e, TPR)
      maxent_t_rcp85 <- c(maxent_t, thr)
      maxent_d_rcp85 <- TPR * (1 - pi)
      
      ###.............. SVM
      
      thr <- threshold(SVM_eval, "no_omission")
      TPR <- SVM_eval@TPR[which( SVM_eval@t == thr)]
      
      ### rcp26
      n_cells <- nrow(na.omit(values(rcp26))) 
      pi <- sum(values(SVM_rcp26 >= thr), na.rm = T) / n_cells
      SVM_e_rcp26 <- c(SVM_e, TPR)
      SVM_t_rcp26 <- c(SVM_t, thr)
      SVM_d_rcp26 <- TPR * (1 - pi)
      
      ### rcp45
      n_cells <- nrow(na.omit(values(rcp45))) 
      pi <- sum(values(SVM_rcp45 >= thr), na.rm = T) / n_cells
      SVM_e_rcp45 <- c(SVM_e, TPR)
      SVM_t_rcp45 <- c(SVM_t, thr)
      SVM_d_rcp45 <- TPR * (1 - pi)
      
      ### rcp60
      n_cells <- nrow(na.omit(values(rcp60))) 
      pi <- sum(values(SVM_rcp60 >= thr), na.rm = T) / n_cells
      SVM_e_rcp60 <- c(SVM_e, TPR)
      SVM_t_rcp60 <- c(SVM_t, thr)
      SVM_d_rcp60 <- TPR * (1 - pi)
      
      ### rcp85
      n_cells <- nrow(na.omit(values(rcp85))) 
      pi <- sum(values(SVM_rcp85 >= thr), na.rm = T) / n_cells
      SVM_e_rcp85 <- c(SVM_e, TPR)
      SVM_t_rcp85 <- c(SVM_t, thr)
      SVM_d_rcp85 <- TPR * (1 - pi)
      
      
      ### Saving partial outputs for the RCPs
      
      bioclim_Pout_rcp26 <- cbind(bioclim_Pout_rcp26, values(bioclim_rcp26))
      bioclim_Pout_rcp45 <- cbind(bioclim_Pout_rcp45, values(bioclim_rcp45))
      bioclim_Pout_rcp60 <- cbind(bioclim_Pout_rcp60, values(bioclim_rcp60))
      bioclim_Pout_rcp85 <- cbind(bioclim_Pout_rcp85, values(bioclim_rcp85))
      
      gower_Pout_rcp26   <- cbind(gower_Pout_rcp26,   values(gower_rcp26))
      gower_Pout_rcp45   <- cbind(gower_Pout_rcp45,   values(gower_rcp45))
      gower_Pout_rcp60   <- cbind(gower_Pout_rcp60,   values(gower_rcp60))
      gower_Pout_rcp85   <- cbind(gower_Pout_rcp85,   values(gower_rcp85))
      
      maxent_Pout_rcp26  <- cbind(maxent_Pout_rcp26,  values(maxent_rcp26))
      maxent_Pout_rcp45  <- cbind(maxent_Pout_rcp45,  values(maxent_rcp45))
      maxent_Pout_rcp60  <- cbind(maxent_Pout_rcp60,  values(maxent_rcp60))
      maxent_Pout_rcp85  <- cbind(maxent_Pout_rcp85,  values(maxent_rcp85))
      
      SVM_Pout_rcp45     <- cbind(SVM_Pout_rcp45,     values(SVM_rcp45))
      SVM_Pout_rcp26     <- cbind(SVM_Pout_rcp26,     values(SVM_rcp26))
      SVM_Pout_rcp60     <- cbind(SVM_Pout_rcp60,     values(SVM_rcp60))
      SVM_Pout_rcp85     <- cbind(SVM_Pout_rcp85,     values(SVM_rcp85))
      
      # CLOSE "i" ----
      # AOGCMs
    }
    
    # CLOSE "j"  ---- 
    # cross-validation
  } 
  
  ###.............. Calculating mean of partial models outputs
  bioclim_Pout_c_mean     <- apply(bioclim_Pout_c,     1, mean)
  bioclim_Pout_rcp26_mean <- apply(bioclim_Pout_rcp26, 1, mean)
  bioclim_Pout_rcp45_mean <- apply(bioclim_Pout_rcp45, 1, mean)
  bioclim_Pout_rcp60_mean <- apply(bioclim_Pout_rcp60, 1, mean)
  bioclim_Pout_rcp85_mean <- apply(bioclim_Pout_rcp85, 1, mean)
  
  gower_Pout_c_mean       <- apply(gower_Pout_c,     1, mean)
  gower_Pout_rcp26_mean   <- apply(gower_Pout_rcp26, 1, mean)
  gower_Pout_rcp45_mean   <- apply(gower_Pout_rcp45, 1, mean)
  gower_Pout_rcp60_mean   <- apply(gower_Pout_rcp60, 1, mean)
  gower_Pout_rcp85_mean   <- apply(gower_Pout_rcp85, 1, mean)
  
  maxent_Pout_c_mean      <- apply(maxent_Pout_c,     1, mean)
  maxent_Pout_rcp26_mean  <- apply(maxent_Pout_rcp26, 1, mean)
  maxent_Pout_rcp45_mean  <- apply(maxent_Pout_rcp45, 1, mean)
  maxent_Pout_rcp60_mean  <- apply(maxent_Pout_rcp60, 1, mean)
  maxent_Pout_rcp85_mean  <- apply(maxent_Pout_rcp85, 1, mean)
  
  SVM_Pout_c_mean         <- apply(SVM_Pout_c,     1, mean)
  SVM_Pout_rcp26_mean     <- apply(SVM_Pout_rcp26, 1, mean)
  SVM_Pout_rcp45_mean     <- apply(SVM_Pout_rcp45, 1, mean)
  SVM_Pout_rcp60_mean     <- apply(SVM_Pout_rcp60, 1, mean)
  SVM_Pout_rcp85_mean     <- apply(SVM_Pout_rcp85, 1, mean)
  
  
  ###.............. Saving data into the output objects
  output_current <- cbind(output_current, Bioclim = bioclim_Pout_c_mean, Gower = gower_Pout_c_mean, Maxent = maxent_Pout_c_mean, SVM = SVM_Pout_c_mean )
  
  output_rcp26   <- cbind(output_rcp26, Bioclim = bioclim_Pout_rcp26_mean, Gower = gower_Pout_rcp26_mean, Maxent = maxent_Pout_rcp26_mean, SVM = SVM_Pout_rcp26_mean )
  
  output_rcp45   <- cbind(output_rcp45, Bioclim = bioclim_Pout_rcp45_mean, Gower = gower_Pout_rcp45_mean, Maxent = maxent_Pout_rcp45_mean, SVM = SVM_Pout_rcp45_mean )
  
  output_rcp60   <- cbind(output_rcp60, Bioclim = bioclim_Pout_rcp60_mean, Gower = gower_Pout_rcp60_mean, Maxent = maxent_Pout_rcp60_mean, SVM = SVM_Pout_rcp60_mean )
  
  output_rcp85   <- cbind(output_rcp85, Bioclim = bioclim_Pout_rcp85_mean, Gower = gower_Pout_rcp85_mean, Maxent = maxent_Pout_rcp85_mean, SVM = SVM_Pout_rcp85_mean )
  
  
  ###.............. Inserting coords to outputs
  coords <- xyFromCell(current, 1:ncell(current))
  
  output_current <- cbind(coords,output_current)
  output_rcp26   <- cbind(coords,output_rcp26)
  output_rcp45   <- cbind(coords,output_rcp45)
  output_rcp60   <- cbind(coords,output_rcp60)
  output_rcp85   <- cbind(coords,output_rcp85)
  
  
  ###.............. Excluding NAs from outputs 
  output_current <- na.omit(output_current)
  output_rcp26   <- na.omit(output_rcp26)
  output_rcp45   <- na.omit(output_rcp45)
  output_rcp60   <- na.omit(output_rcp60)
  output_rcp85   <- na.omit(output_rcp85)
  
  
  ###.............. Evaluation data
  models_e <- data.frame(bioclim = bioclim_e, gower = gower_e, maxent = maxent_e, SVM = SVM_e)
  models_t <- data.frame(bioclim = bioclim_t, gower = gower_t, maxent = maxent_t, SVM = SVM_t)
  models_d <- data.frame(bioclim = bioclim_d, gower = gower_d, maxent = maxent_d, SVM = SVM_d)
  
  models_e_rcp26 <- data.frame(bioclim = bioclim_e_rcp26, gower = gower_e_rcp26, maxent = maxent_e_rcp26, SVM = SVM_e_rcp26)
  models_t_rcp26 <- data.frame(bioclim = bioclim_t_rcp26, gower = gower_t_rcp26, maxent = maxent_t_rcp26, SVM = SVM_t_rcp26)
  models_d_rcp26 <- data.frame(bioclim = bioclim_d_rcp26, gower = gower_d_rcp26, maxent = maxent_d_rcp26, SVM = SVM_d_rcp26)
  
  models_e_rcp45 <- data.frame(bioclim = bioclim_e_rcp45, gower = gower_e_rcp45, maxent = maxent_e_rcp45, SVM = SVM_e_rcp45)
  models_t_rcp45 <- data.frame(bioclim = bioclim_t_rcp45, gower = gower_t_rcp45, maxent = maxent_t_rcp45, SVM = SVM_t_rcp45)
  models_d_rcp45 <- data.frame(bioclim = bioclim_d_rcp45, gower = gower_d_rcp45, maxent = maxent_d_rcp45, SVM = SVM_d_rcp45)
  
  models_e_rcp60 <- data.frame(bioclim = bioclim_e_rcp60, gower = gower_e_rcp60, maxent = maxent_e_rcp60, SVM = SVM_e_rcp60)
  models_t_rcp60 <- data.frame(bioclim = bioclim_t_rcp60, gower = gower_t_rcp60, maxent = maxent_t_rcp60, SVM = SVM_t_rcp60)
  models_d_rcp60 <- data.frame(bioclim = bioclim_d_rcp60, gower = gower_d_rcp60, maxent = maxent_d_rcp60, SVM = SVM_d_rcp60)
  
  models_e_rcp85 <- data.frame(bioclim = bioclim_e_rcp85, gower = gower_e_rcp85, maxent = maxent_e_rcp85, SVM = SVM_e_rcp85)
  models_t_rcp85 <- data.frame(bioclim = bioclim_t_rcp85, gower = gower_t_rcp85, maxent = maxent_t_rcp85, SVM = SVM_t_rcp85)
  models_d_rcp85 <- data.frame(bioclim = bioclim_d_rcp85, gower = gower_d_rcp85, maxent = maxent_d_rcp85, SVM = SVM_d_rcp85)
  
  
  #\o/\o/\o/\o/\o/\o/\o/\o/\o/\o/   YOU SHALL... PASS!  \o/\o/\o/\o/\o/\o/\o/\o/\o/\o/
  
  alarm()
  return(list(c("output_current"  = output_current, 
                "output_rcp26"    = output_rcp26, 
                "output_rcp45"    = output_rcp45, 
                "output_rcp60"    = output_rcp60, 
                "output_rcp85"    = output_rcp85, 
                "TPR_c"           = models_e, 
                "Threshold_c"     = models_t, 
                "Pred_area_c"     = models_d,
                "TPR_rcp26"       = models_e_rcp26, 
                "Threshold_rcp26" = models_t_rcp26, 
                "Pred_area_rcp26" = models_d_rcp26,
                "TPR_rcp45"       = models_e_rcp45, 
                "Threshold_rcp45" = models_t_rcp45, 
                "Pred_area_rcp45" = models_d_rcp45,
                "TPR_rcp60"       = models_e_rcp60, 
                "Threshold_rcp60" = models_t_rcp60, 
                "Pred_area_rcp60" = models_d_rcp60,
                "TPR_rcp85"       = models_e_rcp85, 
                "Threshold_rcp85" = models_t_rcp85, 
                "Pred_area_rcp85" = models_d_rcp85
  )))
  
} # CLOSE "Multiple_ENMs"