multiple_ENMs <- function(occurrence, 
                          background, 
                          # predictor, # biotic variables from XP1.2
                          biovar_current,
                          biovar_rcp26,
                          biovar_rcp45,
                          biovar_rcp60,
                          biovar_rcp85,
                          cross_validation)
{
  
  ###.............. Reading the selected climatic variables
  
  ## Reading the selected current model
  current <- stack(list.files(biovar_current,  pattern = ".grd$", full.names = TRUE))
  
  ###.............. Reading ocurrence and background
  occur <- read.table(occurrence, sep = ";", h = T)
  back  <- read.table(background, sep = ";", h = T)
  
  ###.............. Creating objects for saving the results
  ## Predictive methods
  bioclim_c <- gower_c  <- maxent_c <- SVM_c <- stack()
  bioclim_rcp26 <- gower_rcp26 <- maxent_rcp26 <- SVM_rcp26 <- stack() 
  bioclim_rcp45 <- gower_rcp45 <- maxent_rcp45 <- SVM_rcp45 <- stack()
  bioclim_rcp60 <- gower_rcp60 <- maxent_rcp60 <- SVM_rcp60 <- stack()
  bioclim_rcp85 <- gower_rcp85 <- maxent_rcp85 <- SVM_rcp85 <- stack()
  
  ## Evaluation
  bioclim_e <- gower_e <- maxent_e <- SVM_e <- NULL # True presence rate (TPR)
  bioclim_t <- gower_t <- maxent_t <- SVM_t <- NULL # The highest threshold at which there is no omission
  bioclim_d <- gower_d <- maxent_d <- SVM_d <- NULL # Area predicted as presence
  
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
    # thr <- threshold(bioclim_eval, "no_omission")
    thr <- quantile(extract(bioclim_c[[j]], occur[,1:2]), 0.05) 
    TPR <- evaluate(p = testing[testing[, "pb"] == 1, -1], a = testing[testing[, "pb"] == 0, -1], model = bioclim_model, tr = thr)@TPR
    n_cells <- nrow(na.omit(values(current))) 
    pi <- sum(values(bioclim_c >= thr), na.rm = T) / n_cells # predicted area proportion.
    
    bioclim_e <- c(bioclim_e, TPR)
    bioclim_t <- c(bioclim_t, thr) 
    bioclim_d <- c(bioclim_d, TPR * (1 - pi))
    
    
    ### Gower -------------------------------------------------------------------------------
    
    ## Adjusting models
    gower_model <- domain(trainning[trainning[,"pb"] == 1, -1]) 
    
    ## Predicting
    gower_c <- stack(gower_c, predict(object = gower_model, x = current))
    
    ## Evaluating models
    gower_eval <- evaluate(p = testing[testing[, "pb"] == 1, -1], a = testing[testing[, "pb"] == 0, -1], model = gower_model)
    
    ## Saving evaluation metrics
    thr <- quantile(extract(gower_c[[j]], occur[,1:2]), 0.05) 
    TPR <- evaluate(p = testing[testing[, "pb"] == 1, -1], a = testing[testing[, "pb"] == 0, -1], model = gower_model, tr = thr)@TPR
    n_cells <- nrow(na.omit(values(current))) 
    pi <- sum(values(gower_c >= thr), na.rm = T) / n_cells
    
    gower_e <- c(gower_e, TPR)
    gower_t <- c(gower_t, thr)
    gower_d <- c(gower_d, TPR * (1 - pi))
    
    
    ### Maxent -------------------------------------------------------------------------------
    
    ## Adjusting models
    Sys.setenv(NOAWT = TRUE)
    maxent_model <- maxent(x = trainning[, -1], p = trainning[, 1])
    
    ## Predicting
    maxent_c <- stack(maxent_c, predict(object = maxent_model, x = current))
    
    ## Evaluating models
    maxent_eval <- evaluate(p = testing[testing[, "pb"] == 1, -1], a = testing[testing[, "pb"] == 0, -1], model = maxent_model)
    
    ## Saving evaluation metrics
    thr <- quantile(extract(maxent_c[[j]], occur[,1:2]), 0.05) 
    TPR <- evaluate(p = testing[testing[, "pb"] == 1, -1], a = testing[testing[, "pb"] == 0, -1], model = maxent_model, tr = thr)@TPR
    n_cells <- nrow(na.omit(values(current))) 
    pi <- sum(values(maxent_c >= thr), na.rm = T) / n_cells
    
    maxent_e <- c(maxent_e, TPR) 
    maxent_t <- c(maxent_t, thr)
    maxent_d <- c(maxent_d, TPR * (1 - pi))
    
    
    ### SVM ----------------------------------------------------------------------------------
    
    ## Adjusting models
    SVM_model <- ksvm(pb ~ ., data = trainning)
    
    ## Predicting
    SVM_c <- stack(SVM_c, predict(model = SVM_model, object = current)) 
    
    ## Evaluating models
    SVM_eval <- evaluate(p = testing[testing[, "pb"] == 1, -1], a = testing[testing[, "pb"] == 0, -1], model = SVM_model)
    
    ## Saving evaluation metrics
    thr <- quantile(extract(SVM_c[[j]], occur[,1:2]), 0.05) 
    TPR <- evaluate(p = testing[testing[, "pb"] == 1, -1], a = testing[testing[, "pb"] == 0, -1], model = SVM_model, tr = thr)@TPR
    n_cells <- nrow(na.omit(values(current))) 
    pi <- sum(values(SVM_c >= thr), na.rm = T) / n_cells
    
    SVM_e <- c(SVM_e, TPR)
    SVM_t <- c(SVM_t, thr)
    SVM_d <- c(SVM_d, TPR * (1 - pi))
    
    
    # ***************************************************************************************
    
    ###.............. Making predictions for the RCPs
    
   
    AOGCMs <- c(1, 2, 3)
    for (i in AOGCMs)
    {
      # OPEN "i" ----
      
      ### Reading the selected AOGCMs climatic models
      mdid <- paste0('.',1:3,'.grd$') # models identification.bio02.1 - ".1" stands for model 1.
      
      mdls <- lapply(mdid, function(x){stack(list.files(biovar_rcp26, pattern = x, full.names = TRUE))})
      rcp26  <- mdls[[i]] # now only the jth model will be extracted from the mdls.
      names(rcp26) <- names(current)
      
      mdls <- lapply(mdid, function(x){stack(list.files(biovar_rcp45, pattern = x, full.names = TRUE))})
      rcp45  <- mdls[[i]]
      names(rcp45) <- names(current)
      
      mdls  <- lapply(mdid, function(x){stack(list.files(biovar_rcp60, pattern = x, full.names = TRUE))})
      rcp60  <- mdls[[i]]
      names(rcp60) <- names(current)
      
      mdls  <- lapply(mdid, function(x){stack(list.files(biovar_rcp85, pattern = x, full.names = TRUE))})
      rcp85  <- mdls[[i]]
      names(rcp85) <- names(current)
      
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
      
      # CLOSE "i" ----
      # AOGCMs
    }
    
    # CLOSE "j"  ---- 
    # cross-validation
  } 
  
  ###..............  Saving partial outputs
  bioclim_Pout_c <- na.omit(values(bioclim_c))
  gower_Pout_c   <- na.omit(values(gower_c))
  maxent_Pout_c  <- na.omit(values(maxent_c))
  SVM_Pout_c     <- na.omit(values(SVM_c))
  
  bioclim_Pout_rcp26 <- na.omit(values(bioclim_rcp26))
  bioclim_Pout_rcp45 <- na.omit(values(bioclim_rcp45))
  bioclim_Pout_rcp60 <- na.omit(values(bioclim_rcp60))
  bioclim_Pout_rcp85 <- na.omit(values(bioclim_rcp85))
  
  gower_Pout_rcp26   <- na.omit(values(gower_rcp26))
  gower_Pout_rcp45   <- na.omit(values(gower_rcp45))
  gower_Pout_rcp60   <- na.omit(values(gower_rcp60))
  gower_Pout_rcp85   <- na.omit(values(gower_rcp85))
  
  maxent_Pout_rcp26  <- na.omit(values(maxent_rcp26))
  maxent_Pout_rcp45  <- na.omit(values(maxent_rcp45))
  maxent_Pout_rcp60  <- na.omit(values(maxent_rcp60))
  maxent_Pout_rcp85  <- na.omit(values(maxent_rcp85))
  
  SVM_Pout_rcp45     <- na.omit(values(SVM_rcp45))
  SVM_Pout_rcp26     <- na.omit(values(SVM_rcp26))
  SVM_Pout_rcp60     <- na.omit(values(SVM_rcp60))
  SVM_Pout_rcp85     <- na.omit(values(SVM_rcp85))
  
  # [MA] bioclim_Pout_c é um objeto que tem 20 prediçoes empilhadas em bioclim_c (cross-validation). Por isso eu tinha feito a média (bioclim_Pout_c_mean) e usado a média das 20 prediçoes para salvar em output_current.Não fez sentido para mim salvar bioclim_Pout_c no objeto que vou salvar no HD como meu resultado.
  
  ###.............. Saving the predictive methods into the output objects 
  output_current <- cbind(Bioclim = bioclim_Pout_c,     Gower = gower_Pout_c,     Maxent = maxent_Pout_c,     SVM = SVM_Pout_c )
  
  output_rcp26   <- cbind(Bioclim = bioclim_Pout_rcp26, Gower = gower_Pout_rcp26, Maxent = maxent_Pout_rcp26, SVM = SVM_Pout_rcp26 )
  
  output_rcp45   <- cbind(Bioclim = bioclim_Pout_rcp45, Gower = gower_Pout_rcp45, Maxent = maxent_Pout_rcp45, SVM = SVM_Pout_rcp45 )
  
  output_rcp60   <- cbind(Bioclim = bioclim_Pout_rcp60, Gower = gower_Pout_rcp60, Maxent = maxent_Pout_rcp60, SVM = SVM_Pout_rcp60 )
  
  output_rcp85   <- cbind(Bioclim = bioclim_Pout_rcp85, Gower = gower_Pout_rcp85, Maxent = maxent_Pout_rcp85, SVM = SVM_Pout_rcp85 )
  
  
  # [MA] do jeito que está.. nao estou utilizando essas médias em nenhum outro local. pq no Ensemble vc já fez a média ponderada mas já utilizando os dados do cbind acima. Por isso esse trecho só faz sentido para mim se for para salvar como output final.
  
  ###.............. Calculating mean of partial models outputs 
  # bioclim_Pout_c_mean     <- apply(bioclim_Pout_c,     1, function(x) sum(x*bioclim_d)/sum(bioclim_d))
  # bioclim_Pout_rcp26_mean <- apply(bioclim_Pout_rcp26, 1, function(x) sum(x*rep(bioclim_d, length(AOGCMs)))/sum(rep(bioclim_d, length(AOGCMs))))
  # bioclim_Pout_rcp45_mean <- apply(bioclim_Pout_rcp45, 1, function(x) sum(x*rep(bioclim_d, length(AOGCMs)))/sum(rep(bioclim_d, length(AOGCMs))))
  # bioclim_Pout_rcp60_mean <- apply(bioclim_Pout_rcp60, 1, function(x) sum(x*rep(bioclim_d, length(AOGCMs)))/sum(rep(bioclim_d, length(AOGCMs))))
  # bioclim_Pout_rcp85_mean <- apply(bioclim_Pout_rcp85, 1, function(x) sum(x*rep(bioclim_d, length(AOGCMs)))/sum(rep(bioclim_d, length(AOGCMs))))
  # 
  # gower_Pout_c_mean       <- apply(gower_Pout_c,     1, function(x) sum(x*gower_d)/sum(gower_d))
  # gower_Pout_rcp26_mean   <- apply(gower_Pout_rcp26, 1, function(x) sum(x*rep(gower_d, length(AOGCMs)))/sum(rep(gower_d, length(AOGCMs))))
  # gower_Pout_rcp45_mean   <- apply(gower_Pout_rcp45, 1, function(x) sum(x*rep(gower_d, length(AOGCMs)))/sum(rep(gower_d, length(AOGCMs))))
  # gower_Pout_rcp60_mean   <- apply(gower_Pout_rcp60, 1, function(x) sum(x*rep(gower_d, length(AOGCMs)))/sum(rep(gower_d, length(AOGCMs))))
  # gower_Pout_rcp85_mean   <- apply(gower_Pout_rcp85, 1, function(x) sum(x*rep(gower_d, length(AOGCMs)))/sum(rep(gower_d, length(AOGCMs))))
  # 
  # maxent_Pout_c_mean     <- apply(maxent_Pout_c,     1, function(x) sum(x*maxent_d)/sum(maxent_d))
  # maxent_Pout_rcp26_mean <- apply(maxent_Pout_rcp26, 1, function(x) sum(x*rep(maxent_d, length(AOGCMs)))/sum(rep(maxent_d, length(AOGCMs))))
  # maxent_Pout_rcp45_mean <- apply(maxent_Pout_rcp45, 1, function(x) sum(x*rep(maxent_d, length(AOGCMs)))/sum(rep(maxent_d, length(AOGCMs))))
  # maxent_Pout_rcp60_mean <- apply(maxent_Pout_rcp60, 1, function(x) sum(x*rep(maxent_d, length(AOGCMs)))/sum(rep(maxent_d, length(AOGCMs))))
  # maxent_Pout_rcp85_mean <- apply(maxent_Pout_rcp85, 1, function(x) sum(x*rep(maxent_d, length(AOGCMs)))/sum(rep(maxent_d, length(AOGCMs))))
  # 
  # SVM_Pout_c_mean     <- apply(SVM_Pout_c,     1, function(x) sum(x*SVM_d)/sum(SVM_d))
  # SVM_Pout_rcp26_mean <- apply(SVM_Pout_rcp26, 1, function(x) sum(x*rep(SVM_d, length(AOGCMs)))/sum(rep(SVM_d, length(AOGCMs))))
  # SVM_Pout_rcp45_mean <- apply(SVM_Pout_rcp45, 1, function(x) sum(x*rep(SVM_d, length(AOGCMs)))/sum(rep(SVM_d, length(AOGCMs))))
  # SVM_Pout_rcp60_mean <- apply(SVM_Pout_rcp60, 1, function(x) sum(x*rep(SVM_d, length(AOGCMs)))/sum(rep(SVM_d, length(AOGCMs))))
  # SVM_Pout_rcp85_mean <- apply(SVM_Pout_rcp85, 1, function(x) sum(x*rep(SVM_d, length(AOGCMs)))/sum(rep(SVM_d, length(AOGCMs))))
  # 
  
  ###.............. Generating Ensembles
  dStat_c <- c(bioclim_d, gower_d, maxent_d, SVM_d)
  dStat_fut <- rep(dStat_c, each= length(AOGCMs))
  
  ensemble_c     <- apply(output_current, 1, function(x) sum(x * dStat_c)   / sum(dStat_c))
  ensemble_rcp26 <- apply(output_rcp26,   1, function(x) sum(x * dStat_fut) / sum(dStat_fut))
  ensemble_rcp45 <- apply(output_rcp45,   1, function(x) sum(x * dStat_fut) / sum(dStat_fut))
  ensemble_rcp60 <- apply(output_rcp60,   1, function(x) sum(x * dStat_fut) / sum(dStat_fut))
  ensemble_rcp85 <- apply(output_rcp85,   1, function(x) sum(x * dStat_fut) / sum(dStat_fut))
  
  ###.............. Inserting coords to outputs
  coords <- na.omit(cbind(xyFromCell(current, 1:ncell(current)), values(current)))[,1:2]
  
  output_current <- cbind(coords,output_current)
  output_rcp26   <- cbind(coords,output_rcp26)
  output_rcp45   <- cbind(coords,output_rcp45)
  output_rcp60   <- cbind(coords,output_rcp60)
  output_rcp85   <- cbind(coords,output_rcp85)
  
  ###.............. Saving Ensemble
  FULLensemble       <- data.frame(coords, Ensemble_present= ensemble_c, Ensemble_rcp26= ensemble_rcp26, Ensemble_rcp45= ensemble_rcp45, Ensemble_rcp60= ensemble_rcp60, Ensemble_rcp85= ensemble_rcp85)
  
  ###.............. Evaluation data
  models_e <- data.frame(bioclim = bioclim_e, gower = gower_e, maxent = maxent_e, SVM = SVM_e)
  models_t <- data.frame(bioclim = bioclim_t, gower = gower_t, maxent = maxent_t, SVM = SVM_t)
  models_d <- data.frame(bioclim = bioclim_d, gower = gower_d, maxent = maxent_d, SVM = SVM_d)
  
  
  #\o/\o/\o/\o/\o/\o/\o/\o/\o/\o/   YOU SHALL... PASS!  \o/\o/\o/\o/\o/\o/\o/\o/\o/\o/
  
  return(list(c("output_current"  = output_current, 
                "output_rcp26"    = output_rcp26, 
                "output_rcp45"    = output_rcp45, 
                "output_rcp60"    = output_rcp60, 
                "output_rcp85"    = output_rcp85, 
                "TPR_c"           = models_e, 
                "Threshold_c"     = models_t, 
                "Pred_area_c"     = models_d,
                "Ensemble"        = FULLensemble)))
  
} # CLOSE "Multiple_ENMs"