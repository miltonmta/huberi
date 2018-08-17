###............. Load Packages ----
load_pak <- function(pkg)
{
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

###.............. Processing the current climatic model
read_current <- function (dir)                                
{
  model_raw <- stack(list.files(dir,  pattern = ".bil$", full.names = TRUE))
  e <- extent(-122, -18, -56, 14) 
  model_e <- crop(model_raw, e)                              
  val <- getValues(model_e)                                  
  coord <- xyFromCell(model_e, 1:ncell(model_e))             
  model <- cbind(coord, val)                                 
  model <- na.omit(model)                                    
  return(model)
}

###.............. Processing the AOGCMs climatic models at all RCPs
read_rcp <- function (x)
{
  directories <- list.dirs( x, full.names = TRUE)[-1]
  e <- extent(-122, -18, -56, 14)
  models_list <- list()
  rcp <- NULL
  for (i in 1:length(directories))
  {
    models_raw       <- stack(list.files(directories[i],pattern = ".tif$", full.names = TRUE))
    models_e         <- crop( models_raw , e )
    val              <- values (models_e)
    coord            <- xyFromCell(models_e, 1:ncell(models_e))
    models           <- cbind(coord, val)
    models           <- na.omit(models)
    models_list[[i]] <- models_e
    rcp              <- abind (rcp, models, along = 3)
  }
  return(list("array" = rcp, "rasters" = models_list))
}

###.............. Extracting variables                         ----

create_var <- function(sp,name)
{
  sp_cell <- cellFromXY(current_select, sp[, -1])
  duplicated(sp_cell)
  sp_cell <- unique(sp_cell)
  sp_coord <- xyFromCell(current_select, sp_cell)
  sp_var <- raster::extract(current_select, sp_cell)
  sp_var <- na.omit(cbind(sp_coord, sp_var))
  write.table(sp_var, file = paste0("./data/occurrences/var_", as.factor(name), ".txt"), row.names = F, sep = ";")
  return(sp_var)
}

###.............. Creating background                         ----
create_back <- function(sp, name)
{
  coord <- rasterToPoints(current_select)[, 1:2]
  back_id <- sample(1:nrow(coord), nrow(sp))
  back <- extract(current_select, coord[back_id, ])
  back <- cbind(coord [back_id, ], back)
  write.table(back, paste0("./data/occurrences/back_", as.factor(name), ".txt"), row.names = F, sep = ";") 
  return(back)
}