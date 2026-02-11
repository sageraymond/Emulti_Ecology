
rm(list = ls())
gc()

library("data.table")
library("glmmTMB")
library("doParallel")
library("foreach")
library("doSNOW")
library("crayon")
library("ggplot2")

# groundhog.library(libs, groundhog.day)

# >>> Helper functions ----------------------------------------------------

prepare_cluster <- function(n){
  require("parallel")
  require("foreach")
  require("doSNOW")
  
  nCores <- parallel::detectCores() -1 
  cl <- makeCluster(nCores)
  registerDoSNOW(cl)
  
  # Progress bar
  pb <- txtProgressBar(max = n, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  
  ret <- list(opts, pb, cl)
  names(ret) <- c("options", "progress", "cluster")
  return(ret)
  
  cat("Pass 'x$options' to .opts in foreach;
      'x$progress' to setTxtProgressBar(x$progress, i);
      'x$cluster' to stopCluster(x$cluster) after foreach")
}

# Load guide --------------------------------------------------------------

guide <- readRDS("builds/batchmods_nov25/model_guide.Rds")

# Prepare data ------------------------------------------------------------
dat <- readRDS("builds/prepared_cluster_data.rds")
clust_by_cat <- readRDS("builds/cat_by_clust.rds")

dat1 <- left_join(dat, clust_by_cat, by = "ClusterID")
dat <- dat1
dat <- dat %>%
  dplyr::mutate(compost = Compost/ scatcount,
                garden = Garden / scatcount) %>%
  dplyr::select(-(c(Compost, Garden)))

dat <-setDT(dat)

#get rid cluster ID = NA
dat <- dat[!is.na(ClusterID)] #brings n down to 490 :(

#aNd replace NAs with meaningul zero
dat[is.na(No.Pos), No.Pos := 0]

# Run models -------------------------------------------------------
rerun_all <- TRUE #' [But only if you really hate them all]
if(rerun_all){
  file.remove(list.files("builds/batchmods_oct/models/", full.names = T))
}

# >>> Sequential version --------------------------------------------------

guide

m <- c()
# sub.dat <- c() #' [This might not be necessary]
i <- 1
#' *Note that some models won't run because there is insufficient N for their random effects*
working_guide <- guide[!file.exists(model_path), ] #this is if you want to run everybpdy

warnings <- list()
sub_dat <- c()

#' [Make sure to have all variables scaled inside loop:]
#' [We could be more clever and only scale the variables that are in the formula...But it's a pain in the butt]
#unique(guide[grepl("scaled", var)]$var)
#unique(guide[grepl("scaled", var)]$urbanization_var)

#' [the first time I run this I commented out the tryCatch business so we can see what errors are being thrown]

#for(i in 1:nrow(working_guide)){
#  tryCatch(
#    expr={
      
#      response <- working_guide$response[i]
#      data_filter <- working_guide$data[i]
      
#      if (data_filter == "all_scats") {
#        sub_dat <- dat
#      } else if (data_filter == "pos_scats") {
#        sub_dat <- dat[!is.na(dat[[response]]), ]
#      } else {
#        stop(paste("Unknown data filter:", data_filter, "in row", i))
#      }
      #' [Need to score the new variables too...]
#      sub_dat[, `:=` (Apple_scaled = scale(Apples),
#                      Berries_scaled = scale(Berries),
#                      Birdseed_scaled = scale(Birdseed),
#                      Vegetation_scaled = scale(Vegetation),
#                      Nat_Prey_scaled = scale(NatPrey),
#                      Garbage_scaled = scale (Garbage),
#                      Unk_Anthro_scaled = scale(Unk_Anthro),
#                      proplowincome_scaled = scale(proplowincome),
#                      prophighincome_scaled = scale(prophighincome),
#                      Dogs_scaled = scale(Dogs),
#                      dist_to_road_m_scaled = scale(dist_to_road_m),
#                      dist_to_Bldg_m_scaled = scale(dist_to_Bldg_m),
#                      dist_to_CityCenter_m_scaled = scale(dist_to_CityCenter_m),
#                      dist_to_CityBounary_m_scaled = scale(dist_to_CityBounary_m),
#                      dist_to_NatArea_m_scaled = scale(dist_to_NatArea_m),
#                      dist_to_Ravine_m_scaled = scale(dist_to_Ravine_m),
#                      dist_to_RivEdge_m_scaled = scale(dist_to_RivEdge_m),
#                      cover_Nat_25m_scaled = scale(cover_Nat_25m),
#                      cover_Anth_25m_scaled = scale(cover_Anth_25m),
#                      cover_NatWater_25m_scaled = scale(cover_NatWater_25m),
#                      cover_AnthGrass_25m_scaled = scale(cover_AnthGrass_25m),
 #                     cover_Nat_50m_scaled = scale(cover_Nat_50m),
#                      cover_Anth_50m_scaled = scale(cover_Anth_50m),
#                      cover_NatWater_50m_scaled = scale(cover_NatWater_50m),
#                      cover_AnthGrass_50m_scaled = scale(cover_AnthGrass_50m),
#                      cover_Nat_250m_scaled = scale(cover_Nat_250m),
#                      cover_Anth_250m_scaled = scale(cover_Anth_250m),
 #                     cover_NatWater_250m_scaled = scale(cover_NatWater_250m),
#                      cover_AnthGrass_250m_scaled = scale(cover_AnthGrass_250m),
#                      bldg_dens_25mH_scaled = scale(bldg_dens_25mH),
#                      bldg_dens_50mH_scaled = scale(bldg_dens_50mH),
#                      bldg_dens_250mH_scaled = scale(bldg_dens_250mH),
#                      road_dens_25mH_scaled = scale(road_dens_25mH),
#                      road_dens_50mH_scaled = scale(road_dens_50mH),
#                      road_dens_250mH_scaled = scale(road_dens_250mH),
#                      urban_index_pca_scaled = scale(urban_index_pca),
#                      decay_road_002_scaled = scale(decay_road_002),
#                      decay_road_0005_scaled = scale(decay_road_0005),
#                      decay_road_00015_scaled = scale(decay_road_00015),
#                      decay_bldg_002_scaled = scale(decay_bldg_002),
#                      decay_bldg_0005_scaled = scale(decay_bldg_0005),
#                      decay_bldg_00015_scaled = scale(decay_bldg_00015),
#                      decay_citycenter_002_scaled = scale(decay_citycenter_002),
#                      decay_citycenter_0005_scaled = scale(decay_citycenter_0005),
#                      decay_citycenter_00015_scaled = scale(decay_citycenter_00015),
#                      decay_citybound_002_scaled = scale(decay_citybound_002),
#                      decay_citybound_0005_scaled = scale(decay_citybound_0005),
#                      decay_citybound_00015_scaled = scale(decay_citybound_00015),
#                      decay_Ravine_002_scaled = scale(decay_Ravine_002),
#                      decay_Ravine_0005_scaled = scale(decay_Ravine_0005),
#                      decay_Ravine_00015_scaled = scale(decay_Ravine_00015),
#                      decay_RivEdge_002_scaled = scale(decay_RivEdge_002),
#                      decay_RivEdge_0005_scaled = scale(decay_RivEdge_0005),
#                      decay_RivEdge_00015_scaled = scale(decay_RivEdge_00015),
#                      decay_natarea_002_scaled = scale(decay_natarea_002),
#                      decay_natarea_0005_scaled = scale(decay_natarea_0005),
#                      decay_natarea_00015_scaled = scale(decay_natarea_00015),
#                      Prop.CompletedDegree_scaled = scale(Prop.CompletedDegree),
#                      Prop.NoCertDiplomadegree_scaled = scale(Prop.NoCertDiplomadegree),
#                      Prop.CompletedHighschool_scaled = scale(Prop.CompletedHighschool),
#                      Prop.Over65_scaled = scale(Prop.Over65),
#                      Prop.Under15_scaled = scale(Prop.Under15))]
     
      
       # setdiff(working_guide$var, names(sub_dat))
#      print(paste("Evaluating model", i))
#      print(working_guide[i, model_call])
      
#      m <- eval(parse(text = working_guide[i, ]$model_call))
      
#      if(!is.null(m)) saveRDS(m, working_guide[i, ]$model_path)
#      
#      m <- NULL # just in case...
#      
#    },
#    error=function(e){
#      warnings[[i]] <<- data.table(model_id = working_guide[i, ]$model_id,
#                                  error = e)
#      cat(red("error at"), i, "\r")
#    },
#    warning=function(w){ #' *there were some convergence warnings...spooky*
#      # Let's store them...And well I guess drop those models?
#      warnings[[i]] <<- data.table(model_id = working_guide[i, ]$model_id,
#                                  warning = w)
#      cat(blue(i))
#    }

  #)

  #cat(blue(i), magenta("/"), red(nrow(working_guide)), "\r")
#}

#names(warnings) <- guide$model_id

#working_guide <- guide[!file.exists(model_path), ]
#nrow(working_guide) # these are the models that didn't run

 #  TRY DIFFERENT APPROACH
   warnings <- list()
    
    working_guide <- guide[!file.exists(model_path), ]  # ONLY run thw models that dont exist
    
    for(i in seq_len(nrow(working_guide))) {
      tryCatch({
        # Extract row info
        response_var <- working_guide$response[i]
        predictor_var <- working_guide$var[i]
        model_family <- working_guide$model_family[i]
        formula_str <- working_guide$formula[i]
        model_path <- working_guide$model_path[i]
        model_id <- working_guide$model_id[i]
        
        # Subset data based on guide$exclusion (like "complete.cases(var1, var2)")
        exclusion_expr <- parse(text = working_guide$exclusion[i])
        sub_dat <- dat[eval(exclusion_expr), ]
        
        
   #        Scale variables inside sub_dat 
        sub_dat[, `:=` (
          propApple_scaled = as.vector(scale(propApple)),
          propBerries_scaled = as.vector(scale(propBerries)),
          propSeed_scaled = as.vector(scale(propSeed)),
          propVeg_scaled = as.vector(scale(propVeg)),
          propPrey_scaled = as.vector(scale(propPrey)),
          propGarbage_scaled = as.vector(scale(propGarbage)),
          propAnthro_scaled = as.vector(scale(propAnthro)),
          lowincome_scaled = as.vector(scale(lowincome)),
          highincome_scaled = as.vector(scale(highincome)),
          Dogs1_scaled = as.vector(scale(Dogs1)),
          dist_to_road_m_scaled = as.vector(scale(dist_to_road_m)),
          dist_to_Bldg_m_scaled = as.vector(scale(dist_to_Bldg_m)),
          dist_to_CityCenter_m_scaled = as.vector(scale(dist_to_CityCenter_m)),
          dist_to_CityBounary_m_scaled = as.vector(scale(dist_to_CityBounary_m)),
          dist_to_NatArea_m_scaled = as.vector(scale(dist_to_NatArea_m)),
          dist_to_Ravine_m_scaled = as.vector(scale(dist_to_Ravine_m)),
          dist_to_RivEdge_m_scaled = as.vector(scale(dist_to_RivEdge_m)),
          Prop_Nat_scaled = as.vector(scale(Prop_Nat)),
          Prop_Anth_scaled = as.vector(scale(Prop_Anth)),
          bldg_density_scaled = as.vector(scale(bldg_density)),
          road_density_scaled = as.vector(scale(rd_density)),
          decay_road_002_scaled = as.vector(scale(decay_road_002)),
          decay_road_0005_scaled = as.vector(scale(decay_road_0005)),
          decay_road_00015_scaled = as.vector(scale(decay_road_00015)),
          decay_bldg_002_scaled = as.vector(scale(decay_bldg_002)),
          decay_bldg_0005_scaled = as.vector(scale(decay_bldg_0005)),
          decay_bldg_00015_scaled = as.vector(scale(decay_bldg_00015)),
          decay_citycenter_002_scaled = as.vector(scale(decay_citycenter_002)),
          decay_citycenter_0005_scaled = as.vector(scale(decay_citycenter_0005)),
          decay_citycenter_00015_scaled = as.vector(scale(decay_citycenter_00015)),
          decay_citybound_002_scaled = as.vector(scale(decay_citybound_002)),
          decay_citybound_0005_scaled = as.vector(scale(decay_citybound_0005)),
          decay_citybound_00015_scaled = as.vector(scale(decay_citybound_00015)),
          decay_Ravine_002_scaled = as.vector(scale(decay_Ravine_002)),
          decay_Ravine_0005_scaled = as.vector(scale(decay_Ravine_0005)),
          decay_Ravine_00015_scaled = as.vector(scale(decay_Ravine_00015)),
          decay_RivEdge_002_scaled = as.vector(scale(decay_RivEdge_002)),
          decay_RivEdge_0005_scaled = as.vector(scale(decay_RivEdge_0005)),
          decay_RivEdge_00015_scaled = as.vector(scale(decay_RivEdge_00015)),
          decay_natarea_002_scaled = as.vector(scale(decay_natarea_002)),
          decay_natarea_0005_scaled = as.vector(scale(decay_natarea_0005)),
          decay_natarea_00015_scaled = as.vector(scale(decay_natarea_00015)),
          NoCertDiplomadegree_scaled = as.vector(scale(NoCertDiplomadegree)),
          Over65_scaled = as.vector(scale(Over65)),
          Under15_scaled = as.vector(scale(Under15))
        )]
        
        # Build and run model
        formula_obj <- as.formula(formula_str)
        model <- glmmTMB(formula_obj,
                         family = eval(parse(text = model_family)),
                         data = sub_dat)
        
        # Save the model
        saveRDS(model, model_path)
        cat("Model", model_id, "done\n")
        
      }, error = function(e) {
        warnings[[i]] <<- list(model_id = model_id, error = e$message)
        cat("Error in model", model_id, ":", e$message, "\n")
      })
    }
    
    
    
working_guide <- guide[!file.exists(model_path), ]
nrow(working_guide) # these are the models that didn't run



# >>> Test the models that didn't run -------------------------------------
working_guide[1, ]$model_call

sub_dat <- dat[eval(parse(text = working_guide[1, ]$exclusion))]
eval(parse(text = working_guide[1, ]$model_call))

sub_dat <- dat[eval(parse(text = working_guide[2, ]$exclusion))]
eval(parse(text = working_guide[2, ]$model_call))


# >>> Screen model objects that didn't converge ---------------------------
#' [didn't do this before...but should have ]
m$pdHess

m <- readRDS(guide[1, ]$model_path)
m$sdr$pdHess

guide <- guide[file.exists(model_path), ]

for(i in 1:nrow(guide)){
  m <- readRDS(guide[i, ]$model_path)
  
  guide[i, model_converged := m$sdr$pdHess]
  cat(i, "/", nrow(guide), "\r")
}

#' *should have done this in the loop above...Next time...*
guide[model_converged == FALSE, ]

# OK. 4 failed. I can live with this no problem

# >>> Parallel version ---------------------------------------------------
# Having a weird problem I've never had before: not finding objects in global environment...
# Strange...
#' clust_out <- prepare_cluster(n = nrow(guide))
#' #'
#' success <- foreach(i = 1:nrow(guide),
#'                    .options.snow = clust_out$options,
#'                    .errorhandling = "pass",
#'                    .packages = c("data.table", "glmmTMB")) %dopar% {
#'                     # for some reason foreach is not inheriting dat from env...
#' 
#'                      # return(guide[i, ]$model_call )
#'                      m <- eval(parse(text = guide[i, ]$model_call))
#'                      saveRDS(m, guide[i, ]$model_path)
#'                      #
#'                      setTxtProgressBar(clust_out$progress, i)
#'                    }
#' stopCluster(clust_out$cluster)
#' success
#' 
#' guide[!file.exists(model_path), ]
#' 
#' guide[file.exists(model_path), ]


