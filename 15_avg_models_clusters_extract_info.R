#
# AIM: average models
#Exttract summary info
#
#
#
#

rm(list = ls())
gc()

# Load data and prepare workspace -----------------------------------------
library(data.table)
library(ggplot2)
library(broom.mixed)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(performance)
library(car)
library(pROC)
library(MuMIn)


# 2. Bring in all models and average them---------------------------------------
#data

dat_diet_EM <- readRDS("builds/dat_diet_EM.cluster.rds")
dat_social_EM <- readRDS("builds/dat_social_EM.cluster.rds")
dat_social_CT <- readRDS("builds/dat_social_CT.cluster.rds")
dat_urb_EM <- readRDS("builds/dat_urb_EM.cluster.rds")
dat_urb_CT <- readRDS("builds/dat_urb_CT.cluster.rds")



#Read in model averaged objects
#diet_em_avg <- readRDS("builds/model_avg_objects/diet_em.rds")
#diet_ct_avg <- readRDS("builds/model_avg_objects/diet_ct.rds")
#social_em_avg <- readRDS("builds/model_avg_objects/social_em.rds")
#social_ct_avg <- readRDS("builds/model_avg_objects/social_ct.rds")
#urb_em_avg <- readRDS("builds/model_avg_objects/urb_em.rds")
#urb_ct_avg <- readRDS("builds/model_avg_objects/urb_ct.rds")


#rEAD in
folders <- c(
  "builds/finalmods/cluster/diet/EM",
  "builds/finalmods/cluster/social/EM",
  "builds/finalmods/cluster/social/CT",
  "builds/finalmods/cluster/urb/EM",
  "builds/finalmods/cluster/urb/CT"
)




# Load models into a named list
model_list_nested <- lapply(folders, function(folder) {
  files <- list.files(folder, pattern = "\\.cluster.rds$", full.names = TRUE)
  lapply(files, readRDS)
})
names(model_list_nested) <- gsub("^builds/finalmods/cluster/", "", folders)



# Now average each group separately
model_averages <- lapply(model_list_nested, function(models) {
  model.avg(models, revised.var = TRUE)
})


#3. make a function to pull out the shit that you want
extract_modelavg_info <- function(model_avg_obj, name = NA) {
  
  coefs <- summary(model_avg_obj)$coefmat.full
  coefs <- as.data.table(coefs, keep.rownames = "term")
  setnames(coefs, old = c("Estimate", "Std. Error", "z value", "Pr(>|z|)"),
           new = c("estimate", "std_error", "z_value", "p_value"))
  
  ci <- confint(model_avg_obj)
  ci <- as.data.table(ci, keep.rownames = "term")
  setnames(ci, c("2.5 %", "97.5 %"), c("conf_low", "conf_high"))
  
  out <- merge(coefs, ci, by = "term", all.x = TRUE)
  
  out[, model_set := name]
  out[, n_models := length(model_avg_obj$models)]  # Number of models
  
  
  base_model <- model_avg_obj$models[[1]]$model
  
  return(out[])
}
#And apply it to my models

summary_list <- rbindlist(
  lapply(names(model_averages), function(nm) {
    extract_modelavg_info(model_averages[[nm]], name = nm)
  }),
  fill = TRUE
)

summary_list #This has beta, p, and Cis. This is a good start


#4. build functions to etxcract vif, roc auc, r2 values

pseudo_adj_r2 <- function(model) {
  tryCatch({
    # Number of observations
    n <- nobs(model)
    
    # Number of fixed-effect predictors (subtract intercept)
    p <- length(fixef(model)$cond) - 1
    
    # Log-likelihood of full model
    ll_full <- logLik(model)
    
    # Log-likelihood of null model (intercept only)
    null_model <- update(model, . ~ 1)
    ll_null <- logLik(null_model)
    
    # McFadden pseudo-R²
    r2 <- 1 - as.numeric(ll_full / ll_null)
    
    # Adjusted pseudo-R²
    r2_adj <- 1 - (1 - r2) * (n - 1) / (n - p - 1)
    
    return(r2_adj)
  }, error = function(e) {
    message("Failed to calculate adjusted pseudo-R²: ", e$message)
    return(NA_real_)
  })
}

# dDefine the main function
average_r2_by_folder <- function(model_list_nested) {
  r2_summary <- lapply(names(model_list_nested), function(folder_name) {
    models <- model_list_nested[[folder_name]]
    
    # Keep only valid glmmTMB models
    valid_models <- Filter(function(m) !is.null(m) && inherits(m, "glmmTMB"), models)
    
    if (length(valid_models) == 0) {
      warning(paste("No valid models in folder:", folder_name))
      return(data.table(
        folder = folder_name,
        avg_adj_R2 = NA_real_,
        n_models = 0
      ))
    }
    
    # Compute adjusted pseudo-R² for each model
    r2_values <- sapply(valid_models, pseudo_adj_r2)
    
    # Average across models
    avg_adj_r2 <- mean(r2_values, na.rm = TRUE)
    
    data.table(
      folder = folder_name,
      avg_adj_R2 = avg_adj_r2,
      n_models = length(valid_models)
    )
  })
  
  r2_summary_dt <- rbindlist(r2_summary)
  return(r2_summary_dt)
}

#Apply
r2_summary <- average_r2_by_folder(model_list_nested)
print(r2_summary)

#Now VIF
average_vif_by_folder <- function(model_list_nested) {
  vif_summary <- lapply(names(model_list_nested), function(folder_name) {
    models <- model_list_nested[[folder_name]]
    
    valid_models <- Filter(function(m) !is.null(m) && inherits(m, "glmmTMB"), models)
    
    if (length(valid_models) == 0) {
      warning(paste("No valid models in folder:", folder_name))
      return(data.table(
        folder = folder_name,
        avg_vif = NA_real_,
        n_models = 0
      ))
    }
    
    vif_values <- sapply(valid_models, function(mod) {
      tryCatch({
        collinearity <- performance::check_collinearity(mod)
        mean(collinearity$VIF, na.rm = TRUE)
      }, error = function(e) {
        message(sprintf("Failed to calculate VIF for model in '%s': %s", folder_name, e$message))
        NA_real_
      })
    })
    
    data.table(
      folder = folder_name,
      avg_vif = mean(vif_values, na.rm = TRUE),
      n_models = length(valid_models)
    )
  })
  
  rbindlist(vif_summary)
}

# Apply
vif_summary <- average_vif_by_folder(model_list_nested)
print(vif_summary)



#OK. now TRY to compare using LRT


#data liost for null model
data_list <- list(
  "cluster/diet/EM" = dat_diet_EM,
  "cluster/social/EM" = dat_social_EM,
  "cluster/social/CT" = dat_social_CT,
  "cluster/urb/EM" = dat_urb_EM,
  "cluster/urb/CT" = dat_urb_CT
)


#Build null
build_null_update <- function(full_model) {
  tryCatch({
    f <- formula(full_model)
    rand_terms <- lme4::findbars(f)
    response <- deparse(f[[2]])
    
    if (length(rand_terms) == 0) {
      null_formula <- as.formula(paste(response, "~ 1"))
    } else {
      rand_part <- paste0(sapply(rand_terms, function(x) paste0("(", deparse(x), ")")), collapse = " + ")
      null_formula <- as.formula(paste(response, "~", rand_part))
    }
    
    update(full_model, formula. = null_formula, evaluate = TRUE)
    
  }, error = function(e) {
    message("Failed to build null model: ", e$message)
    return(NULL)
  })
}

# Compute average LRT p-values per folder
average_lrt_by_folder <- function(model_list_nested) {
  rbindlist(lapply(names(model_list_nested), function(folder) {
    models <- model_list_nested[[folder]]
    
    # Keep only glmmTMB models
    valid_models <- Filter(function(m) !is.null(m) && inherits(m, "glmmTMB"), models)
    
    if (length(valid_models) == 0) {
      warning("No valid models in folder: ", folder)
      return(data.table(folder = folder, n_models = 0, n_success = 0, avg_LRT_p = NA_real_))
    }
    
    # Run LRT vs null
    p_vals <- sapply(valid_models, function(m) {
      null_mod <- build_null_update(m)
      if (!is.null(null_mod)) {
        tryCatch({
          comp <- anova(m, null_mod, test = "LRT")
          comp$`Pr(>Chisq)`[2]
        }, error = function(e) {
          message("LRT failed for a model in ", folder, ": ", e$message)
          NA_real_
        })
      } else NA_real_
    })
    
    data.table(
      folder = folder,
      n_models = length(valid_models),
      n_success = sum(!is.na(p_vals)),
      avg_LRT_p = mean(p_vals, na.rm = TRUE)
    )
  }))
}

# Example usage
lrt_summary <- average_lrt_by_folder(model_list_nested)
print(lrt_summary)


#OK. Now get ready to save things
summary_list <- as.data.frame(summary_list)
summary_list <- summary_list %>%
  dplyr::select(term, estimate, p_value, conf_low, conf_high, model_set)

#write.csv(summary_list, file = "figures/model_average_summary_cluster_Nov26.csv")

#Now combine vif, lrt, and r2
metricssummary <- merge(r2_summary, vif_summary, by = "folder", all = TRUE)
metricssummary <- merge(metricssummary, lrt_summary, by = "folder", all = TRUE)

#Now I am going to save these things as CSV files
metricssummary <- as.data.frame(metricssummary)
metricssummary <- metricssummary %>%
  dplyr::select(c(folder, avg_adj_R2, avg_vif, avg_LRT_p))

#write.csv(metricssummary, file = "figures/model_average_performance_cluster_Nov26.csv")




#Also save these model averaged objects

diet_em_avg <- model_averages[[1]]
social_em_avg <- model_averages[[2]]
social_ct_avg <- model_averages[[3]]
urb_em_avg <- model_averages[[4]]
urb_ct_avg <- model_averages[[5]]

saveRDS(diet_em_avg, file = file.path("builds/model_avg_objects_cluster/diet_em.rds"))
saveRDS(social_em_avg, file = file.path("builds/model_avg_objects_cluster/social_em.rds"))
saveRDS(social_ct_avg, file = file.path("builds/model_avg_objects_cluster/social_ct.rds"))
saveRDS(urb_em_avg, file = file.path("builds/model_avg_objects_cluster/urb_em.rds"))
saveRDS(urb_ct_avg, file = file.path("builds/model_avg_objects_cluster/urb_ct.rds"))

# Save model averages
saveRDS(model_averages, file = "builds/model_avg_objects_cluster/model_averages.rds")

# Save original model list
saveRDS(model_list_nested, file = "builds/model_avg_objects_cluster/model_list_nested.rds")

#Dave dataframes
saveRDS(dat_diet_EM, file = file.path("builds/dat_for_modavg_cluster/diet_EM.rds"))
saveRDS(dat_social_EM, file = file.path("builds/dat_for_modavg_cluster/social_EM.rds"))
saveRDS(dat_social_CT, file = file.path("builds/dat_for_modavg_cluster/social_CT.rds"))
saveRDS(dat_urb_EM, file = file.path("builds/dat_for_modavg_cluster/urb_EM.rds"))
saveRDS(dat_urb_CT, file = file.path("builds/dat_for_modavg_cluster/urb_CT.rds"))
