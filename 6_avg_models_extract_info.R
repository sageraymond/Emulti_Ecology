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
library(readr)

# 2. Bring in all models and average them---------------------------------------
#data

dat_diet_EM <- readRDS("builds/dat_for_modavg/dat_diet_EM.rds")
dat_diet_CT <- readRDS("builds/dat_for_modavg/dat_diet_CT.rds")
dat_social_EM <- readRDS("builds/dat_for_modavg/dat_social_EM.rds")
dat_social_CT <- readRDS("builds/dat_for_modavg/dat_social_CT.rds")
dat_urb_EM <- readRDS("builds/dat_for_modavg/dat_urb_EM.rds")
dat_urb_CT <- readRDS("builds/dat_for_modavg/dat_urb_CT.rds")


#Read in model averaged objects
#diet_em_avg <- readRDS("builds/model_avg_objects/diet_em.rds")
#diet_ct_avg <- readRDS("builds/model_avg_objects/diet_ct.rds")
#social_em_avg <- readRDS("builds/model_avg_objects/social_em.rds")
#social_ct_avg <- readRDS("builds/model_avg_objects/social_ct.rds")
#urb_em_avg <- readRDS("builds/model_avg_objects/urb_em.rds")
#urb_ct_avg <- readRDS("builds/model_avg_objects/urb_ct.rds")


#rEAD in
folders <- c(
  "builds/finalmods/diet/EM",
  "builds/finalmods/diet/CT",
  "builds/finalmods/social/EM",
  "builds/finalmods/social/CT",
  "builds/finalmods/urb/EM",
  "builds/finalmods/urb/CT"
)




# Load models into a named list
model_list_nested <- lapply(folders, function(folder) {
  files <- list.files(folder, pattern = "\\.rds$", full.names = TRUE)
  lapply(files, readRDS)
})

names(model_list_nested) <- gsub("^builds/finalmods/", "", folders)

m <- read_rds("builds/finalmods/diet/CT/1.rds")
m <- read_rds("builds/finalmods/diet/EM/1.rds")
m <- read_rds("builds/finalmods/social/CT/1.rds")
m <- read_rds("builds/finalmods/social/EM/1.rds")
m <- read_rds("builds/finalmods/urb/CT/1.rds")
m <- read_rds("builds/finalmods/urb/EM/1.rds")

r2(m, tolerance = 1e-8)



lengths(model_list_nested)
# Now average each group separately
models <- model_list_nested[[5]]
# names(models) <- paste("model_", seq(1:length(models)))

length(models)
model_averages <- lapply(model_list_nested, function(models) {
  model.avg(models, revised.var = TRUE)
})
names(model_averages) <- names(model_list_nested)

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

average_r2_by_folder <- function(model_list_nested) {
  r2_summary <- lapply(names(model_list_nested), function(folder_name) {
    models <- model_list_nested[[folder_name]]
    
    # Remove NULL or non-model elements
    valid_models <- Filter(function(m) !is.null(m) && inherits(m, "glmmTMB"), models)
    
    if(length(valid_models) == 0) {
      warning(paste("No valid models in folder:", folder_name))
      return(data.table(
        folder = folder_name,
        avg_R2_marginal = NA_real_,
        avg_R2_conditional = NA_real_,
        n_models = 0
      ))
    }
    
    r2_values <- lapply(valid_models, function(mod) {
      r2_vals <- tryCatch({
        r2(mod)
      }, error = function(e) {
        message(sprintf("Failed to calculate R2 for model in folder '%s': %s", folder_name, e$message))
        data.frame(R2_marginal = NA_real_, R2_conditional = NA_real_)
      })
      return(r2_vals)
    })
    
    r2_dt <- rbindlist(r2_values, fill = TRUE)
    
    avg_marginal <- mean(r2_dt$R2_marginal, na.rm = TRUE)
    avg_conditional <- mean(r2_dt$R2_conditional, na.rm = TRUE)
    
    data.table(
      folder = folder_name,
      avg_R2_marginal = avg_marginal,
      avg_R2_conditional = avg_conditional,
      n_models = length(valid_models)
    )
  })
  
  r2_summary_dt <- rbindlist(r2_summary)
  return(r2_summary_dt)
}


calc_vif <- function(model) {
  # Try to get VIF for fixed effects only
  tryCatch({
    # extract fixed effects terms
    fixef_vars <- names(lme4::fixef(model))
    # car::vif expects lm or glm; glmmTMB might work directly or use workaround
    # If glmmTMB, use a workaround by refitting with lme4::glmer or use vif from performance package
    # Let's try performance::check_collinearity (supports glmmTMB)
    vif_df <- performance::check_collinearity(model)
    mean(vif_df$VIF, na.rm = TRUE)
  }, error = function(e) NA)
}

calc_roc_auc <- function(model) {
  # Only works for binary logistic models with predictions and true response
  tryCatch({
    # Extract response and predicted values
    dat <- model.frame(model)
    response <- model.response(dat)
    # Only for binary response (factor or 0/1)
    if (!all(response %in% c(0,1))) return(NA)
    # Get predicted probabilities
    preds <- predict(model, type = "response")
    roc_obj <- pROC::roc(response, preds, quiet = TRUE)
    as.numeric(pROC::auc(roc_obj))
  }, error = function(e) NA)
}

#data liost for null model
data_list <- list(
  "diet/EM" = dat_diet_EM,
  "diet/CT" = dat_diet_CT,
  "social/EM" = dat_social_EM,
  "social/CT" = dat_social_CT,
  "urb/EM" = dat_urb_EM,
  "urb/CT" = dat_urb_CT
)


#Build null
build_null_model <- function(model, data) {
  tryCatch({
    # Extract formula and get LHS
    f <- formula(model)
    lhs <- deparse(f[[2]])
    
    # Find random effect terms
    rand_terms <- lme4::findbars(f)
    if (length(rand_terms) == 0) {
      message("Skipping null model: no random effects.")
      return(NULL)
    }
    
    # Build RHS for null model
    rhs_rand <- paste0(sapply(rand_terms, function(x) paste0("(", deparse(x), ")")), collapse = " + ")
    null_formula <- as.formula(paste(lhs, "~", rhs_rand))
    
    # Get model family
    fam <- model$call$family
    if (is.name(fam)) fam <- eval(fam)
    
    # Fit null model
    glmmTMB(null_formula, data = data, family = fam)
  }, error = function(e) {
    message("Failed to build null model: ", e$message)
    return(NULL)
  })
}

build_null_update <- function(full_model) {
  tryCatch({
    # Find random effect terms
    rand_terms <- lme4::findbars(formula(full_model))
    if (length(rand_terms) == 0) {
      message("Skipping null model: no random effects.")
      return(NULL)
    }
    
    # Construct formula using same data environment
    response <- deparse(formula(full_model)[[2]])
    rand_part <- paste0(sapply(rand_terms, function(x) paste0("(", deparse(x), ")")), collapse = " + ")
    new_formula <- as.formula(paste(response, "~", rand_part))
    
    # Use update() to rebuild on same data object
    updated_model <- update(full_model, formula. = new_formula, evaluate = TRUE)
    updated_model
  }, error = function(e) {
    message("Failed to build null model via update: ", e$message)
    return(NULL)
  })
}

#get lrt p value
lrt_summary <- rbindlist(lapply(names(model_list_nested), function(folder) {
  models <- model_list_nested[[folder]]
  
  p_vals <- sapply(models, function(m) {
    null_mod <- build_null_update(m)
    if (!is.null(null_mod)) {
      tryCatch({
        comp <- anova(m, null_mod, test = "LRT")
        comp$`Pr(>Chisq)`[2]
      }, error = function(e) {
        message("LRT failed in ", folder, ": ", e$message)
        NA_real_
      })
    } else {
      NA_real_
    }
  })
  
  data.table(
    folder = folder,
    n_models = length(models),
    n_success = sum(!is.na(p_vals)),
    avg_LRT_p = mean(p_vals, na.rm = TRUE)
  )
}))


#
extract_model_metrics <- function(model, data) {
  r2_vals <- tryCatch({
    r2_out <- performance::r2(model)
    list(marginal = r2_out$R2_marginal, conditional = r2_out$R2_conditional)
  }, error = function(e) list(marginal = NA, conditional = NA))
  
  vif_val <- calc_vif(model)
  auc_val <- calc_roc_auc(model)
  n_obs <- tryCatch(nobs(model), error = function(e) NA)
  
  
  data.table(
    R2_marginal = r2_vals$marginal,
    R2_conditional = r2_vals$conditional,
    VIF = vif_val,
    ROC_AUC = auc_val,
    n_obs = n_obs
  )
}

#OK and get ready to apply to your models

summarize_all_metrics <- function(model_list_nested, data_list, lrt_summary = NULL) {
  result <- rbindlist(lapply(names(model_list_nested), function(folder_name) {
    models <- model_list_nested[[folder_name]]
    data <- data_list[[folder_name]]
    
    metrics_list <- lapply(models, function(mod) extract_model_metrics(mod, data))
    metrics_dt <- rbindlist(metrics_list)
    
    avg_metrics <- metrics_dt[, lapply(.SD, mean, na.rm = TRUE)]
    avg_metrics[, n_models := .N, by = NULL]
    avg_metrics[, folder := folder_name]
    
    avg_metrics
  }))
  
  # Merge with lrt_summary if provided
  if (!is.null(lrt_summary)) {
    result <- merge(result, lrt_summary[, .(folder, avg_LRT_p)], by = "folder", all.x = TRUE)
  } else {
    result[, avg_LRT_p := NA_real_]
  }
  
  setcolorder(result, c("folder", "n_models", "n_obs", "R2_marginal", "R2_conditional", "VIF", "ROC_AUC", "avg_LRT_p"))
  return(result[])
}


# apply to models
metrics_summary <- summarize_all_metrics(model_list_nested, data_list, lrt_summary)
metrics_summary


#Now I am going to save these things as CSV files
metrics_summary <- as.data.frame(metrics_summary)
metrics_summary <- metrics_summary %>%
  dplyr::select(-(n_models))

#write.csv(metrics_summary, file = "figures/model_average_performance_Nov21.csv")

summary_list <- as.data.frame(summary_list)
summary_list <- summary_list %>%
  dplyr::select(term, estimate, p_value, conf_low, conf_high, model_set)

#write.csv(summary_list, file = "figures/model_average_summary_Nov21.csv")



#Also save these model averaged objects
names(model_averages)
diet_em_avg <- model_averages[["diet/EM"]]
diet_ct_avg <- model_averages[["diet/CT"]]
social_em_avg <- model_averages[["social/EM"]]
social_ct_avg <- model_averages[["social/CT"]]
urb_em_avg <- model_averages[["urb/EM"]]
urb_ct_avg <- model_averages[["urb/CT"]]


saveRDS(diet_em_avg, file = file.path("builds/model_avg_objects/diet_em.rds"))
saveRDS(diet_ct_avg, file = file.path("builds/model_avg_objects/diet_ct.rds"))
saveRDS(social_em_avg, file = file.path("builds/model_avg_objects/social_em.rds"))
saveRDS(social_ct_avg, file = file.path("builds/model_avg_objects/social_ct.rds"))
saveRDS(urb_em_avg, file = file.path("builds/model_avg_objects/urb_em.rds"))
saveRDS(urb_ct_avg, file = file.path("builds/model_avg_objects/urb_ct.rds"))

# Save model averages
saveRDS(model_averages, file = "builds/model_avg_objects/model_averages.rds")

# Save original model list
saveRDS(model_list_nested, file = "builds/model_avg_objects/model_list_nested.rds")

