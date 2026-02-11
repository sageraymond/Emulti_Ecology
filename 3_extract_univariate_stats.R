#
# AIM: compare models. 
#Compare univariate (linear) to null
#Compare univvariate (quadratic) to null
#Compare linear and quadratic
#
#
#
#

rm(list = ls())
gc()


library("data.table")
library("glmmTMB")
library("doParallel")
library("foreach")
library("doSNOW")
library("crayon")
library("ggplot2")
library("tidyr")
library(broom.mixed)
library(performance)

# Load data and prepare workspace -----------------------------------------

master_guide <- readRDS("builds/batchmods_nov/model_guide.Rds")
master_guide[!file.exists(model_path), ]

# master_guide <- master_guide[file.exists(model_path), ]

# 

# >>> Cast wide by model_type ---------------------------------------------
master_guide_wide <- dcast(master_guide,
                           ... ~ model_type,
                           value.var = c("model_call", "model_path", "model_id",
                                         "formula"))

master_guide_wide[model_complexity_comparison_ID == "model_complexity_id_361"]
master_guide_wide

master_guide_wide[model_complexity_comparison_ID == "model_complexity_id_370"]
master_guide_wide[duplicated(model_complexity_comparison_ID)]
if(nrow(master_guide_wide[duplicated(model_complexity_comparison_ID)]) > 0){ 
  cat(yellow("Stupid Lundy!"))
}else{ 
  cat(blue("good Lundy!"))
}
master_guide_wide


# >>> Tests ---------------------------------------------------------------
#' [ALL of these should have NULL models and univariate]

master_guide_wide[is.na(model_id_null_model)]
#' *hopefully 0 rows*


master_guide_wide[is.na(model_id_univariate)]
#' *hopefully 0 rows*

#' [If var != an urban variable, there should be an urbanization model id)]


#' [Should only be 1 row per model_complexity_comparison_ID]
master_guide_wide[duplicated(model_complexity_comparison_ID), ]

#' [Now check for ones that don't exist and exclude]
master_guide_wide <- master_guide_wide[file.exists(model_path_null_model) &
                                         file.exists(model_path_univariate), ]

# >>> Compare -------------------------------------------------------------
comps <- list()  # make list

for (i in 1:nrow(master_guide_wide)) {
  
  model_path_null_model <- master_guide_wide[i, model_path_null_model]
  model_path_univariate <- master_guide_wide[i, model_path_univariate]
  model_id <- master_guide_wide[i, model_complexity_comparison_ID]
  
  # Skip if model paths are missing or file(s) don't exist
  if (is.na(model_path_null_model) || is.na(model_path_univariate)) next
  if (!file.exists(model_path_null_model) || !file.exists(model_path_univariate)) next
  
  # Read models
  m1 <- readRDS(model_path_null_model)
  m2 <- readRDS(model_path_univariate)
  
  # Compare  models
  out <- anova(m1, m2)
  
  # Store comparison
  comps[[i]] <- data.table(
    model_complexity_comparison_ID = model_id,
    null_uni_chisq = out$Chisq[2],
    null_uni_p = out$`Pr(>Chisq)`[2]
  )
  
  cat(i, "/", nrow(master_guide_wide), "\r")
}

# Combine all results into one data.table
comps.dt <- rbindlist(comps, fill = TRUE)

# where did comparison fail?
comps.dt[is.na(null_uni_chisq)] #logically those nine dudes that didnt converge
# Ideally 0 rows


# >>> Merge into master guide wide ----------------------------------------

master_guide_wide.mrg <- merge(master_guide_wide,
                               comps.dt,
                               by = "model_complexity_comparison_ID",
                               all.x = T)
master_guide_wide.mrg



# Encapsulate some functions --------------------------------------
#Tidy function for model summary info
tidy_glmmTMB <- function(m) {
  
  # Get tidy output for fixed effects
  tidy_df <- tidy(m, effects = "fixed", conf.int = TRUE)
  setDT(tidy_df)
  
  # Extract model-level stats
  n_obs <- nobs(m)
  aic_val <- AIC(m)
  
  # Try to get R² values
  r2_vals <- tryCatch({
    r2(m)
  }, error = function(e) data.frame(R2_marginal = NA, R2_conditional = NA))
  
  R2_marginal <- r2_vals$R2_marginal[1]
  R2_conditional <- r2_vals$R2_conditional[1]
  
  # Initialize output
  out <- data.table(
    n_obs = n_obs,
    AIC = aic_val,
    R2_marginal = R2_marginal,
    R2_conditional = R2_conditional
  )
  
  # Extract main and quadratic terms
  if (nrow(tidy_df) >= 1) {
    # Assume the first fixed effect after intercept is the main predictor
    main_row <- tidy_df[term != "(Intercept)"][1]
    if (!is.null(main_row)) {
      out[, c("term_main", "beta_main", "p_main", "lwr_main", "upr_main") :=
            list(main_row$term, main_row$estimate, main_row$p.value,
                 main_row$conf.low, main_row$conf.high)]
    }
    
    # Check for a quadratic term (typically I(x^2))
    quad_row <- tidy_df[grepl("^I\\(.*\\^2\\)$", term)]
    if (nrow(quad_row) > 0) {
      quad_row <- quad_row[1]
      out[, c("term_quad", "beta_quad", "p_quad", "lwr_quad", "upr_quad") :=
            list(quad_row$term, quad_row$estimate, quad_row$p.value,
                 quad_row$conf.low, quad_row$conf.high)]
    }
  }
  
  return(out)
}


#' *This will be for model coefficients*
#tidy_glmmTMB <- function(m){
#  m.tidy <- tidy(m)
#  setDT(m.tidy)
  
#  if(length(unique(m.tidy$component)) > 1){
#    m.tidy[, key := paste(component, term, sep = ".")]
#  }else{
#    m.tidy[, key := term]
#  }
#  m.tidy
  
#  cis <- confint(m) |> as.data.frame()
#  cis$key <- row.names(cis)
#  setDT(cis)
#  setnames(cis, c("2.5 %", "97.5 %"), c("lwr_CI", "upper_CI"))
#  cis
  
#  cis[key %in% m.tidy$key, ]
#  m.tidy.mrg <- merge(m.tidy,
#                      cis,
#                      by = "key",
#                      all.x = T,
#                      all.y = F)
#  m.tidy.mrg
#  return(m.tidy.mrg)
  
#}


# Extract coefficients with CIs --------------------------------------------------------
ms <- lapply(master_guide_wide.mrg$model_path_univariate,
             FUN=readRDS)
ms.tidy <- lapply(ms,
                  FUN=tidy_glmmTMB)
ms.tidy
names(ms.tidy) <- master_guide_wide.mrg$model_id_univariate

ms.tidy <- rbindlist(ms.tidy, idcol = "model_id_univariate", fill = TRUE)
ms.tidy
unique(ms.tidy$beta_quad) #Good

#Inrtegrate into model guide
master_guide_wide.mrg1 <- merge(master_guide_wide.mrg,
                                ms.tidy,
                               by = "model_id_univariate",
                               all.x = T)
master_guide_wide.mrg1

#OK. That's all I need for this point. 

#saveRDS(master_guide_wide.mrg1, "builds/batchmods_nov/model_guide_with_model_stats.Rds")
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ---------------------------------

#Now pull out a csv that will be useful
master_guide_wide.mrg2 <- master_guide_wide.mrg1

master_guide_wide.mrg2 <- master_guide_wide.mrg2[, .(response,
                           term_main,
                           null_uni_chisq,
                           null_uni_p,
                           n_obs,
                           AIC,
                           R2_marginal,
                           R2_conditional,
                           beta_main,
                           p_main,
                           lwr_main,
                           upr_main,
                           term_quad,
                           beta_quad,
                           p_quad,
                           lwr_quad,
                           upr_quad)]

#Switch to dplyr. sorry Lundy
master_guide_wide.mrg2 <- as.data.frame(master_guide_wide.mrg2)

master_guide_wide.mrg2 <- master_guide_wide.mrg2 %>%
  dplyr::rename("Outcome" = response,
                "Predictor" = term_main,
                "LRT Chi" = null_uni_chisq,
                "LRT P" = null_uni_p,
                "N" = n_obs,
                "beta" = beta_main,
                "Marg. R2" = R2_marginal,
                "Con. R2" = R2_conditional,
                "P" = p_main,
                "Low CI" = lwr_main,
                "Upp. CI" = upr_main)

#write.csv(master_guide_wide.mrg2, "data/AllUnivariateNov21.csv")


#it will also be useful to have a table that gives me lrt results and delta AIC for lin vs quad models

