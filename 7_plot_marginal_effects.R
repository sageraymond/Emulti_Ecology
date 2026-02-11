#
# AIM: 
#get effect sizes by unscaled everybody and extracting coefficients
#Plot maginal effects using unscaled data. 
#
#

rm(list = ls())
gc()

# Load data and prepare workspace -----------------------------------------
library(data.table)
library(ggplot2)
library(ggeffects)
library(stats)
library(ggpubr)
library(broom)
library(dplyr)

model_averages <- readRDS("builds/model_avg_objects/model_averages.rds")
model_list_nested <- readRDS("builds/model_avg_objects/model_list_nested.rds")

#Read in model averaged objects
diet_em_avg <- readRDS("builds/model_avg_objects/diet_em.rds")
diet_ct_avg <- readRDS("builds/model_avg_objects/diet_ct.rds")
social_em_avg <- readRDS("builds/model_avg_objects/social_em.rds")
social_ct_avg <- readRDS("builds/model_avg_objects/social_ct.rds")
urb_em_avg <- readRDS("builds/model_avg_objects/urb_em.rds")
urb_ct_avg <- readRDS("builds/model_avg_objects/urb_ct.rds")

#read in the data 
diet_em_dat <- readRDS("builds/New_data_for_modavg/dat_diet_EM.rds")
diet_ct_dat <- readRDS("builds/New_data_for_modavg/dat_diet_CT.rds")
social_em_dat <- readRDS("builds/New_data_for_modavg/dat_social_EM.rds")
social_ct_dat <- readRDS("builds/New_data_for_modavg/dat_social_CT.rds")
urb_em_dat <- readRDS("builds/New_data_for_modavg/dat_urb_EM.rds")
urb_ct_dat <- readRDS("builds/New_data_for_modavg/dat_urb_CT.rds")



# Encapsulate some functions --------------------------------------

# Create marginal prediction manually
predict_marginal_effect_unscaled <- function(model_list, model_data, unscaled_var, scaled_var, n_points = 50) {
  # Step 1: Make sequence over *unscaled* variable
  pred_seq_unscaled <- seq(min(model_data[[unscaled_var]], na.rm = TRUE),
                           max(model_data[[unscaled_var]], na.rm = TRUE),
                           length.out = n_points)
  
  # Step 2: Scale it using the same transformation used in modeling
  mean_val <- mean(model_data[[unscaled_var]], na.rm = TRUE)
  sd_val <- sd(model_data[[unscaled_var]], na.rm = TRUE)
  pred_seq_scaled <- (pred_seq_unscaled - mean_val) / sd_val
  
  # Step 3: Build grid using *scaled* values for prediction
  base_data <- as.data.frame(lapply(model_data, function(x) {
    if (is.numeric(x)) median(x, na.rm = TRUE)
    else if (is.factor(x)) levels(x)[1]
    else x[1]
  }))
  
  grid <- base_data[rep(1, n_points), ]
  grid[[scaled_var]] <- pred_seq_scaled  # use scaled for prediction
  
  # Step 4: Predict using each model
  preds_list <- lapply(model_list, function(mod) {
    tryCatch({
      pred <- predict(mod, newdata = grid, type = "response", se.fit = TRUE)
      data.frame(fit = pred$fit, se.fit = pred$se.fit)
    }, error = function(e) {
      data.frame(fit = rep(NA, n_points), se.fit = rep(NA, n_points))
    })
  })
  
  pred_matrix <- do.call(cbind, lapply(preds_list, function(p) p$fit))
  se_matrix   <- do.call(cbind, lapply(preds_list, function(p) p$se.fit))
  
  predicted <- rowMeans(pred_matrix, na.rm = TRUE)
  se_avg <- apply(se_matrix, 1, function(x) sqrt(mean(x^2, na.rm = TRUE)))
  
  out <- data.frame(
    unscaled_predictor = pred_seq_unscaled,
    predicted_response = predicted,
    conf_low = predicted - 1.96 * se_avg,
    conf_high = predicted + 1.96 * se_avg
  )
  
  return(out)
}


#plotting function
plot_marginal_effect_unscaled <- function(
    marginal_df,
    xlab = "Unscaled Predictor",
    model_name = NULL,
    line_color = "blue",
    ci_fill = "lightblue",
    line_size = 1.2
) {
  ggplot(marginal_df, aes(x = unscaled_predictor, y = predicted_response)) +
    geom_ribbon(aes(ymin = conf_low, ymax = conf_high), fill = ci_fill, alpha = 0.3) +  # <- draw first
    geom_line(color = line_color, size = line_size) +                                    # <- draw second
    labs(
      x = xlab,
      y = "Predicted Response",
      title = paste("Marginal Effect of", xlab, "in", model_name)
    ) +
    plot_theme()
}


#Define themes for Diet, social, and whatever
plot_theme <- function(base_size = 12, base_family = "") {
  theme_classic() +
    theme(
      plot.title = element_text(size = 12, face = "plain"),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 12),
      #panel.grid.major = element_line(color = "gray80"),
      #panel.grid.minor = element_blank(),
      #legend.position = "bottom",
      #legend.title = element_text(size = 10),
      #legend.text = element_text(size = 9),
      legend = element_blank())
}

diet_plot_style <- function(line_color = "#50C5B7", ci_fill = "#50C5B7", line_size = 2) {
  list(
    geom_ribbon(aes(ymin = conf_low, ymax = conf_high), fill = ci_fill, alpha = 0.3),
    geom_line(color = line_color, size = line_size),
    plot_theme()
  )
}

diet_plot_style_2 <- function(line_color = "#50C5B7", ci_fill = "#50C5B7", line_size = 2) {
  list(
    geom_ribbon(aes(ymin = conf_low, ymax = conf_high), fill = ci_fill, alpha = 0.3),
    geom_line(color = line_color, size = line_size, linetype = "dashed"),
    plot_theme()
  )
}

social_plot_style <- function(line_color = "#9CEC5B", ci_fill = "#9CEC5B", line_size = 2) {
  list(
    geom_ribbon(aes(ymin = conf_low, ymax = conf_high), fill = ci_fill, alpha = 0.3),
    geom_line(color = line_color, size = line_size),
    plot_theme()
  )
}

social_plot_style_2 <- function(line_color = "#9CEC5B", ci_fill = "#9CEC5B", line_size = 2) {
  list(
    geom_ribbon(aes(ymin = conf_low, ymax = conf_high), fill = ci_fill, alpha = 0.3),
    geom_line(color = line_color, size = line_size, linetype = "dashed"),
    plot_theme()
  )
}


urb_plot_style <- function(line_color = "#6184D8", ci_fill = "#6184D8", line_size = 2) {
  list(
    geom_ribbon(aes(ymin = conf_low, ymax = conf_high), fill = ci_fill, alpha = 0.3),
    geom_line(color = line_color, size = line_size),
    plot_theme()
  )
}

urb_plot_style_2 <- function(line_color = "#6184D8", ci_fill = "#6184D8", line_size = 2) {
  list(
    geom_ribbon(aes(ymin = conf_low, ymax = conf_high), fill = ci_fill, alpha = 0.3),
    geom_line(color = line_color, size = line_size, linetype = "dashed"),
    plot_theme()
  )
}



#Try it
model_name <- "urb/CT"
model_data <- readRDS("builds/New_data_for_modavg/dat_urb_CT.rds")
model_list <- model_list_nested[[model_name]]

# Unscaled and scaled variable names
unscaled_var <- "dist_to_NatArea_m"
scaled_var <- "dist_to_NatArea_m_scaled"

marginal_df <- predict_marginal_effect_unscaled(
  model_list = model_list,
  model_data = model_data,
  unscaled_var = unscaled_var,
  scaled_var = scaled_var
)

plot_marginal_effect_unscaled(marginal_df, xlab = "YYY (unscaled)", model_name = model_name)


#OK, that seens to be working, which is freking exciting
#let's see if we can loop this

#Sytart by kind of mapping everything in a list
# A named list mapping model identifiers to their model list, dataset, and predictors
model_info <- list(
  "diet/EM" = list(
    model_list = model_list_nested[["diet/EM"]],
    data = diet_em_dat,
    vars = list(Berries = "Berries_scaled",
                Nat_Prey = "Nat_Prey_scaled")
  ),
  "diet/CT" = list(
    model_list = model_list_nested[["diet/CT"]],
    data = diet_ct_dat,
    vars = list(Berries = "Berries_scaled",
                Unk_Anthro = "Unk_Anthro_scaled")
  ),
  "social/EM" = list(
    model_list = model_list_nested[["social/EM"]],
    data = social_em_dat,
    vars = list(prophighincome = "prophighincome_scaled", Prop.Under15 = "Prop.Under15_scaled")
  ),
  "social/CT" = list(
    model_list = model_list_nested[["social/CT"]],
    data = social_ct_dat,
    vars = list(Prop.Over65 = "Prop.Over65_scaled", prophighincome = "prophighincome_scaled",
                proplowincome = "proplowincome_scaled")
  ),
  "urb/EM" = list(
    model_list = model_list_nested[["urb/EM"]],
    data = urb_em_dat,
    vars = list(cover_Anth_250m = "cover_Anth_250m_scaled",
                decay_road_0005 = "decay_road_0005_scaled",
                dist_to_Ravine_m = "dist_to_Ravine_m_scaled",
                decay_natarea_0005 = "decay_natarea_0005_scaled",
                road_dens_250mH = "road_dens_250mH_scaled")
  ),
  "urb/CT" = list(
    model_list = model_list_nested[["urb/CT"]],
    data = urb_ct_dat,
    vars = list(bldg_dens_25mH = "bldg_dens_25mH_scaled",
                decay_Ravine = "decay_Ravine_002_scaled",
                dist_to_NatArea_m = "dist_to_NatArea_m_scaled",
                road_dens_250mH = "road_dens_250mH_scaled")
  )
)


# Create a list to hold plots (optional, for later saving or display)
all_plots <- list()

# Loop over each model
for (model_name in names(model_info)) {
  info <- model_info[[model_name]]
  model_list <- info$model_list
  model_data <- info$data
  vars <- info$vars
  
  # Loop over each variable pair
  for (unscaled_var in names(vars)) {
    scaled_var <- vars[[unscaled_var]]
    
    # Predict
    marginal_df <- predict_marginal_effect_unscaled(
      model_list = model_list,
      model_data = model_data,
      unscaled_var = unscaled_var,
      scaled_var = scaled_var
    )
    
    # Plot
    plot_obj <- plot_marginal_effect_unscaled(
      marginal_df,
      xlab = unscaled_var,
      model_name = model_name
    )
    
    # Store it
    plot_key <- paste(model_name, unscaled_var, sep = "_")
    all_plots[[plot_key]] <- plot_obj
    
  }
}


# View one plot
bb <- all_plots[["diet/EM_Berries"]]
bb
bb + theme_classic()
remove(bb)


#OK. its's great that this is working. I'm sure i could do this more elegantly, but now
#pull out plots for everyboyd. i will just do so manually


EM_diet_Berries <- all_plots[["diet/EM_Berries"]] + diet_plot_style() + ggtitle("K. Infection - Berries") + ylim(0,1) + xlab("Berries (%)")
EM_diet_Nat_Prey <- all_plots[["diet/EM_Nat_Prey"]] + diet_plot_style() + ggtitle("L. Infection - Nat. Prey") + ylim(0,1) + xlab("Nat. Prey (%)")

CT_diet_Berries <- all_plots[["diet/CT_Berries"]] +diet_plot_style_2() + ggtitle("M. Shedding - Berries") + ylim(0,10) +xlab("Berries(%)")
CT_diet_Unk_Anthro <- all_plots[["diet/CT_Unk_Anthro"]] + diet_plot_style_2() + ggtitle("O. Shedding - Unk. Anthro.") + ylim(0,10) + xlab("Unk. Anthro. (%)")

blank_plot <- ggplot() + theme_void()

#DietPlots <- ggarrange(EM_diet_Berries, EM_diet_Nat_Prey,
 #                      CT_diet_Berries, CT_diet_Apples, CT_diet_Unk_Anthro, blank_plot,
  #                     nrow = 3, ncol = 2)
#ggsave("figures/DietMargEffects.png", DietPlots, width = 8, height = 12, dpi = 700,  bg = "white") 


EM_social_prophighincome <- all_plots[["social/EM_prophighincome"]] + social_plot_style() + ggtitle("P. Infection - High Income") + ylim(0,1) + xlab("Prop. High Income")
EM_social_Prop.Under15 <- all_plots[["social/EM_Prop.Under15"]] + social_plot_style() + ggtitle("Q. Infection - Children") + ylim(0,1) + xlab("Prop. Children (< 15)")

CT_social_prophighincome <- all_plots[["social/CT_prophighincome"]] + social_plot_style_2() + ggtitle("R. Shedding - High Income") + ylim(0,7) + xlab("Prop. High Income") 
CT_social_Prop.Over65 <- all_plots[["social/CT_Prop.Over65"]] + social_plot_style_2() + ggtitle("T. Shedding - Seniors") + ylim(0,7) + xlab("Prop. Seniors (> 65)")
CT_social_proplowincome <- all_plots[["social/CT_proplowincome"]] + social_plot_style_2() + ggtitle("S. Shedding - Low Income") + ylim(0,7) + xlab("Prop. Low Income")

#SocialPlots <- ggarrange(EM_social_prophighincome, EM_social_Prop.Under15,
 #                        CT_social_prophighincome, CT_social_proplowincome, CT_social_Prop.Over65, blank_plot,
  #                     nrow = 3, ncol = 2)
#ggsave("figures/SocialMargEffects.png", SocialPlots, width = 8, height = 9, dpi = 700,  bg = "white") 


#For dec terms, make sure you trasnform to m
EM_urb_cover_Anth_250m <- all_plots[["urb/EM_cover_Anth_250m"]] + urb_plot_style() + ggtitle("A. Infection - Anthro. Cover") + ylim(0,1) + xlab("Anthro. Land Cover (250m; %)")
EM_urb_decay_road_0005 <- all_plots[["urb/EM_decay_road_0005"]] + urb_plot_style() + ggtitle("E. Infection - Dist. to Road") + ylim(0,1) + xlab("Dist. to Road (m)") +
  scale_x_continuous(
  labels = function(decay_road_0005_scaled) {
    rounded_dist <- -log(1 - decay_road_0005_scaled) / 0.005
    round(rounded_dist)
  },
  limits = c(0, 0.95),
  breaks = scales::pretty_breaks(n = 5)
)
EM_urb_dist_to_Ravine_m <- all_plots[["urb/EM_dist_to_Ravine_m"]] + urb_plot_style() + ggtitle("D. Infection - Dist. to Ravine") + ylim(0,1) + xlab("Dist. to Ravine (m)")
EM_urb_decay_natarea_0005 <- all_plots[["urb/EM_decay_natarea_0005"]] + urb_plot_style() + ggtitle("C. Infection - Dist. to Nat. Area") + ylim(0,1) + xlab("Dist. to Nat. Area (m)") +
scale_x_continuous(
    labels = function(decay_natarea_0005_scaled) {
      rounded_dist <- -log(1 - decay_natarea_0005_scaled) / 0.005
      round(rounded_dist)
    },
    limits = c(0, 0.95),
    breaks = scales::pretty_breaks(n = 5)
  )

EM_urb_road_dens_250mH <- all_plots[["urb/EM_road_dens_250mH"]] + urb_plot_style() + ggtitle("B. Infection - Road Density") + ylim(0,1) + xlab("Road Density (m/ Ha)")

CT_urb_bldg_dens_25mH <- all_plots[["urb/CT_bldg_dens_25mH"]] + urb_plot_style_2() +ggtitle("G. Shedding - Building Density") + ylim(0,4) + xlab("Building Density (25-m)")
CT_urb_road_dens_250mH <- all_plots[["urb/CT_road_dens_250mH"]] + urb_plot_style_2() +ggtitle("G. Shedding - Road Density") + ylim(0,4) + xlab("Road Density (250-m)")
CT_urb_cover_Anth_50m <- all_plots[["urb/CT_cover_Anth_50m"]] + urb_plot_style_2() +  ggtitle("F. Shedding - Anthro. Land Cover") + ylim(0,7) + xlab("Anthro. Land Cover (50m; %)")
CT_urb_decay_RivEdge_00015 <- all_plots[["urb/CT_decay_RivEdge_00015"]] + urb_plot_style_2()+ggtitle("I. Shedding - Dist. to River") + ylim(0,7) + xlab("Dist. to River (m)") +
  scale_x_continuous(
    labels = function(decay_RivEdge_00015_scaled) {
      rounded_dist <- -log(1 - decay_RivEdge_00015_scaled) / 0.0015
      round(rounded_dist)
    },
    limits = c(0, 0.95),
    breaks = scales::pretty_breaks(n = 5)
  )

CT_urb_dist_to_NatArea_m <- all_plots[["urb/CT_dist_to_NatArea_m"]] + urb_plot_style_2() +ggtitle("H. Shedding - Dist. to Nat. Area") + ylim(0,7) + xlab("Dist to Nat. Area (m)")

CT_urb_decay_citybound_00015 <- all_plots[["urb/CT_decay_citybound_00015"]] + urb_plot_style_2() +ggtitle("J. Shedding - Dist. City Boundary") + ylim(0,7) + xlab("Dist to City Boundary (m)") +
  scale_x_continuous(
    labels = function(decay_citybound_00015_scaled) {
      rounded_dist <- -log(1 - decay_citybound_00015_scaled) / 0.0015
      round(rounded_dist)
    },
    limits = c(0, 0.95),
    breaks = scales::pretty_breaks(n = 5)
  )


#UrbPlots <- ggarrange(EM_urb_cover_Anth_250m, EM_urb_road_dens_250mH, EM_urb_decay_natarea_0005, EM_urb_dist_to_Ravine_m, EM_urb_decay_road_0005,
 #                     CT_urb_cover_Anth_50m, CT_urb_bldg_dens_250mH, CT_urb_dist_to_NatArea_m, CT_urb_decay_RivEdge_00015, CT_urb_decay_citybound_00015,
  #                       nrow = 5, ncol = 2)
#ggsave("figures/UrbMargEffects.png", UrbPlots, width = 8, height = 15, dpi = 700,  bg = "white") 


AllPlots <- ggarrange(EM_urb_cover_Anth_250m, EM_urb_road_dens_250mH, EM_urb_decay_natarea_0005, EM_urb_dist_to_Ravine_m, EM_urb_decay_road_0005,
                     CT_urb_cover_Anth_50m, CT_urb_bldg_dens_250mH, CT_urb_dist_to_NatArea_m, CT_urb_decay_RivEdge_00015, CT_urb_decay_citybound_00015,
                     EM_diet_Berries, EM_diet_Nat_Prey, CT_diet_Berries, CT_diet_Apples, CT_diet_Unk_Anthro,
                     EM_social_prophighincome, EM_social_Prop.Under15, CT_social_prophighincome, CT_social_proplowincome, CT_social_Prop.Over65,
                    nrow = 4, ncol = 5)
#ggsave("figures/Fig3.pdf", AllPlots, width = 15, height = 10, dpi = 700,  bg = "white") 



                
#Now majke a csv with the relevant info
# Helper to unscale a coefficient
unscale_beta <- function(beta, unscaled_var, data) {
  sd_val <- sd(data[[unscaled_var]], na.rm = TRUE)
  return(beta / sd_val)
}

# Master function to extract unscaled betas, effect sizes, and CIs
extract_effects <- function(mod_avg_obj, data, outcome_type = c("logistic", "gaussian")) {
  outcome_type <- match.arg(outcome_type)
  
  # Extract estimates
  coefs <- coef(mod_avg_obj, full = TRUE)
  coefs_df <- data.frame(term = names(coefs), estimate = coefs, row.names = NULL)
  
  # Get confidence intervals
  ci <- suppressWarnings(confint(mod_avg_obj, full = TRUE))
  ci_df <- data.frame(term = rownames(ci), conf.low = ci[, 1], conf.high = ci[, 2])
  
  # Join
  df <- dplyr::left_join(coefs_df, ci_df, by = "term")
  
  # Remove intercept
  df <- dplyr::filter(df, !grepl("Intercept", term, ignore.case = TRUE))
  
  # Identify term types
  df <- df %>%
    mutate(
      clean_term = gsub("cond\\((.*)\\)", "\\1", term),         # remove cond(...)
      is_quadratic = grepl("\\^2", clean_term),
      base_var = case_when(
        is_quadratic ~ gsub(".*\\((.*)_scaled\\^2\\)", "\\1", clean_term),
        TRUE ~ gsub("_scaled", "", clean_term)
      )
    )
  
  # Get SDs for all vars
  df <- df %>%
    mutate(
      sd_val = purrr::map_dbl(base_var, ~ if (.x %in% names(data)) sd(data[[.x]], na.rm = TRUE) else NA_real_)
    )
  
  # Unscale
  df <- df %>%
    mutate(
      beta_unscaled = ifelse(is_quadratic,
                             estimate / (sd_val^2),   # <-- key formula
                             estimate / sd_val),
      conf.low_unscaled = ifelse(is_quadratic,
                                 conf.low / (sd_val^2),
                                 conf.low / sd_val),
      conf.high_unscaled = ifelse(is_quadratic,
                                  conf.high / (sd_val^2),
                                  conf.high / sd_val),
      effect_size = case_when(
        outcome_type == "logistic" ~ exp(beta_unscaled),
        outcome_type == "gaussian" ~ beta_unscaled
      ),
      effect_low = case_when(
        outcome_type == "logistic" ~ exp(conf.low_unscaled),
        outcome_type == "gaussian" ~ conf.low_unscaled
      ),
      effect_high = case_when(
        outcome_type == "logistic" ~ exp(conf.high_unscaled),
        outcome_type == "gaussian" ~ conf.high_unscaled
      )
    )
  
  return(dplyr::select(df,
                       term = clean_term, base_var, is_quadratic,
                       beta_unscaled, effect_size, effect_low, effect_high
  ))
}


# Named list of models and their metadata
models_info <- list(
  diet_em = list(model = diet_em_avg, data = diet_em_dat, type = "logistic"),
  diet_ct = list(model = diet_ct_avg, data = diet_ct_dat, type = "gaussian"),
  social_em = list(model = social_em_avg, data = social_em_dat, type = "logistic"),
  social_ct = list(model = social_ct_avg, data = social_ct_dat, type = "gaussian"),
  urb_em = list(model = urb_em_avg, data = urb_em_dat, type = "logistic"),
  urb_ct = list(model = urb_ct_avg, data = urb_ct_dat, type = "gaussian")
)

# Loop through and collect output
results_list <- lapply(names(models_info), function(model_name) {
  mod_info <- models_info[[model_name]]
  df <- extract_effects(mod_info$model, mod_info$data, outcome_type = mod_info$type)
  df$model <- model_name
  df
})

# Combine and write to CSV
final_results <- bind_rows(results_list) %>%
  dplyr::select(model, term, beta_unscaled, effect_size, effect_low, effect_high)

final_results <- final_results %>%
  dplyr::filter(term != "(Int)")

#write.csv(final_results, "figures/unscaled_effects.csv", row.names = FALSE)


#OK I think the above wouldnt have worked for shedding intensity. but i like that code
#So here is the way to do it jsut for CT models
extract_gaussian_scaled_effects_avg <- function(model_avg, data) {
  
  # coefficients
  coefs <- coef(model_avg)
  
  # confidence intervals (matrix)
  ci <- confint(model_avg)
  
  # remove intercept
  intercept_idx <- grep("Intercept", names(coefs))
  if (length(intercept_idx) > 0) {
    coefs <- coefs[-intercept_idx]
    ci <- ci[-intercept_idx, , drop = FALSE]
  }
  
  # clean term names
  terms <- names(coefs)
  clean_terms <- gsub("cond\\((.*)\\)", "\\1", terms)
  
  # detect base variable name
  base_var <- gsub("_scaled.*", "", clean_terms)
  is_quad <- grepl("\\^2", clean_terms)
  
  # SD of original variables
  SDs <- sapply(base_var, function(v) sd(data[[v]], na.rm = TRUE))
  
  # unscale betas
  beta_scaled <- as.numeric(coefs)
  beta_unscaled <- ifelse(
    is_quad,
    beta_scaled / (SDs^2),
    beta_scaled / SDs
  )
  
  # CI unscaled
  ci_low <- ci[,1]
  ci_high <- ci[,2]
  ci_low_unscaled <- ifelse(is_quad, ci_low / (SDs^2), ci_low / SDs)
  ci_high_unscaled <- ifelse(is_quad, ci_high / (SDs^2), ci_high / SDs)
  
  # percent change
  percent_change <- 100 * (exp(beta_unscaled) - 1)
  percent_low <- 100 * (exp(ci_low_unscaled) - 1)
  percent_high <- 100 * (exp(ci_high_unscaled) - 1)
  
  data.frame(
    term = clean_terms,
    base_var = base_var,
    is_quadratic = is_quad,
    beta_unscaled = beta_unscaled,
    percent_change = percent_change,
    percent_low = percent_low,
    percent_high = percent_high
  )
}


gauss_res_urb <- extract_gaussian_scaled_effects_avg(urb_ct_avg, urb_ct_dat)
gauss_res_diet <- extract_gaussian_scaled_effects_avg(diet_ct_avg, diet_ct_dat)
gauss_res_social <- extract_gaussian_scaled_effects_avg(social_ct_avg, social_ct_dat)

gauss_res_urb
gauss_res_diet
gauss_res_social
    