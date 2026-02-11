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

model_averages <- readRDS("builds/model_avg_objects_cluster/model_averages.rds")
model_list_nested <- readRDS("builds/model_avg_objects_cluster/model_list_nested.rds")

#Read in model averaged objects
diet_em_avg <- readRDS("builds/model_avg_objects_cluster/diet_em.rds")
social_em_avg <- readRDS("builds/model_avg_objects_cluster/social_em.rds")
social_ct_avg <- readRDS("builds/model_avg_objects_cluster/social_ct.rds")
urb_em_avg <- readRDS("builds/model_avg_objects_cluster/urb_em.rds")
urb_ct_avg <- readRDS("builds/model_avg_objects_cluster/urb_ct.rds")

#read in the data 
diet_em_dat <- readRDS("builds/dat_for_modavg_cluster/diet_EM.rds")
social_em_dat <- readRDS("builds/dat_for_modavg_cluster/social_EM.rds")
social_ct_dat <- readRDS("builds/dat_for_modavg_cluster/social_CT.rds")
urb_em_dat <- readRDS("builds/dat_for_modavg_cluster/urb_EM.rds")
urb_ct_dat <- readRDS("builds/dat_for_modavg_cluster/urb_CT.rds")



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

#let's see if we can loop this

#Sytart by kind of mapping everything in a list
# A named list mapping model identifiers to their model list, dataset, and predictors
model_info <- list(
  "diet/EM" = list(
    model_list = model_list_nested[["diet/EM"]],
    data = diet_em_dat,
    vars = list(propAnthro = "propAnthro_scaled",
                propSeed = "propSeed_scaled",
                compost = "compost_scaled",
                propApple = "propApple_scaled",
                propGarbage = "propGarbage_scaled")
  ),
  
  "social/EM" = list(
    model_list = model_list_nested[["social/EM"]],
    data = social_em_dat,
    vars = list(lowincome = "lowincome_scaled", 
                Over65 = "Over65_scaled")
  ),
  "social/CT" = list(
    model_list = model_list_nested[["social/CT"]],
    data = social_ct_dat,
    vars = list(highincome = "highincome_scaled",
                NoCertDiplomadegree = "NoCertDiplomadegree_scaled")
  ),
  "urb/EM" = list(
    model_list = model_list_nested[["urb/EM"]],
    data = urb_em_dat,
    vars = list(decay_road_0005 = "decay_road_0005_scaled",
                decay_citybound_00015 = "decay_citybound_00015_scaled")
  ),
  "urb/CT" = list(
    model_list = model_list_nested[["urb/CT"]],
    data = urb_ct_dat,
    vars = list(dist_to_RivEdge_m = "dist_to_RivEdge_m_scaled",
                rd_density = "rd_density_scaled")
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
bb <- all_plots[["diet/EM_propSeed"]]
bb
bb + theme_classic()
remove(bb)


#OK. its's great that this is working. I'm sure i could do this more elegantly, but now
#pull out plots for everyboyd. i will just do so manually
EM_diet_Anthro <- all_plots[["diet/EM_propAnthro"]] + diet_plot_style() + ggtitle("F. Infection - Unk. anthro") + ylim(0,0.9) + xlab("Unk. anthro (prop.)")
EM_diet_Seed <- all_plots[["diet/EM_propSeed"]] + diet_plot_style() + ggtitle("G. Infection - Birdseed") +  ylim(0,0.9) +  xlab("Birdseed (prop.)")
EM_diet_Garbage <- all_plots[["diet/EM_propGarbage"]] + diet_plot_style() + ggtitle("I. Infection - Garbage") + ylim(0,0.9) + xlab("Garbage (prop.)")
EM_diet_Compost <- all_plots[["diet/EM_compost"]] + diet_plot_style() + ggtitle("E. Infection - Compost") +  ylim(0,0.9) + xlab("Compost prevalence")
EM_diet_Apples <- all_plots[["diet/EM_propApple"]] + diet_plot_style() + ggtitle("H. Infection - Apples") + ylim(0,0.9) + xlab("Apples (prop.)")

blank_plot <- ggplot() + theme_void()

EM_social_lowincome <- all_plots[["social/EM_lowincome"]] + social_plot_style() + ggtitle("J. Infection - Low income") + xlab("Prop. Low income") + ylim(0,1)
EM_social_Over65 <- all_plots[["social/EM_Over65"]] + social_plot_style() + ggtitle("K. Infection - Seniors") + xlab("Prop. seniors (> 65)") + ylim(0,1)

CT_social_highincome <- all_plots[["social/CT_highincome"]] + social_plot_style_2() + ggtitle("M. Shedding - High income") + xlab("Prop. High income") +ylim(0,6)
CT_social_Education <- all_plots[["social/CT_NoCertDiplomadegree"]] + social_plot_style_2() + ggtitle("L. Shedding - No post secondary") + xlab("Prop. no post secondary") +ylim(0,6)


EM_urb_decay_road_0005 <- all_plots[["urb/EM_decay_road_0005"]] + urb_plot_style() + ggtitle("A. Infection -\nDist. to road") + xlab("Dist. to road (m)") + ylim(0,1) +
  scale_x_continuous(
  labels = function(decay_road_0005_scaled) {
    rounded_dist <- -log(1 - decay_road_0005_scaled) / 0.005
    round(rounded_dist)
  },
  limits = c(0, 0.95),
  breaks = scales::pretty_breaks(n = 5)
)

EM_urb_decay_citybound_00015 <- all_plots[["urb/EM_decay_citybound_00015"]] + urb_plot_style() + ggtitle("B. Infection -\nDist. to city boundary") +  xlab("Dist. to city boundary (m)") + ylim(0,1) +
scale_x_continuous(
    labels = function(decay_citybound_00015_scaled) {
      rounded_dist <- -log(1 - decay_citybound_00015_scaled) / 0.0015
      round(rounded_dist)
    },
    limits = c(0.965, 0.99),
  #  breaks = scales::pretty_breaks(n = 5)
  )

CT_urb_rd_density <- all_plots[["urb/CT_rd_density"]] + urb_plot_style_2() +ggtitle("C. Shedding - Road density") + xlab("Road density (m/Ha)") + ylim(0,6)
CT_urb_dist_rivedge <- all_plots[["urb/CT_dist_to_RivEdge_m"]] + urb_plot_style_2() +  ggtitle("D. Shedding - Dist. to river")  + xlab("Dist. to river (m)") + ylim(0,6)

Row1 <- ggarrange(EM_urb_decay_road_0005, EM_urb_decay_citybound_00015, CT_urb_rd_density, CT_urb_dist_rivedge, nrow = 1)
Row2 <- ggarrange(EM_diet_Compost, EM_diet_Anthro, EM_diet_Seed, EM_diet_Apples, EM_diet_Garbage, nrow = 1)
Row3 <- ggarrange(EM_social_lowincome, EM_social_Over65, CT_social_Education, CT_social_highincome, nrow = 1)

AllRows <- ggarrange(Row1, Row2, Row3, nrow = 3)

#ggsave("figures/Fig4_Nov26.pdf", AllRows, width = 14, height = 8, dpi = 700,  bg = "white") 



                
#Now majke a csv with the relevant info
# Helper to unscale a coefficient
unscale_beta <- function(beta, unscaled_var, data) {
  sd_val <- sd(data[[unscaled_var]], na.rm = TRUE)
  return(beta / sd_val)
}

# Master function to extract unscaled betas, effect sizes, and CIs
extract_effects <- function(mod_avg_obj, data, outcome_type = c("binomial", "gaussian")) {
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
        outcome_type == "binomial" ~ exp(beta_unscaled),
        outcome_type == "gaussian" ~ beta_unscaled
      ),
      effect_low = case_when(
        outcome_type == "binomial" ~ exp(conf.low_unscaled),
        outcome_type == "gaussian" ~ conf.low_unscaled
      ),
      effect_high = case_when(
        outcome_type == "binomial" ~ exp(conf.high_unscaled),
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
  diet_em = list(model = diet_em_avg, data = diet_em_dat, type = "binomial"),
  social_em = list(model = social_em_avg, data = social_em_dat, type = "binomial"),
  social_ct = list(model = social_ct_avg, data = social_ct_dat, type = "gaussian"),
  urb_em = list(model = urb_em_avg, data = urb_em_dat, type = "binomial"),
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

#write.csv(final_results, "figures/unscaled_effects_clusters_Nov26.csv", row.names = FALSE)
