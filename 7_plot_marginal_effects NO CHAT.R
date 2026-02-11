
rm(list = ls())
gc()

# Load data and prepare workspace -----------------------------------------
library(data.table)
library(ggplot2)
library(ggeffects)
library(stats)
library(ggpubr)
library(broom)
library(stringr)
library(MuMIn)
library(glmmTMB)

model_averages <- readRDS("builds/model_avg_objects/model_averages.rds")

names(model_averages)
model_averages[["diet/CT"]] |> summary() #correct predictors
model_averages[["diet/EM"]] |> summary() # correct predictors
model_averages[["social/CT"]] |> summary() #correct preds
model_averages[["social/EM"]] |> summary() #correct preds
model_averages[["urb/CT"]] |> summary() #correct preds
model_averages[["urb/EM"]] |> summary() #correct preds

#read in the data 
files <- list.files("builds/New_data_for_modavg", pattern = ".rds",
                    full.names = T)
files
dat.list <- lapply(files, readRDS)
names(dat.list) <- list.files("builds/New_data_for_modavg", pattern = ".rds",
                              full.names = F)

# The names of the data MUSt match the names of the models in the model list
names(model_averages)
names(dat.list) <- gsub(".rds", "", names(dat.list))
names(dat.list) <- word(names(dat.list), 2, 3, sep = "_")
names(dat.list) <- gsub("_", "/", names(dat.list))
names(dat.list)
setdiff(names(model_averages), names(dat.list))
setdiff(names(dat.list), names(model_averages))
#' [THESE MUST BE LENGTH 0]


# Write a helper function
create_all_grids <- function(m){
  # This creates predictions for each variable, controlling all others at the mean
  # it assumes all variables are continuous.
  dat <- m$x |> as.data.frame()
  setDT(dat)
  
  dat <- dat[, !c("(Intercept)"), with = F]
  
  # Deal with I()^2 where present
  dat <- dat[, !grepl("I\\(", names(dat)), with = F]
  
  # loop through columns and calculate mean of each
  mean_list <- lapply(dat, function(x){
    mean(x)
  })
  
  # Now cycle through, sequencing along one variable at a time and controlling for the mean of the others
  grid <- list()
  for(k in 1:ncol(dat)){
    means <- mean_list[-k] |> unlist()
    means <- data.table(t(means))
    names(means) <- names(mean_list)[-k]
    
    grid[[k]] <- data.table(X = seq(min(dat[, ..k]), max(dat[, ..k]), 
                                    length.out = 100 ),
                            variable_predicted = names(dat)[k],
                            means)
    setnames(grid[[k]], "X", names(dat)[k])
  }
  grid <- rbindlist(grid, use.names = TRUE)
  return(grid)
}

unscale_all_variables <- function(original_dat, new_dat){
  # original data is the unscaled data
  # new_dat is the scaled data we want to unscale
  # OK this is a pain in the fucking ass
  new <- copy(new_dat)
  original <- copy(original_dat)
  
  # This all assumes that "_scaled" designates scaled columns.
  scaled_vars <- names(new)[names(new) %in% names(original)]
  unscaled_vars <- gsub("_scaled", "", scaled_vars)
  if(!all(unscaled_vars %in% names(original_dat))) print("FUCK YOU LUNDY")
  
  moments <- lapply(unscaled_vars,
         function(x){
           # print(x)
           y <- original[, x] |> unlist()
           center <- mean(y)
           sd <- sd(y)
           return(data.frame(center, sd))
         })
  names(moments) <- scaled_vars
  moments <- rbindlist(moments, idcol = "scaled_vars")
  
  new.unscaled <- lapply(1:nrow(moments),
                             function(j){
                               
                               nm <- moments[j, ]$scaled_vars
                               y <- unlist(new[, nm, with = F]) * moments[j, ]$sd + moments[j, ]$center
                               y <- unname(y)
                               y
                               out <- data.table(Y = y)
                               
                               setnames(out, "Y", gsub("_scaled", "", nm))
                               return(out) 
                })
  new <- data.table(new,
                        do.call(cbind, new.unscaled)) #data.table::cbindlist(new.unscaled)
  return(new)
  
}



# Predictions -------------------------------------------------------------


# Loop through, extract predictions, and then unscale.
m <- c()
dat <- c()
grid <- c()
out <- c()
pred <- list()
i <- 1

for(i in 1:length(model_averages)){
  
  m <- model_averages[[i]]
  names(model_averages)[i]
  names(dat.list)
  dat <- dat.list[[names(model_averages)[i]]]

  grid <- create_all_grids(m)
  out <- predict(m, newdata = grid[, !c("variable_predicted")],
                        type = "response", se.fit = TRUE, re.form = NA)
  pred[[i]] <- data.table(grid,
                     pred = out$fit,
                     se = out$se)
  
  # Now unscale
  pred[[i]] <- unscale_all_variables(original_dat = dat,
                                     new_dat = pred[[i]])
  pred[[i]][, variable_predicted := gsub("_scaled", "", variable_predicted)]

}

names(pred) <- names(model_averages)
pred

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ---------------------------------------
# Plot --------------------------------------------------------------------

pred
# CT is shedding intensity, EM is infection status.

# First let's assign aesthetics to columns within the dataset.
# linetype:
# CT = "dashed"
# EM = "solid"

# >>> Prepare aesthetics --------------------------------------------------
nms <- names(pred)
x <- 1
pred <- lapply(1:length(nms), function(x){
  pred[[x]][, response := ifelse(grepl("EM",nms[x]),
                                 "Infection status", "Shedding intensity (log)")]
  pred[[x]][, predictor := fcase(grepl("diet", nms[x]), "Host diet",
                                 grepl("social", nms[x]), "Neighbourhood demography",
                                 grepl("urb", nms[x]), "Urbanization",
                                 default = NA)]
  pred[[x]][, line_type := ifelse(response == "Infection status", "solid", "dashed")]
  pred[[x]][, fill_type := fcase(predictor == "Host diet", "#50C5B7",
                                 predictor == "Neighbourhood demography", "#9CEC5B",
                                 predictor == "Urbanization", "#6184D8",
                                 default = NA)]
  pred[[x]][, line_color := fcase(predictor == "Host diet", "#50C5B7",
                                 predictor == "Neighbourhood demography", "#9CEC5B",
                                 predictor == "Urbanization", "#6184D8",
                                 default = NA)]
  return(pred[[x]])
})
pred

lapply(pred, function(x) x$line_type)
lapply(pred, function(x) x$response)

raw_vars <- lapply(1:length(pred), function(x) unique(pred[[x]]$variable_predicted)) |> 
  unlist()

var_key <- data.table(variable_predicted = unique(raw_vars))
var_key[, variable_label := fcase(variable_predicted == "Berries", "Berries (%)",
                                  variable_predicted == "Nat_Prey", "Nat. prey (%)",
                                  variable_predicted == "Unk_Anthro", "Unk. anthro (%)",
                                  variable_predicted == "prophighincome", "Prop high income",
                                  variable_predicted == "Prop.Under15", "Prop. children",
                                  variable_predicted == "Prop.Over65", "Prop. seniors",
                                  variable_predicted == "proplowincome", "Prop. low income",
                                  variable_predicted == "cover_Anth_250m", "Anthro cover (%)",
                                  variable_predicted == "decay_road_0005", "Dist. to road (m)",
                                  variable_predicted == "dist_to_Ravine_m", "Dist. to ravine (m)",
                                  variable_predicted == "decay_natarea_0005", "Dist. to nat. area (m)",
                                  variable_predicted == "road_dens_250mH", "Road density",
                                  variable_predicted == "dist_to_NatArea_m", "Dist. to nat. area (m)",
                                  variable_predicted == "bldg_dens_25mH", "Building density",
                                  variable_predicted == "decay_Ravine_002", "Dist. to ravine (m)")]

pred <- lapply(pred, function(x){
  merge(x, var_key, all.x = T, all.y = F, by = "variable_predicted")
})
pred

# Now melt so that predictors are LONG
# There are some errors, I think some of the dataset names are scrambled because variables are missing
#
pred.mlt <- lapply(1:length(pred), function(x){
  x.mlt <- melt(pred[[x]][, !grepl("_scaled", names(pred[[x]])), with = F],
       measure.vars = unique(pred[[x]]$variable_predicted),
       value.name = "predictor_value",
       variable.name = "predictor_variable")
  x.mlt <- x.mlt[variable_predicted == predictor_variable]
  x.mlt$variable_predicted <- NULL
  return(x.mlt)
})


lapply(pred.mlt, function(x) unique(x$predictor_variable )) |> unlist() |> unique()

# Transform decay values
pred.mlt <- lapply(1:length(pred.mlt),
               function(x){
                 pred.mlt[[x]][, decay_predictor := ifelse(grepl("decay", predictor_variable),
                                                       "yes", "no")]
                 
                 pred.mlt[[x]][decay_predictor == "yes", decay_alpha_char := word(predictor_variable, -1, sep = "_")]
                 pred.mlt[[x]][decay_predictor == "yes", decay_alpha_char := paste0(substr(decay_alpha_char, 1, 1),
                                                       ".",
                                                       substr(decay_alpha_char, 2, nchar(decay_alpha_char) ) )]
                 pred.mlt[[x]][, decay_alpha := as.numeric(decay_alpha_char)]
                 #
                 pred.mlt[[x]][, trans_predictor_value := ifelse(decay_predictor == "yes",
                                                                   -log(1 - predictor_value) / decay_alpha, 
                                                                   predictor_value)]
                 # pred.mlt[[x]]$decay_alpha_char <- NULL
               })


test <- "12345"
paste0(substr(test, 1, 1), ".", substr(test, 2, nchar(test)))

pred.mlt[[5]]

# >>> Plot ----------------------------------------------------------------
dt <- pred.mlt[[1]]
plot_model <- function(dt, n_row = 1){
  ggplot(data = dt) +
    geom_ribbon(aes(x = trans_predictor_value, ymin = pred - se, ymax = pred + se,
                    fill = fill_type), alpha = .5) +
    geom_path(aes(x = trans_predictor_value, y = pred, color = line_color,
                  linetype = line_type), size = 2) +
    scale_fill_identity() +
    scale_linetype_identity() +
    scale_x_continuous(limits = c(0, NA)) + 
    scale_color_identity() +
    facet_wrap(~variable_label, scales = "free_x", nrow = n_row) +
    xlab("Predictor value") +
    ylab(unique(dt$response)) +
    theme_classic() +
    theme(plot.title = element_text(size = 12, face = "plain"),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 12)) +
    theme(strip.background = element_blank(),
          strip.text = element_text(size = 12))
}


figs <- lapply(pred.mlt, plot_model, n_row = 1)
figs
figs[1] #em + host diet
figs[2] #ct + host dite
figs[3] #em + social
figs[4] #ct + social
figs[5] #em + urb
figs[6] #ct + urb

blank_plot <- ggplot() + theme_void()

diet <- ggarrange(figs[[1]], figs[[2]], blank_plot, nrow = 1, widths = c(2,2,1))
social <- ggarrange(figs[[3]], figs[[4]], nrow = 1, widths = c(2,3))
urb1 <- ggarrange(figs[[6]], blank_plot, nrow = 1, widths = c(4,1))

finalplot <-ggarrange(figs[[5]], urb1, diet, social, nrow = 4)
#ggsave("figures/Fig3_Nov23.pdf", finalplot, width = 15, height = 10, dpi = 700,  bg = "white") 

