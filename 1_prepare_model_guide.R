#' *Prepare model guide for univariates only!*
#
#
#
# Prepare workspace ---------------------------------------
#
#
#

rm(list = ls())
gc()

#Load libraries
library(data.table)


# ~~~~~~~~~~~~~~~~~~~~~~~~ ------------------------------------------------
#  Create model guide for univariate models
# Do so for both outcome variables

#0.Read in data ---------------------------------------------------------------
dat <- readRDS("builds/prepared_data.rds")

#1. Prep guides------------------------------------------------------------------
guide <- CJ(response = c("Emulti", "CT"),
            var = c("Apples", "Berries", "Birdseed", "Vegetation", "Nat_Prey",
                    "Garbage", "Unk_Anthro",
                    "proplowincome", "prophighincome",
                    "Prop.NoCertDiplomadegree", 
                    "Prop.Over65", "Prop.Under15", "Dogs",
                    "dist_to_road_m", "dist_to_Bldg_m", "dist_to_CityCenter_m",
                    "dist_to_CityBounary_m", "dist_to_NatArea_m", "dist_to_Ravine_m",
                    "cover_Nat_25m", "cover_Anth_25m", 
                    "cover_Nat_50m", "cover_Anth_50m", 
                    "cover_Nat_250m", "cover_Anth_250m",
                    "bldg_dens_25mH", "bldg_dens_50mH", "bldg_dens_250mH",
                    "road_dens_25mH", "road_dens_50mH", "road_dens_250mH",
                    "urban_index_pca",
                    "decay_road_002", "decay_road_0005", "decay_road_00015",
                    "decay_bldg_002", "decay_bldg_0005", "decay_bldg_00015",
                    "decay_citycenter_002", "decay_citycenter_0005", "decay_citycenter_00015",
                    "decay_citybound_002", "decay_citybound_0005", "decay_citybound_00015",
                    "decay_Ravine_002", "decay_Ravine_0005", "decay_Ravine_00015",
                    "decay_RivEdge_002", "decay_RivEdge_0005", "decay_RivEdge_00015",
                    "decay_natarea_002", "decay_natarea_0005", "decay_natarea_00015",
                    "year"))


head(guide)

#
setdiff(guide$var, names(dat)) #this tells me if anything in guide is missing from dat
#REMEMBER that things will still need to be scaled at some point. On the fly


# >>> Create formulas -----------------------------------------------------

#Create null model formulas
guide[, null_model_formula := paste(response, "~ 1")]

#Create univariate formulas
guide[, univariate_formula := paste(response, "~", var)]
head(guide)

# Specify the data you need to use for each model
guide[, data := ifelse(response == "CT", "pos_scats", "all_scats")]
guide

# >>> Specify model family -----------------------------------------
unique(guide$response)

guide[response == "CT", model_family := "gaussian()"]
guide[response == "Emulti", model_family := "binomial(link = 'logit')"]

guide[is.na(model_family), ] #everybody has a family. :)

#Create quadratic models
numeric_vars <- names(dat)[sapply(dat, is.numeric)] #get names of numeroic dudes

guide2 <- guide[var %in% numeric_vars]
guide2$quad <- "quadratic"
guide$quad <- "linear"

#USE SCARY REGEX
guide2[, univariate_formula := gsub(
  pattern = "^(\\s*[^~]+\\s*~\\s*)(\\w+)$",
  replacement = "\\1\\2 + I(\\2^2)",
  x = univariate_formula
)]

#merge guides back together
guides_merged <- rbind(guide, guide2)
guide <- guides_merged

# >>> Add model comparison IDs --------------------------------------------

guide[, model_complexity_comparison_ID := paste0("model_complexity_id_",
                                                 seq(1:.N))]

guide[model_complexity_comparison_ID == "model_complexity_id_96"]

# >>> Add an exclusion formula to make sure models are comparable --------------------------------------------------
#the goal here is to make a column that IDs the models that should be compared

guide[, exclusion := paste0("complete.cases(",
                            var, 
                            ", ", response, ")")]

unique(guide$exclusion)


# >>> Make the data long --------------------------------------------------

guide.long <- melt(guide,
                   measure.vars = c("null_model_formula",
                                    "univariate_formula"),
                   variable.name = "model_type",
                   value.name = "formula")

guide.long

# Drop the NA models:
guide.long <- guide.long[!is.na(formula), ]
guide.long

guide.long[, model_type := gsub("_formula", "", model_type)]

# >>> Add an individual model ID ------------------------------------------

guide.long[, model_id := paste0("model_", seq(1:.N))]
guide.long

guide.long[, model_path := paste0("builds/batchmods_nov/models/", model_id, ".Rds")]
guide.long

guide.long[duplicated(model_path), ] #No duplicates! yay


# >>> Add random effects to formulas --------------------------------------
guide.long[,formula := paste0(formula, " + (1|ClusterID)")]
unique(guide.long$formula)

# >>> Test that formulas are correctly formed -----------------------------

for(i in 1:nrow(guide.long)){
  as.formula(guide.long[i, ]$formula)
} # This will errror out if there's a syntactical mistake

x <- c()
for(i in 1:nrow(guide.long)){
  x <- dat[eval(parse(text = guide.long[i, ]$exclusion)), ]
  if(nrow(x) == 0){
    print("Watch it buddy")
  }
}



# >>> Create ELEGANT executable call in guide -----------------------------
guide.long[, model_call := paste0("glmmTMB(", 
                                  formula, ", ",
                                  "family=", model_family, ", ",
                                  "data = sub_dat)")]

guide.long[, model_call := paste0("glmmTMB(", 
                                    formula, ", ",
                                    "family=", model_family, ", ",
                                    "data = sub_dat)")]
guide.long

guide.long[1, ]$model_call
# eval(parse(text = guide.long[1, ]$model_call))

#


# ~~~~~~~~~~~~~~~~~~~~~~~~ ------------------------------------------------
#  Save guide ------------------------------------
saveRDS(guide.long, "builds/batchmods_nov/model_guide.Rds")

