
rm(list = ls())
gc()

library("data.table")
library("glmmTMB")
library("doParallel")
library("foreach")
library("doSNOW")
library("crayon")
library("ggplot2")
library(tidyverse)

# groundhog.library(libs, groundhog.day)

# Load data ------------------------------------------------------------
#dat <- read_csv("data/ClustDat_Oct23.csv")

#dat <- dat %>%
#  dplyr::mutate(bldg_density = n_buildings / area_m2,
#                rd_density = road_length_m / area_m2,
#                Prop_Nat = Nat / area_m2,
#                Prop_Anth = Anthro / area_m2)
#dat$ClusterID <- as.factor(dat$ClusterID)

#dat2 <- readRDS("builds/prepared_data.rds")

#dat2 <- dat2 %>% group_by(ClusterID, Emulti) %>% summarise(n = n())
#dat2 <- dat2 %>%
#  dplyr::filter(Emulti == 1)
#dat2 <- dat2 %>%
#  dplyr::select(-(Emulti)) %>%
#  dplyr::rename("No.Pos" = n)
#dat2 <- as.data.frame(dat2)

#dat_final <- left_join(dat, dat2, by = "ClusterID")
#dat_final <- as.data.frame(dat_final)

#saveRDS(dat_final, "builds/prepared_cluster_data.rds")


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


#target model looks like this:
#model <- glm(cbind(infected, total - infected) ~ predictors, #infected: no infected scats in cluster
#total: total scats in the cluster
             #data = data,
             #family = binomial)



# Make model guide
#1. Prep guides------------------------------------------------------------------
guide <- CJ(response = c("cbind(No.Pos, scatcount - No.Pos)", "meanSI"),
            var = c("propApple", "propBerries", "propSeed", "propVeg", "propPrey",
                    "propGarbage", "propAnthro",
                    "lowincome", "highincome",
                    "NoCertDiplomadegree", 
                    "Over65", "Under15", "Dogs1",
                    "dist_to_road_m", "dist_to_Bldg_m", "dist_to_CityCenter_m",
                    "dist_to_CityBounary_m", "dist_to_NatArea_m", "dist_to_Ravine_m",
                    "dist_to_RivEdge_m",
                    "Prop_Nat", "Prop_Anth",
                    "bldg_density", "rd_density",
                    "decay_road_002", "decay_road_0005", "decay_road_00015",
                    "decay_bldg_002", "decay_bldg_0005", "decay_bldg_00015",
                    "decay_citycenter_002", "decay_citycenter_0005", "decay_citycenter_00015",
                    "decay_citybound_002", "decay_citybound_0005", "decay_citybound_00015",
                    "decay_Ravine_002", "decay_Ravine_0005", "decay_Ravine_00015",
                    "decay_RivEdge_002", "decay_RivEdge_0005", "decay_RivEdge_00015",
                    "decay_natarea_002", "decay_natarea_0005", "decay_natarea_00015",
                    "compost", "garden"
                    ))


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

# >>> Specify model family -----------------------------------------
unique(guide$response)

guide[response == "meanSI", model_family := "gaussian()"]
guide[response == "cbind(No.Pos, scatcount - No.Pos)", model_family := "binomial"]

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

guide.long[, model_path := paste0("builds/batchmods_nov25/models/", model_id, ".Rds")]
guide.long

guide.long[duplicated(model_path), ] #No duplicates! yay


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
                                  "family =", model_family, ", ",
                                  "data = sub_dat")]

guide.long

guide.long[1, ]$model_call



# eval(parse(text = guide.long[1, ]$model_call))

#


# ~~~~~~~~~~~~~~~~~~~~~~~~ ------------------------------------------------
#  Save guide ------------------------------------
saveRDS(guide.long, "builds/batchmods_nov25/model_guide.Rds")





