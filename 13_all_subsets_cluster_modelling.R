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

# Load data and prepare workspace -----------------------------------------
library(data.table)
library(dplyr)
library(MuMIn)
library(glmmTMB)

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


#Start with DIET + Emulti
#From3 .
dat_diet_EM <- dat %>%
  dplyr::select(No.Pos, scatcount, propAnthro, propApple, propGarbage, propSeed, compost) %>%
  dplyr::mutate(propAnthro_scaled = scale(propAnthro),
                propApple_scaled = scale(propApple),
                propGarbage_scaled = scale(propGarbage),
                propSeed_scaled = scale(propSeed),
                compost_scaled = scale(compost))

dat_diet_EM <- na.omit(dat_diet_EM) #n = 27

#Check for corelated dudes
numeric_cols <- sapply(dat_diet_EM, is.numeric)
cor(dat_diet_EM[, numeric_cols, with = FALSE], method = "pearson", use = "pairwise.complete.obs")

#No correlation concerns

# fit model with all parameters
diet_EM<-glmmTMB(cbind(No.Pos, scatcount - No.Pos) ~ propAnthro_scaled + propApple_scaled + propSeed_scaled +
                   propGarbage_scaled + compost_scaled,
                data = dat_diet_EM, 
                family = binomial,
                na.action = "na.fail")

#All subsets
diet_EM_results <- dredge(diet_EM, 
                          trace = 2)

diet_EM_results #4 top models; first four

#From3 . Shedding intensity: apples Q, berries L, birdseed Q, Unk Anthro L:
#Nothing for shedding intensity

#OK next up is social stuff

#From3 Infection: prop children (l), prop high income (q), prop low income (l),
dat_social_EM <- dat %>%
  dplyr::select(No.Pos, scatcount, lowincome, NoCertDiplomadegree, Over65) %>%
  dplyr::mutate(Over65_scaled = scale(Over65),
                NoCertDiplomadegree_scaled = scale(NoCertDiplomadegree),
                lowincome_scaled = scale(lowincome))

dat_social_EM <- na.omit(dat_social_EM) #n = 18

#Check for corelated dudes
numeric_cols <- sapply(dat_social_EM, is.numeric)
cor(dat_social_EM[, numeric_cols, with = FALSE], method = "pearson", use = "pairwise.complete.obs")

# no cert is correlated with everything

# fit model with all parameters
social_EM<-glmmTMB(cbind(No.Pos, scatcount - No.Pos) ~ Over65_scaled + NoCertDiplomadegree_scaled + 
                     lowincome_scaled +I(NoCertDiplomadegree_scaled^2) + I(Over65_scaled^2),
                   data = dat_social_EM, 
                   family = binomial,
                   na.action = "na.fail")


social_EM_results <- dredge(
  social_EM,
  trace = 2)
social_EM_results #2 good ones, 9 and 10

#From3 . Shedding intensity: no post secondary (L), seniors (L), children (L), high income (q), low income (q)
dat_social_CT <- dat %>%
  dplyr::select(meanSI, NoCertDiplomadegree, highincome) %>%
  dplyr::mutate(NoCertDiplomadegree_scaled = scale(NoCertDiplomadegree),
                highincome_scaled = scale(highincome))

dat_social_CT <- na.omit(dat_social_CT) #n = 16

#Check for corelated dudes
numeric_cols <- sapply(dat_social_CT, is.numeric)
cor(dat_social_CT[, numeric_cols, with = FALSE], method = "pearson", use = "pairwise.complete.obs")
#they are correlated


# fit model with all parameters
social_CT<-glmmTMB(meanSI ~ NoCertDiplomadegree_scaled + highincome_scaled,
                 data = dat_social_CT, 
                 family = gaussian(), 
                 na.action = "na.fail")

#All subsets
social_CT_results <- dredge(social_CT, 
                          trace = 2)
social_CT_results #models 1 amd 3 are good


#Finally urbanization
# anth cover (quad), city boundary distance (l), Nat. Area Dec.002 (l), nat cover (q)
#Riv. Dec.0005 (quad), Road Dec.002 (l), road density (l)
dat_urb_EM <- dat %>%
  dplyr::select(No.Pos, scatcount, Anthro, decay_bldg_002,
                decay_citybound_00015, decay_citycenter_00015,
                decay_road_0005, rd_density) %>%
  dplyr::mutate(Anthro_scaled = scale(Anthro),
                decay_bldg_002_scaled = scale(decay_bldg_002),
                decay_citybound_00015_scaled = scale(decay_citybound_00015),
                decay_citycenter_00015_scaled = scale(decay_citycenter_00015),
                decay_road_0005_scaled = scale(decay_road_0005),
                rd_density_scaled = scale(rd_density))

dat_urb_EM <- na.omit(dat_urb_EM) #n = 23

#Check for corelated dudes
numeric_cols <- sapply(dat_urb_EM, is.numeric)
#dec bldg and dec rd


# fit model with all parameters
urb_EM<-glmmTMB(cbind(No.Pos, scatcount - No.Pos) ~ Anthro_scaled  + decay_bldg_002_scaled + decay_citybound_00015_scaled +
                  decay_citycenter_00015_scaled + decay_road_0005_scaled + rd_density_scaled +
                  I(Anthro_scaled^2) + I(decay_bldg_002_scaled^2) + I(decay_citybound_00015_scaled^2) +
                  I(decay_road_0005_scaled^2) + I(rd_density_scaled^2),
                data = dat_urb_EM, 
                family = binomial,
                na.action = "na.fail")


urb_EM_results <- dredge(urb_EM,
  trace = 2,
  m.lim = c(1, 2))

urb_EM_results #top 1 is good, nbopdy else

#From3 . Shedding intensity: 
# all lin
# cc dist, riv dist, road density, Ravine Dec.002
dat_urb_CT <- dat %>%
  dplyr::select(meanSI, rd_density, dist_to_CityCenter_m, dist_to_RivEdge_m,
                decay_Ravine_002) %>%
  dplyr::mutate(rd_density_scaled = scale(rd_density),
                dist_to_CityCenter_m_scaled = scale(dist_to_CityCenter_m),
                dist_to_RivEdge_m_scaled = scale(dist_to_RivEdge_m),
                decay_Ravine_002_scaled = scale(decay_Ravine_002))


dat_urb_CT <- na.omit(dat_urb_CT) #n = 19

#Check for corelated dudes
numeric_cols <- sapply(dat_urb_CT, is.numeric)
cor(dat_urb_CT[, numeric_cols, with = FALSE], method = "pearson", use = "pairwise.complete.obs")

#riv and rav are correlated

# fit model with all parameters
urb_CT<-glmmTMB(meanSI ~ rd_density_scaled + dist_to_CityCenter_m_scaled + 
                  dist_to_RivEdge_m_scaled  + decay_Ravine_002_scaled,
                   data = dat_urb_CT, 
                   family = gaussian(), 
                   na.action = "na.fail")
#All subsets
urb_CT_results <- dredge(urb_CT, 
                            trace = 2,
                         m.lim = c(1,2))

urb_CT_results #top 2 are good

  
#Now I want to save within 2AIC models in Rds manner so they are easy to go find
saveRDS(dat_diet_EM, "builds/dat_diet_EM.cluster.rds")
saveRDS(dat_social_CT, "builds/dat_social_CT.cluster.rds")
saveRDS(dat_social_EM, "builds/dat_social_EM.cluster.rds")
saveRDS(dat_urb_CT, "builds/dat_urb_CT.cluster.rds")
saveRDS(dat_urb_EM, "builds/dat_urb_EM.cluster.rds")


#
diet_EM_results #good models = 1, 2, 3, 4
model_indices <- c(1, 2, 3, 4) 

for (i in seq_along(model_indices)) {
  idx <- model_indices[i]
  
  # Extract the model
  model <- get.models(diet_EM_results, subset = idx)[[1]]
  
  # Save it to .rds file
  saveRDS(model, file = file.path("C:/Users/sager/OneDrive/Desktop/school/MSc/SurveyInfo/Recon_Maps/Scat2/2025ScatProject/builds/finalmods/cluster/diet/EM", paste0(i, ".cluster.rds")))
  
  cat("Saved model", idx, "as", i, ".rds\n")
}



#
social_EM_results #good models = 9,10
model_indices <- c(9,10) 

for (i in seq_along(model_indices)) {
  idx <- model_indices[i]
  
  # Extract the model
  model <- get.models(social_EM_results, subset = idx)[[1]]
  
  # Save it to .rds file
  saveRDS(model, file = file.path("C:/Users/sager/OneDrive/Desktop/school/MSc/SurveyInfo/Recon_Maps/Scat2/2025ScatProject/builds/finalmods/cluster/social/EM", paste0(i, ".cluster.rds")))
  
  cat("Saved model", idx, "as", i, ".rds\n")
}



#
social_CT_results #good models = 1,3
model_indices <- c(1,3) 

for (i in seq_along(model_indices)) {
  idx <- model_indices[i]
  
  # Extract the model
  model <- get.models(social_CT_results, subset = idx)[[1]]
  
  # Save it to .rds file
  saveRDS(model, file = file.path("C:/Users/sager/OneDrive/Desktop/school/MSc/SurveyInfo/Recon_Maps/Scat2/2025ScatProject/builds/finalmods/cluster/social/CT", paste0(i, ".cluster.rds")))
  
  cat("Saved model", idx, "as", i, ".rds\n")
}




#
urb_EM_results #good models = 1

model_indices <- c(1) 

for (i in seq_along(model_indices)) {
  idx <- model_indices[i]
  
  # Extract the model
  model <- get.models(urb_EM_results, subset = idx)[[1]]
  
  # Save it to .rds file
  saveRDS(model, file = file.path("C:/Users/sager/OneDrive/Desktop/school/MSc/SurveyInfo/Recon_Maps/Scat2/2025ScatProject/builds/finalmods/cluster/urb/EM", paste0(i, ".cluster.rds")))
  
  cat("Saved model", idx, "as", i, ".rds\n")
}


#
urb_CT_results #good models = 1,2

model_indices <- c(1,2) 

for (i in seq_along(model_indices)) {
  idx <- model_indices[i]
  
  # Extract the model
  model <- get.models(urb_CT_results, subset = idx)[[1]]
  
  # Save it to .rds file
  saveRDS(model, file = file.path("C:/Users/sager/OneDrive/Desktop/school/MSc/SurveyInfo/Recon_Maps/Scat2/2025ScatProject/builds/finalmods/cluster/urb/CT", paste0(i, ".cluster.rds")))
  
  cat("Saved model", idx, "as", i, ".rds\n")
}



