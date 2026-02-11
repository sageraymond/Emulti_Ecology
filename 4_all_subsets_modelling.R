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

dat <- readRDS("builds/prepared_data1.rds")
dat <- as.data.frame(dat)

#Start by just looking at category
dat_category_EM <- dat %>% dplyr::select(ClusterID, Emulti, Category, year)
dat_category_EM <- na.omit(dat_category_EM)

Category_EM<-glmmTMB(Emulti ~ (1|ClusterID) + Category,
                 data = dat_category_EM, 
                 family = binomial, 
                 na.action = "na.fail")

Year_EM<-glmmTMB(Emulti ~ (1|ClusterID) + year,
                     data = dat_category_EM, 
                     family = binomial, 
                     na.action = "na.fail")


dat_category_CT <- dat %>% dplyr::select(ClusterID, CT, Category, year)
dat_category_CT <- na.omit(dat_category_CT)

Category_CT<-glmmTMB(CT ~ (1|ClusterID) + Category,
                     data = dat_category_CT, 
                     family = gaussian, 
                     na.action = "na.fail")

Year_CT<-glmmTMB(CT ~ (1|ClusterID) + year,
                 data = dat_category_CT, 
                 family = gaussian, 
                 na.action = "na.fail")

summary(Category_EM)
summary(Category_CT)
confint(Category_CT)


summary(Year_EM)
summary(Year_CT)

#Start with DIET + Emulti
#From3 . Infection: berries L, garbage Q, nat prey L, veg Q. 
dat_diet_EM <- dat %>%
  dplyr::select(ClusterID, Emulti, Berries, Nat_Prey) %>%
  dplyr::mutate(Berries_scaled = scale(Berries),
                Nat_Prey_scaled = scale(Nat_Prey))

dat_diet_EM <- na.omit(dat_diet_EM) #n = 490

#Check for corelated dudes
cor(dat_diet_EM[sapply(dat_diet_EM, is.numeric)], method = c("pearson"), use = "pairwise.complete.obs")
#No correlation concerns

# fit model with all parameters
diet_EM<-glmmTMB(Emulti ~ (1|ClusterID) + Berries_scaled  + 
                   Nat_Prey_scaled,
                data = dat_diet_EM, 
                family = binomial, 
                na.action = "na.fail")

#All subsets
diet_EM_results <- dredge(diet_EM, 
                          trace = 2)

diet_EM_results

#From3 . Shedding intensity: apples Q, berries L, birdseed Q, Unk Anthro L:
dat_diet_CT <- dat %>%
  dplyr::select(ClusterID, CT, Berries, Unk_Anthro) %>%
  dplyr::mutate(Berries_scaled = scale(Berries),
                Unk_Anthro_scaled = scale(Unk_Anthro))

dat_diet_CT <- na.omit(dat_diet_CT) #n = 162

#Check for corelated dudes
cor(dat_diet_CT[sapply(dat_diet_CT, is.numeric)], method = c("pearson"), use = "pairwise.complete.obs")
#No correlation concerns

# fit model with all parameters
diet_CT<-glmmTMB(CT ~ (1|ClusterID) + Berries_scaled + 
                    Unk_Anthro_scaled ,
                 data = dat_diet_CT, 
                 family = gaussian(), 
                 na.action = "na.fail")

#All subsets
diet_CT_results <- dredge(diet_CT, 
                          trace = 2) #3 top models


#OK next up is social stuff

#From3 Infection: prop children (l), prop high income (q), prop low income (l),
dat_social_EM <- dat %>%
  dplyr::select(ClusterID, Emulti, Prop.Under15, prophighincome, proplowincome) %>%
  dplyr::mutate(Prop.Under15_scaled = scale(Prop.Under15),
                prophighincome_scaled = scale(prophighincome),
                proplowincome_scaled = scale(proplowincome))

dat_social_EM <- na.omit(dat_social_EM) #n = 480

#Check for corelated dudes
cor(dat_social_EM[sapply(dat_social_EM, is.numeric)], method = c("pearson"), use = "pairwise.complete.obs")
# low and high are correlated, which is unsuprising

# fit model with all parameters
social_EM<-glmmTMB(Emulti ~ (1|ClusterID) + Prop.Under15_scaled + prophighincome_scaled + 
                     proplowincome_scaled +  
                   I(prophighincome_scaled^2),
                 data = dat_social_EM, 
                 family = binomial, 
                 na.action = "na.fail")


social_EM_results <- dredge(
  social_EM,
  trace = 2,
  subset = dc(`I(prophighincome_scaled^2)`, prophighincome_scaled))
social_EM_results

#From3 . Shedding intensity:  seniors (L), high income (q), low income (q)
dat_social_CT <- dat %>%
  dplyr::select(ClusterID, CT, Prop.Over65,
                proplowincome, prophighincome) %>%
  dplyr::mutate(Prop.Over65_scaled = scale(Prop.Over65),
                proplowincome_scaled = scale(proplowincome),
                prophighincome_scaled = scale(prophighincome))

dat_social_CT <- na.omit(dat_social_CT) #n = 159

#Check for corelated dudes
cor(dat_social_CT[sapply(dat_social_CT, is.numeric)], method = c("pearson"), use = "pairwise.complete.obs")
#low and high income

# fit model with all parameters
social_CT<-glmmTMB(CT ~ (1|ClusterID)  + Prop.Over65_scaled + 
                      proplowincome_scaled + prophighincome_scaled +
                   I(proplowincome_scaled^2) + I(prophighincome_scaled^2),
                 data = dat_social_CT, 
                 family = gaussian(), 
                 na.action = "na.fail")

#All subsets
social_CT_results <- dredge(social_CT, 
                          trace = 2)
social_CT_results


#Finally urbanization
#From3 anth cover 250 (q), Nat. Area Dec. (0005) (q), Rd. Dec. (0005) (q), dst ravine (q), rd dens 250 (q), 
dat_urb_EM <- dat %>%
  dplyr::select(ClusterID, Emulti, cover_Anth_250m, decay_natarea_0005, 
                decay_road_0005, dist_to_Ravine_m, road_dens_250mH) %>%
  dplyr::mutate(cover_Anth_250m_scaled = scale(cover_Anth_250m),
                decay_natarea_0005_scaled = scale(decay_natarea_0005),
                decay_road_0005_scaled = scale(decay_road_0005),
                dist_to_Ravine_m_scaled = scale(dist_to_Ravine_m),
                road_dens_250mH_scaled = scale(road_dens_250mH))

dat_urb_EM <- na.omit(dat_urb_EM) #n = 490

#Check for corelated dudes
cor(dat_urb_EM[sapply(dat_urb_EM, is.numeric)], method = c("pearson"), use = "pairwise.complete.obs")
# anth cover and road density
#dec road and road density


# fit model with all parameters
urb_EM<-glmmTMB(Emulti ~ (1|ClusterID) + cover_Anth_250m_scaled + decay_natarea_0005_scaled + 
                  decay_road_0005_scaled +  dist_to_Ravine_m_scaled + road_dens_250mH_scaled +
                     I(decay_natarea_0005_scaled^2) +
                I(cover_Anth_250m_scaled^2) +
                  I(decay_road_0005_scaled^2) +
                  I(dist_to_Ravine_m_scaled^2) +
                  I(road_dens_250mH_scaled^2) , 
                   data = dat_urb_EM, 
                   family = binomial, 
                   na.action = "na.fail")


urb_EM_results <- dredge(urb_EM,
  trace = 2)

urb_EM_results

#From3 . Shedding intensity: 
#Bldg. Dens. (250) Q,  Anth. Cover (50) Q, For. Cover (50) Q, Bldg. Dec. (0005) Q, 
#C.B. Dec. (00015) Q, C.C. Dec. (00015) L, Ravine Dec. (00015) L, River Dec. (00015) Q, 
#Dist. Nat Area L,  Rd. Dens. (50) Q, Urbanization PCA Q
dat_urb_CT <- dat %>%
  dplyr::select(ClusterID, CT, bldg_dens_25mH,
                decay_Ravine_002, dist_to_NatArea_m,
                road_dens_250mH) %>%
  dplyr::mutate(dist_to_NatArea_m_scaled = scale(dist_to_NatArea_m),
                bldg_dens_25mH_scaled = scale(bldg_dens_25mH),
                decay_Ravine_002_scaled = scale(decay_Ravine_002),
                road_dens_250mH_scaled = scale(road_dens_250mH))


dat_urb_CT <- na.omit(dat_urb_CT) #n = 162

#Check for corelated dudes
cor(dat_urb_CT[sapply(dat_urb_CT, is.numeric)], method = c("pearson"), use = "pairwise.complete.obs")

# fit model with all parameters
urb_CT<-glmmTMB(CT ~ (1|ClusterID) + bldg_dens_25mH_scaled  + 
                decay_Ravine_002_scaled  + dist_to_NatArea_m_scaled +
                road_dens_250mH_scaled,
                   data = dat_urb_CT, 
                   family = gaussian(), 
                   na.action = "na.fail")
#All subsets
urb_CT_results <- dredge(urb_CT, 
                            trace = 2,
                         m.lim = c(0,12))

  
#Now I want to save within 2AIC models in Rds manner so they are easy to go find
saveRDS(dat_diet_CT, "builds/dat_diet_CT.rds")
saveRDS(dat_diet_EM, "builds/dat_diet_EM.rds")
saveRDS(dat_social_CT, "builds/dat_social_CT.rds")
saveRDS(dat_social_EM, "builds/dat_social_EM.rds")
saveRDS(dat_urb_CT, "builds/dat_urb_CT.rds")
saveRDS(dat_urb_EM, "builds/dat_urb_EM.rds")


#
diet_EM_results #good models = 1, 2, 3
model_indices <- c(1, 2, 3) 

for (i in seq_along(model_indices)) {
  idx <- model_indices[i]
  
  # Extract the model
  model <- get.models(diet_EM_results, subset = idx)[[1]]
  
  # Save it to .rds file
  saveRDS(model, file = file.path("C:/Users/sager/OneDrive/Desktop/school/MSc/SurveyInfo/Recon_Maps/Scat2/2025ScatProject/builds/finalmods/diet/EM", paste0(i, ".rds")))
  
  cat("Saved model", idx, "as", i, ".rds\n")
}


#And CT
diet_CT_results #good models = 1, 2, 3
model_indices <- c(1, 2, 3) 

for (i in seq_along(model_indices)) {
  idx <- model_indices[i]
  
  # Extract the model
  model <- get.models(diet_CT_results, subset = idx)[[1]]
  
  # Save it to .rds file
  saveRDS(model, file = file.path("C:/Users/sager/OneDrive/Desktop/school/MSc/SurveyInfo/Recon_Maps/Scat2/2025ScatProject/builds/finalmods/diet/CT", paste0(i, ".rds")))
  
  cat("Saved model", idx, "as", i, ".rds\n")
}




#
social_EM_results #good models = 1, 2, 3
model_indices <- c(1, 2, 3) 

for (i in seq_along(model_indices)) {
  idx <- model_indices[i]
  
  # Extract the model
  model <- get.models(social_EM_results, subset = idx)[[1]]
  
  # Save it to .rds file
  saveRDS(model, file = file.path("C:/Users/sager/OneDrive/Desktop/school/MSc/SurveyInfo/Recon_Maps/Scat2/2025ScatProject/builds/finalmods/social/EM", paste0(i, ".rds")))
  
  cat("Saved model", idx, "as", i, ".rds\n")
}



#
social_CT_results #good models = 2, 3, 11, 12, 13, 14
model_indices <- c(2, 3, 11, 12, 13, 14) 

for (i in seq_along(model_indices)) {
  idx <- model_indices[i]
  
  # Extract the model
  model <- get.models(social_CT_results, subset = idx)[[1]]
  
  # Save it to .rds file
  saveRDS(model, file = file.path("C:/Users/sager/OneDrive/Desktop/school/MSc/SurveyInfo/Recon_Maps/Scat2/2025ScatProject/builds/finalmods/social/CT", paste0(i, ".rds")))
  
  cat("Saved model", idx, "as", i, ".rds\n")
}




#
urb_EM_results #good models = 5, 9, 13, 15, 34, 37, 58, 63, 67,80,81,86,88,93,97,100,101
# anth cover and road density
#dec road and road density

model_indices <- c(5, 9, 13, 15, 34, 37, 58, 63, 67,80,81,86,88,93,97,100,101) 

for (i in seq_along(model_indices)) {
  idx <- model_indices[i]
  
  # Extract the model
  model <- get.models(urb_EM_results, subset = idx)[[1]]
  
  # Save it to .rds file
  saveRDS(model, file = file.path("C:/Users/sager/OneDrive/Desktop/school/MSc/SurveyInfo/Recon_Maps/Scat2/2025ScatProject/builds/finalmods/urb/EM", paste0(i, ".rds")))
  
  cat("Saved model", idx, "as", i, ".rds\n")
}


#
urb_CT_results #good models = 1 through 14, minus 2

model_indices <- c(1,3,4,5,6,7,8,9,10,11,12,13,14) 

for (i in seq_along(model_indices)) {
  idx <- model_indices[i]
  
  # Extract the model
  model <- get.models(urb_CT_results, subset = idx)[[1]]
  
  # Save it to .rds file
  saveRDS(model, file = file.path("C:/Users/sager/OneDrive/Desktop/school/MSc/SurveyInfo/Recon_Maps/Scat2/2025ScatProject/builds/finalmods/urb/CT", paste0(i, ".rds")))
  
  cat("Saved model", idx, "as", i, ".rds\n")
}



