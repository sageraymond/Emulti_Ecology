#
# AIM: plot model coefficients (scaled)
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
library(glmmTMB)

dat <- readRDS("builds/prepared_data1.rds")


#1.Build category models------------------------------------------------------
Category_EM<-glmmTMB(Emulti ~ (1|ClusterID) + Category,
                     data = dat, 
                     family = binomial, 
                     na.action = "na.omit")
datem <- model.frame(Category_EM)
Category_EM1<-glmmTMB(Emulti ~ (1|ClusterID) + Category,
                     data = datem, 
                     family = binomial, 
                     na.action = "na.omit")
Category_EM_Null<-glmmTMB(Emulti ~ (1|ClusterID),
                      data = datem, 
                      family = binomial, 
                      na.action = "na.omit")

anova(Category_EM1, Category_EM_Null) # 0.0185

Category_CT<-glmmTMB(CT ~ (1|ClusterID) + Category,
                     data = dat, 
                     family = gaussian, 
                     na.action = "na.omit")
datct <- model.frame(Category_CT)
Category_CT1<-glmmTMB(CT ~ (1|ClusterID) + Category,
                      data = datct, 
                      family = gaussian, 
                      na.action = "na.omit")
Category_CT_Null<-glmmTMB(CT ~ (1|ClusterID),
                          data = datct, 
                          family = gaussian, 
                          na.action = "na.omit")

anova(Category_CT1, Category_CT_Null) #0.8185

summary(Category_EM)
confint(Category_EM)
nobs(Category_EM)
r2(Category_EM)

summary(Category_CT)
confint(Category_CT)
nobs(Category_CT)
r2(Category_CT)


#Calculate effect sizes manually : P
ref_level <- levels(dat$Category)[1]

newdat <- data.frame(
  Category = levels(dat$Category),
  ClusterID = NA   
)

pred <- predict(Category_CT, newdata = newdat, se.fit = TRUE)

pred_orig <- exp(pred$fit) - 1

ref_val <- pred_orig[newdat$Category == ref_level]
percent_change <- 100 * (pred_orig - ref_val) / ref_val

#annnnd Cis
crit <- 1.96
lower_log <- pred$fit - crit * pred$se.fit
upper_log <- pred$fit + crit * pred$se.fit

#back transform the puppies!
lower_orig <- exp(lower_log) - 1
upper_orig <- exp(upper_log) - 1

lower_pc <- 100 * (lower_orig - ref_val) / ref_val
upper_pc <- 100 * (upper_orig - ref_val) / ref_val

# Plonk in table
effect_table <- data.frame(
  Category = levels(dat$Category),
  PercentChange = percent_change,
  CI_lower = lower_pc,
  CI_upper = upper_pc
)

effect_table



#2. Make a plot for category--------------------------------------------------
datcat <- as.data.frame(dat)

datcatEM <- datcat %>%
  dplyr::select(Emulti, Category) %>%
  dplyr::group_by(Category, Emulti) %>%
  dplyr::summarise(count = n(), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = Emulti, values_from = count, values_fill = 0) %>%
  dplyr::rename(
    not_infected = `0`,
    infected = `1`
  ) %>%
  dplyr::mutate(
    total = not_infected + infected,
    prop_infected = infected / total
  )
datcatEM <- na.omit(datcatEM)
  

datcatCT <- datcat %>%
  dplyr::select(CT_old, Category) %>%
  dplyr::group_by(Category) %>%
  dplyr::summarise(meanCT = mean(CT_old, na.rm = TRUE),
                   sdCT = sd(CT_old, na.rm = TRUE))

datcatCT <- na.omit(datcatCT)

CatEMPlot <- ggplot(data = datcatEM,
       mapping = aes(x=Category, y=prop_infected)) +
  geom_bar(stat = "identity", colour = "black", position = "stack", fill= "#533A71") +
  geom_text(aes(label = scales::percent(prop_infected, accuracy = 1)), 
            vjust = -0.5, size = 3.5) +
  theme_classic() +
  theme(axis.text.x = element_text(colour = "black", face = "plain", size = 12),
        axis.text.y = element_text(colour = "black", face = "plain", size = 12),
        axis.title.x = element_text(colour = "black", face = "plain", size = 12),
        axis.title.y = element_text(colour = "black", face = "plain", size = 12),
        legend.text = element_text(colour = "black", face = "plain", size = 12),
        legend.title = element_blank()) +
  labs(x="Category", y="Infection Rate") + 
  ggtitle("A. Infection Rate by Category")


CatCTPlot <- 
  ggplot(data = datcatCT,
                    mapping = aes(x=Category, y=meanCT)) +
  geom_bar(stat = "identity", colour = "black", position = "stack", fill= "#533A71") +
  geom_errorbar(aes(ymin = meanCT - sdCT, ymax = meanCT + sdCT), 
                width = 0.2, colour = "black") +
  geom_text(aes(label = round(meanCT, 1)), vjust = -1.5, hjust = 1.5, size = 3.5) +
theme_classic() +
  theme(axis.text.x = element_text(colour = "black", face = "plain", size = 12),
        axis.text.y = element_text(colour = "black", face = "plain", size = 12),
        axis.title.x = element_text(colour = "black", face = "plain", size = 12),
        axis.title.y = element_text(colour = "black", face = "plain", size = 12),
        legend.text = element_text(colour = "black", face = "plain", size = 12),
        legend.title = element_blank()) +
  labs(x="Category", y="Mean Shedding Intensity") + 
  ggtitle("B. Shedding Intensity by Category")

Fig1 <- ggarrange(CatEMPlot, CatCTPlot, nrow = 2)
#ggsave("figures/Fig2.pdf", Fig1, width = 6, height = 8, dpi = 700,  bg = "white") 

#I will later pull out all the info for these dudes


# 2. Bring in all models and extract coefficients--------------------------------------------
folder_path_diet_EM <- "C:/Users/sager/OneDrive/Desktop/school/MSc/SurveyInfo/Recon_Maps/Scat2/2025ScatProject/builds/finalmods/diet/EM"

# 
model_files_diet <- list.files(path = folder_path_diet_EM, pattern = "\\.rds$", full.names = TRUE)

# Read models and tidy
tidy_models_diet_EM <- lapply(seq_along(model_files_diet), function(i) {
  m <- readRDS(model_files_diet[i])
  t <- tryCatch({
    out <- tidy(m, effects = "fixed", conf.int = TRUE)
    out$model_id <- paste0("Model_", i)
    out
  }, error = function(e) NULL)
  return(t)
})

# Combine into one data.table
all_coefs_diet_EM <- rbindlist(tidy_models_diet_EM, fill = TRUE)
all_coefs_diet_EM
all_coefs_diet_EM$cat <- "Diet"
all_coefs_diet_EM$type <- "Infection"

#and for ct
folder_path_diet_CT <- "C:/Users/sager/OneDrive/Desktop/school/MSc/SurveyInfo/Recon_Maps/Scat2/2025ScatProject/builds/finalmods/diet/CT"

# 
model_files_diet <- list.files(path = folder_path_diet_CT, pattern = "\\.rds$", full.names = TRUE)

# Read models and tidy
tidy_models_diet_CT <- lapply(seq_along(model_files_diet), function(i) {
  m <- readRDS(model_files_diet[i])
  t <- tryCatch({
    out <- tidy(m, effects = "fixed", conf.int = TRUE)
    out$model_id <- paste0("Model_", i)
    out
  }, error = function(e) NULL)
  return(t)
})

# Combine into one data.table
all_coefs_diet_CT <- rbindlist(tidy_models_diet_CT, fill = TRUE)
all_coefs_diet_CT
all_coefs_diet_CT$cat <- "Diet"
all_coefs_diet_CT$type <- "Shedding Intensity"

folder_path_social_EM <- "C:/Users/sager/OneDrive/Desktop/school/MSc/SurveyInfo/Recon_Maps/Scat2/2025ScatProject/builds/finalmods/social/EM"

# 
model_files_social <- list.files(path = folder_path_social_EM, pattern = "\\.rds$", full.names = TRUE)

# Read models and tidy
tidy_models_social_EM <- lapply(seq_along(model_files_social), function(i) {
  m <- readRDS(model_files_social[i])
  t <- tryCatch({
    out <- tidy(m, effects = "fixed", conf.int = TRUE)
    out$model_id <- paste0("Model_", i)
    out
  }, error = function(e) NULL)
  return(t)
})

# Combine into one data.table
all_coefs_social_EM <- rbindlist(tidy_models_social_EM, fill = TRUE)
all_coefs_social_EM
all_coefs_social_EM$cat <- "Social"
all_coefs_social_EM$type <- "Infection"

#and for ct
folder_path_social_CT <- "C:/Users/sager/OneDrive/Desktop/school/MSc/SurveyInfo/Recon_Maps/Scat2/2025ScatProject/builds/finalmods/social/CT"

# 
model_files_social <- list.files(path = folder_path_social_CT, pattern = "\\.rds$", full.names = TRUE)

# Read models and tidy
tidy_models_social_CT <- lapply(seq_along(model_files_social), function(i) {
  m <- readRDS(model_files_social[i])
  t <- tryCatch({
    out <- tidy(m, effects = "fixed", conf.int = TRUE)
    out$model_id <- paste0("Model_", i)
    out
  }, error = function(e) NULL)
  return(t)
})

# Combine into one data.table
all_coefs_social_CT <- rbindlist(tidy_models_social_CT, fill = TRUE)
all_coefs_social_CT
all_coefs_social_CT$cat <- "Social"
all_coefs_social_CT$type <- "Shedding Intensity"


folder_path_urb_EM <- "C:/Users/sager/OneDrive/Desktop/school/MSc/SurveyInfo/Recon_Maps/Scat2/2025ScatProject/builds/finalmods/urb/EM"

# 
model_files_urb <- list.files(path = folder_path_urb_EM, pattern = "\\.rds$", full.names = TRUE)

# Read models and tidy
tidy_models_urb_EM <- lapply(seq_along(model_files_urb), function(i) {
  m <- readRDS(model_files_urb[i])
  t <- tryCatch({
    out <- tidy(m, effects = "fixed", conf.int = TRUE)
    out$model_id <- paste0("Model_", i)
    out
  }, error = function(e) NULL)
  return(t)
})

# Combine into one data.table
all_coefs_urb_EM <- rbindlist(tidy_models_urb_EM, fill = TRUE)
all_coefs_urb_EM
all_coefs_urb_EM$cat <- "Urb."
all_coefs_urb_EM$type <- "Infection"

#and for ct
folder_path_urb_CT <- "C:/Users/sager/OneDrive/Desktop/school/MSc/SurveyInfo/Recon_Maps/Scat2/2025ScatProject/builds/finalmods/urb/CT"

# 
model_files_urb <- list.files(path = folder_path_urb_CT, pattern = "\\.rds$", full.names = TRUE)

# Read models and tidy
tidy_models_urb_CT <- lapply(seq_along(model_files_urb), function(i) {
  m <- readRDS(model_files_urb[i])
  t <- tryCatch({
    out <- tidy(m, effects = "fixed", conf.int = TRUE)
    out$model_id <- paste0("Model_", i)
    out
  }, error = function(e) NULL)
  return(t)
})

# Combine into one data.table
all_coefs_urb_CT <- rbindlist(tidy_models_urb_CT, fill = TRUE)
all_coefs_urb_CT
all_coefs_urb_CT$cat <- "Urb."
all_coefs_urb_CT$type <- "Shedding Intensity"


#Bind em all
modeldat <- rbind(all_coefs_diet_EM, all_coefs_diet_CT,
                  all_coefs_social_EM, all_coefs_social_CT,
                  all_coefs_urb_EM, all_coefs_urb_CT)


modeldat <- modeldat[term != "(Intercept)"]
unique(modeldat$term)

modeldat[, term := fcase(
  term == "Berries_scaled", "Berries (%)",
  term == "Garbage_scaled", "Garbage (%)",
  term == "I(Garbage_scaled^2)", "Garbage (%; squared)",
  term == "Vegetation_scaled", "Vegetation (%)",
  term == "Nat_Prey_scaled", "Nat. Prey (%)",
  term == "I(Vegetation_scaled^2)", "Vegetation (%; squared)",
  term == "Birdseed_scaled", "Birdseed (%)",
  term == "I(Birdseed_scaled^2)", "Birdseed (%; squared)",
  term == "Unk_Anthro_scaled", "Unk. Anthro. (%)",
  term == "Apples_scaled", "Apples (%)",
  term == "prophighincome_scaled", "Prop. High Income",
  term == "I(prophighincome_scaled^2)", "Prop. High Income (squared)",
  term == "Prop.Under15_scaled", "Prop. Children",
  term == "cover_Anth_250m_scaled", "Anthro. Cover (%; 250-m)",
  term == "decay_road_0005_scaled", "Decay. Dist to Road",
  term == "I(decay_road_0005_scaled^2)", "Decay. Dist to Road (squared)",
  term == "I(cover_Anth_250m_scaled^2)", "Anthro. Cover (%; 250-m; squared)",
  term == "road_dens_250mH_scaled", "Road Density (250-m)",
  term == "decay_natarea_0005_scaled", "Dec. Dist to Nat. Area",
  term == "I(decay_natarea_0005_scaled^2)", "Dec. Dist to Nat. Area (squared)",
  term == "decay_Ravine_002_scaled", "Dec. Dist. to Ravine",
  term == "dist_to_Ravine_m_scaled", "Dist. to Ravine (m)",
  term == "I(dist_to_Ravine_m_scaled^2)", "Dist to Ravine (m; squared)",
  term == "bldg_dens_25mH_scaled", "Building Density (25-m)",
  term == "bldg_dens_250mH_scaled", "Building Density (250-m)",
  term == "I(bldg_dens_250mH_scaled^2)", "Building Density (250-m; squared)",
  term == "cover_Anth_50m_scaled", "Anthro. Cover (%; 50-m)",
  term == "I(cover_Anth_50m_scaled^2)", "Anthro. Cover (%; 50-m; squared)",
  term == "decay_RivEdge_00015_scaled", "Dec. Dist. to River",
  term == "dist_to_NatArea_m_scaled", "Dist. to Nat. Area (m)",
  term == "I(decay_RivEdge_00015_scaled^2)", "Dec. Dist to River (squared)",
  term == "decay_citybound_00015_scaled", "Dec. Dist to City Boundary",
  term == "I(proplowincome_scaled^2)", "Prop. Low Income (squared)",
  term == "Prop.Over65_scaled", "Prop. Seniors",
  term == "proplowincome_scaled", "Prop. Low Income",
  default = term
)]

modeldat$Plot <- paste(modeldat$cat, modeldat$type) 


#3. Plot these dudes---------------------------------------------------------
# Start with e multi + diet
DietEM <- 
  ggplot(data = modeldat[Plot == "Diet Infection"],
                  aes(x = term, y = estimate, group = model_id)) +
  geom_point(position = position_dodge(width = 0.5), size = 5, 
             color = "#50C5B7",
             shape = 21,
             fill = "#50C5B7") +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high),
                width = 0.2, position = position_dodge(width = 0.5), size = 1,
                color = "#50C5B7") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_classic() +
  labs(title = "C. Infection + Diet",
       y = "Coefficient Estimate", x = "Predictor") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, colour = "black", face = "plain", size = 12),
        axis.text.y = element_text(colour = "black", face = "plain", size = 12),
        axis.title.y = element_text(colour = "black", face = "plain", size = 12),
        legend.text = element_text(colour = "black", face = "plain", size = 12),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.position = "none",
        strip.text = element_text(size = 12, face = "plain"))

DietCT <- 
  ggplot(data = modeldat[Plot == "Diet Shedding Intensity"],
         aes(x = term, y = estimate, group = model_id)) +
  geom_point(position = position_dodge(width = 0.5), size = 5, 
             color = "#50C5B7",
             shape = 25,
             fill = "#50C5B7") +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high),
                width = 0.2, position = position_dodge(width = 0.5), size = 1,
                color = "#50C5B7") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_classic() +
  labs(title = "D. Shedding Intensity + Diet",
       y = "Coefficient Estimate", x = "Predictor") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, colour = "black", face = "plain", size = 12),
        axis.text.y = element_text(colour = "black", face = "plain", size = 12),
        axis.title.y = element_text(colour = "black", face = "plain", size = 12),
        legend.text = element_text(colour = "black", face = "plain", size = 12),
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none",
        strip.text = element_text(size = 12, face = "plain"))

# urbanization
UrbEM <- 
  ggplot(data = modeldat[Plot == "Urb. Infection"],
         aes(x = term, y = estimate, group = model_id)) +
  geom_point(position = position_dodge(width = 0.5), size = 5, 
             color = "#6184d8",
             shape = 21,
             fill = "#6184d8") +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high),
                width = 0.2, position = position_dodge(width = 0.5), size = 1,
                color = "#6184d8") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_classic() +
  labs(title = "A. Infection + Urbanization",
       y = "Coefficient Estimate", x = "Predictor") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, colour = "black", face = "plain", size = 12),
        axis.text.y = element_text(colour = "black", face = "plain", size = 12),
        axis.title.y = element_text(colour = "black", face = "plain", size = 12),
        legend.text = element_text(colour = "black", face = "plain", size = 12),
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none",
        strip.text = element_text(size = 12, face = "plain"))

UrbCT <- 
  ggplot(data = modeldat[Plot == "Urb. Shedding Intensity"],
         aes(x = term, y = estimate, group = model_id)) +
  geom_point(position = position_dodge(width = 0.5), size = 5, 
             color = "#6184d8",
             shape = 25,
             fill = "#6184d8") +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high),
                width = 0.2, position = position_dodge(width = 0.5), size = 1,
                color = "#6184d8") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_classic() +
  labs(title = "B. Shedding Intensity + Urbanization",
       y = "Coefficient Estimate", x = "Predictor") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, colour = "black", face = "plain", size = 12),
        axis.text.y = element_text(colour = "black", face = "plain", size = 12),
        axis.title.y = element_text(colour = "black", face = "plain", size = 12),
        legend.text = element_text(colour = "black", face = "plain", size = 12),
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none",
        strip.text = element_text(size = 12, face = "plain"))

# SOcial
SocialEM <- 
  ggplot(data = modeldat[Plot == "Social Infection"],
         aes(x = term, y = estimate, group = model_id)) +
  geom_point(position = position_dodge(width = 0.5), size = 5, 
             color = "#9CEC5B",
             shape = 21,
             fill = "#9CEC5B") +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high),
                width = 0.2, position = position_dodge(width = 0.5), size = 1,
                color = "#9CEC5B") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_classic() +
  labs(title = "E. Infection + Demography",
       y = "Coefficient Estimate", x = "Predictor") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, colour = "black", face = "plain", size = 12),
        axis.text.y = element_text(colour = "black", face = "plain", size = 12),
        axis.title.y = element_text(colour = "black", face = "plain", size = 12),
        legend.text = element_text(colour = "black", face = "plain", size = 12),
        legend.title = element_blank(),
        legend.position = "none",
        strip.text = element_text(size = 12, face = "plain"))

SocialCT <- 
  ggplot(data = modeldat[Plot == "Social Shedding Intensity"],
         aes(x = term, y = estimate, group = model_id)) +
  geom_point(position = position_dodge(width = 0.5), size = 5, 
             color = "#9CEC5B",
             shape = 25,
             fill = "#9CEC5B") +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high),
                width = 0.2, position = position_dodge(width = 0.5), size = 1,
                color = "#9CEC5B") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_classic() +
  labs(title = "F. Shedding Intensity + Demography",
       y = "Coefficient Estimate", x = "Predictor") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, colour = "black", face = "plain", size = 12),
        axis.text.y = element_text(colour = "black", face = "plain", size = 12),
        axis.title.y = element_text(colour = "black", face = "plain", size = 12),
        legend.text = element_text(colour = "black", face = "plain", size = 12),
        legend.title = element_blank(),
        legend.position = "none",
        strip.text = element_text(size = 12, face = "plain"))



AB <- ggarrange(UrbEM, UrbCT, nrow = 1)
CD <- ggarrange(DietEM, DietCT, nrow = 1)
EF <- ggarrange(SocialEM, SocialCT, nrow = 1)
FigS1 <- annotate_figure(
  ggarrange(AB, CD, EF, nrow = 3, heights = c(5, 4, 4)),
  left = text_grob("                                    ", rot = 110)  
)
#ggsave("figures/FigS1_Nov20.pdf", FigS1, width = 8, height = 12, dpi = 700,  bg = "white") 
