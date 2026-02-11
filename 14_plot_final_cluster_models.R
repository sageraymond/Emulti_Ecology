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

# 2. Bring in all models and extract coefficients--------------------------------------------
folder_path_diet_EM <- "C:/Users/sager/OneDrive/Desktop/school/MSc/SurveyInfo/Recon_Maps/Scat2/2025ScatProject/builds/finalmods/cluster/diet/EM"

# 
model_files_diet <- list.files(path = folder_path_diet_EM, pattern = "\\.cluster.rds$", full.names = TRUE)

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



folder_path_social_EM <- "C:/Users/sager/OneDrive/Desktop/school/MSc/SurveyInfo/Recon_Maps/Scat2/2025ScatProject/builds/finalmods/cluster/social/EM"

# 
model_files_social <- list.files(path = folder_path_social_EM, pattern = "\\.cluster.rds$", full.names = TRUE)

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
folder_path_social_CT <- "C:/Users/sager/OneDrive/Desktop/school/MSc/SurveyInfo/Recon_Maps/Scat2/2025ScatProject/builds/finalmods/cluster/social/CT"

# 
model_files_social <- list.files(path = folder_path_social_CT, pattern = "\\.cluster.rds$", full.names = TRUE)

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


folder_path_urb_EM <- "C:/Users/sager/OneDrive/Desktop/school/MSc/SurveyInfo/Recon_Maps/Scat2/2025ScatProject/builds/finalmods/cluster/urb/EM"

# 
model_files_urb <- list.files(path = folder_path_urb_EM, pattern = "\\.cluster.rds$", full.names = TRUE)

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
folder_path_urb_CT <- "C:/Users/sager/OneDrive/Desktop/school/MSc/SurveyInfo/Recon_Maps/Scat2/2025ScatProject/builds/finalmods/cluster/urb/CT"

# 
model_files_urb <- list.files(path = folder_path_urb_CT, pattern = "\\.cluster.rds$", full.names = TRUE)

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
modeldat <- rbind(all_coefs_diet_EM,
                  all_coefs_social_EM, all_coefs_social_CT,
                  all_coefs_urb_EM, all_coefs_urb_CT)


modeldat <- modeldat[term != "(Intercept)"]
unique(modeldat$term)

modeldat[, term := fcase(
  term == "compost_scaled", "Compost scats",
  term == "propAnthro_scaled", "Unk. anthro (prop.)",
  term == "propApple_scaled", "Apples (prop.)",
  term == "propSeed_scaled", "Birdseed (prop.)",
  term == "propGarbage_scaled", "Garbage (prop.)",
  term == "lowincome_scaled", "Prop. low income",
  term == "Over65_scaled", "Prop. seniors",
  term == "I(Over65_scaled^2)", "Prop. seniors (quad.)",
  term == "NoCertDiplomadegree_scaled", "Prop. no post secondary",
  term == "highincome_scaled", "Prop. high income",
  term == "rd_density_scaled", "Road density (m/ Ha)",
  term == "dist_to_RivEdge_m_scaled", "Dist. to river (m)",
  term == "decay_citybound_00015_scaled", "Decay. dist. to city boundary",
  term == "decay_road_0005_scaled", "Decay. dist to road",
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
  labs(title = "C. Infection +\nDiet",
       y = "Coefficient Estimate", x = "Predictor") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, colour = "black", face = "plain", size = 12),
        axis.text.y = element_text(colour = "black", face = "plain", size = 12),
        axis.title.y = element_text(colour = "black", face = "plain", size = 12),
        legend.text = element_text(colour = "black", face = "plain", size = 12),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
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
  labs(title = "A. Infection +\nUrbanization",
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
  labs(title = "B. Shedding Intensity +\nUrbanization",
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
  labs(title = "D. Infection +\nDemography",
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
  labs(title = "E. Shedding Intensity +\nDemography",
       y = "Coefficient Estimate", x = "Predictor") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, colour = "black", face = "plain", size = 12),
        axis.text.y = element_text(colour = "black", face = "plain", size = 12),
        axis.title.y = element_text(colour = "black", face = "plain", size = 12),
        legend.text = element_text(colour = "black", face = "plain", size = 12),
        legend.title = element_blank(),
        legend.position = "none",
        strip.text = element_text(size = 12, face = "plain"))

blank <- ggplot() +
  theme_void()

ABC <- ggarrange(UrbEM, UrbCT, DietEM, nrow = 1)
EFblank <- ggarrange(SocialEM, SocialCT, blank, nrow = 1)
FigS2 <- annotate_figure(
  ggarrange(ABC, EFblank, nrow = 2, heights = c(4, 4)),
  left = text_grob("                                    ", rot = 110)  
)
#ggsave("figures/FigS2_Nov26.pdf", FigS2, width = 10, height = 8, dpi = 700,  bg = "white") 
