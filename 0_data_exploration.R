
# Prepare workspace ---------------------------------------
#

rm(list = ls())
gc()

#Load libs
library(dplyr)
library(glmmTMB)
library(DHARMa)
library(performance)
library(tidyr)
library(tibble)

# Load data ---------------------------------------------------------------
#1. Clean data a bit and save as RDS --------------------------------------
dat <- readRDS("builds/prepared_data1.rds")
count <- dat %>% group_by(ClusterID, Category) %>%
  summarise(n = n())
cat_by_clust <- count %>%
  pivot_wider(
    names_from = Category,  
    values_from = n,        
    values_fill = 0         
  )

cat_by_clust <- as.data.frame(cat_by_clust)
cat_by_clust <- cat_by_clust %>%
  dplyr::select(c(Compost, Garden, ClusterID))

#saveRDS(cat_by_clust, "builds/cat_by_clust.rds")


#dat$CT_old <- dat$CT
#dat$CT <- log(dat$CT_old + 1)

#dat <- read_csv("C:/Users/sager/OneDrive/Desktop/school/MSc/SurveyInfo/Recon_Maps/Scat2/2025ScatProject/data/rawdata/FinalData_Aug16.csv")

#dat <- dat[!is.na(dat$x) & !is.na(dat$y), ]
#dat$pos <- numFactor(dat$x, dat$y)
#dat$group <- factor(1)

#
#str(dat)

#Clean up RDS
#dat$Category <- factor(dat$Category,
 #                      levels = c("Control", "Compost", "Park", 
  #                                "Garden", "Pond", "Ecotone"))

#dat$Emulti <- factor(dat$Emulti, levels = c(0,1))
#dat <- dat %>%
 # dplyr::mutate(CT = 40 - CP)
#dat$CT #looks good

#dat$year <- ifelse(dat$year == "2012", "2015", dat$year)
#dat$year <- factor(dat$year, levels = c(2015, 2017, 2018, 2019, 2020, 2023, 2025))

#dat$CT #looks good

#dat$Apples_pres <- factor(dat$Apples_pres, levels = c(0,1))
#dat$Berries_pres <- factor(dat$Berries_pres, levels = c(0,1))
#dat$Birdseed_pres <- factor(dat$Birdseed_pres, levels = c(0,1))
#dat$Vegetation_pres <- factor(dat$Vegetation_pres, levels = c(0,1))
#dat$Nat_Prey_pres <- factor(dat$Nat_Prey_pres, levels = c(0,1))
#dat$Garbage_pres <- factor(dat$Garbage_pres, levels = c(0,1))
#dat$Unk_Anthro_pres <- factor(dat$Unk_Anthro_pres, levels = c(0,1))

#dat$ClusterID <- as.factor(dat$ClusterID)

#realise i need to correct social vars by pop/ area
#dat <- dat %>%
 #dplyr::mutate(Prop.CompletedDegree = No.CompletedDegree / totalpop,
  #              Prop.NoCertDiplomadegree = No.NoCertDiplomadegree / totalpop,
   #             Prop.CompletedHighschool  =No.CompletedHighschool / totalpop,
    #            Prop.Over65 = No.Over65 / totalpop,
     #           Prop.Under15 = No.Under15 / totalpop) %>%
  #dplyr::select(-c(No.CompletedDegree,
   #                No.NoCertDiplomadegree,
    #               No.CompletedHighschool,
     #              No.Over65,
      #             No.Under15))

#dat$pos <- numFactor(dat$x, dat$y)
#dat$group <- factor(dat$group)      # or just: dat$group <- factor(1)


#saveRDS(dat, "builds/prepared_data1.rds")


#2. Figure out appropriate distribution for CT values---------------------------
#need to figure out how to model CT
n <- glmmTMB(CT ~ Nat_Prey + (1|ClusterID),
             data = dat,
             family = gaussian())
summary(n)
nobs(n) #good; it considered only infected scats

resids <- DHARMa::simulateResiduals(n)
plot(resids)
# it looks pretty good, but it's not totally happy

hist(dat$CT) #it loooks pretty gaussian!


m.disp <- glmmTMB(CT ~ Nat_Prey + (1|ClusterID),
             dispformula = ~ Nat_Prey,
             data = dat,
             family = gaussian())
summary(m.disp)
anova(m, m.disp) #haha, nope!


hist(log(dat$CT +1)) #that looks nicer. 

dat1 <- dat
dat1$CT <- ifelse(dat1$CT == 0, 0.001, dat1$CT)
  
m <- glmmTMB(CT ~ Nat_Prey + (1|ClusterID),
             data = dat1,
             family = lognormal())
summary(m)
nobs(m) #good; it considered only infected scats

resids <- DHARMa::simulateResiduals(m)
plot(resids)
# Oh that's much worse

m.disp <- glmmTMB(CT ~ Nat_Prey + (1|ClusterID),
                  dispformula = ~ Nat_Prey,
                  data = dat1,
                  family = lognormal())
summary(m.disp)
anova(m, m.disp) #here it does improve performance. Hmmm 

#OK. back to gaussian. Check it a bit
m <- glmmTMB(CT ~ Nat_Prey + (1|ClusterID),
             data = dat,
             family = gaussian())
#ZI?
check_overdispersion(m) #no overdispersion...

#Gaussian seems like an option.

#What about log transforming those puppies?
m <- glmmTMB(log(CT + 1) ~ Nat_Prey + (1|ClusterID),
             data = dat,
             family = gaussian())
summary(m)
nobs(m) #good; it considered only infected scats

resids <- DHARMa::simulateResiduals(m)
plot(resids)
# it looks pretty good, but it's not totally happy

AIC(m, n)
anova(m,n)



#What about poisson or NB though? Can't do poisson but didn't consider nbinom
m <- glmmTMB(CT ~ Nat_Prey + (1|ClusterID),
             data = dat1,
             family = nbinom2())
summary(m)
nobs(m) #good; it considered only infected scats

resids <- DHARMa::simulateResiduals(m)
plot(resids)



#3. Extract some summary stuff for first paragraph-------------------------------
dat %>%
  group_by(year) %>% summarise (n = n())

dat %>% 
  group_by(Category) %>% summarise (n = n())

dat %>% 
  group_by(Emulti) %>% summarise (n = n())

min(dat$CT, na.rm = TRUE)
max(dat$CT, na.rm = TRUE)


#does infection fvary among years
chidat <- dat1 %>% group_by(Emulti, year) %>% summarise(n = n())

ctable <- chidat %>%
  pivot_wider(names_from = year, values_from = n, values_fill = 0) %>%
  column_to_rownames("Emulti") %>%
  as.matrix()

chisq_result <- chisq.test(ctable)

chisq_result #X-squared = 4.4347, df = 5, p-value = 0.4887

#does shedding intensity vary among years?
anova_result <- aov(CT ~ year, data = dat1)
summary(anova_result)

dat1 %>%
  group_by(year) %>%
  summarise(mean = mean(CT, na.rm = TRUE))


#nOW consider spatial scat stuff without cluster
mmm <- glmmTMB(Emulti ~ 1 + exp(pos + 0 | group), data=dat, family = binomial)
summary(mmm)
AIC(mmm)

mm <- glmmTMB(Emulti ~ 1, data = dat, family = binomial)
summary(mm)
AIC(mm)

r2(mmm)


#OK, so the spatial piece definitiely improves things, which is unsurprising

re <- ranef(mmm)$cond$group

eee <- glmmTMB(CT ~ 1 + exp(pos + 0 | group), data=dat, family = gaussian, na.action = "na.omit")
summary(eee)
AIC(eee)

#OK, so this suggests that there is no spatial pattern in CT

