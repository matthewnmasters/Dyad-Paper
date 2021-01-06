#Combining only the info we need for the draft Kassandra sent

###############################################################################################################################################################
# The collaborators feel it is appropriate to move forward with Stavrova's approach using the robust SE method after comparing it to frailty models           #
# We believe that the PH Assumption is satisfied by looking at two time points of 18 years.                                                                   #
# Decided not to look at discrete time or time-varying approaches                                                                                             #
###############################################################################################################################################################

library(survival) #for survival analysis
library(tidyverse) #for possible data manipulation
library(lme4) #linear mixed effects
library(lmerTest) #adds p values to lme4, and is simply a "wrapper"
library(psych) #for Kappa
library(gmodels) #for some crosstabs

setwd("S:/USER/MMASTERS")
finalData<- readRDS(file = "finaldata.Rds") 

finalDataPH<- finalData %>%
  mutate(DB82 = as.factor(DB82)) %>% #Diabetes at baseline: 0,1
  mutate(SEX = factor(SEX, levels=c(1,0))) %>% #Flip it so it models being male (0) vs female (1)
  mutate(actIsoColl = ifelse(BSINDEX_3COMP %in% c(2,3), 1, 0)) %>% #collapse actor isolation
  mutate(spouseIsoColl = ifelse(SpousalIsoScore %in% c(2,3), 1, 0)) %>% #collapse partner isolation
  mutate(finalGEO = as.character(urbanRural)) %>% #imputing missing values where possible and create a mismatch level
  mutate(finalGEO = ifelse(urbanRural=="MISSING", as.character(spouseUR), finalGEO)) %>%
  mutate(finalGEO = ifelse(urbanRural=="RURAL", ifelse(spouseUR=="URBAN", "MISMATCH", finalGEO),finalGEO)) %>%
  mutate(finalGEO = ifelse(urbanRural=="URBAN", ifelse(spouseUR=="RURAL", "MISMATCH", finalGEO),finalGEO)) %>%
  mutate(finalGEO = factor(finalGEO, levels = c("URBAN", "RURAL", "MISMATCH", "MISSING"))) %>%
  mutate(ISOGROUP = factor(ISOGROUP, levels = c( "NOISO", "MAN", "FEM", "BOTH"))) %>% #create dyad type variables for potential analysis stratified by sex
  mutate(NewOCC = ifelse(OCCNOW ==0, ifelse(OCCLIFE==7, 3, OCCNOW), OCCNOW)) %>%
  mutate(NewOCC = factor(NewOCC)) %>%
  mutate(SERegion = ifelse(STATE %in% c("TX", "LA", "OK", "AR", "MS", "TN", "KY", "AL", "GA", "FL", "SC", "NC", "VA", "WV", "MD", "DE", "PR", "DC", "VI"),1,0)) %>%
  mutate(WRegion = ifelse(STATE %in% c("CA", "OR", "WA", "ID", "NV", "AK", "HI", "GU", "UT", "AZ", "NM", "CO", "WY", "MT" ),1,0)) %>%
  mutate(MWRegion = ifelse(STATE %in% c("ND", "SD", "NE", "KS", "MN", "IA", "MO", "WI", "IL", "MI", "IN", "OH" ),1,0)) %>%
  mutate(NERegion = ifelse(STATE %in% c("ME", "NH", "MA", "VT", "RI", "CT", "NY", "NJ", "PA" ),1,0)) %>%
  mutate(CensusRegion = ifelse(SERegion ==1, "SouthEast", ifelse(WRegion ==1, "West", ifelse(MWRegion==1, "MidWest", ifelse(NERegion ==1, "NorthEast", "MISSING"))))) %>%
  mutate(CensusRegion = factor(CensusRegion, levels= c("SouthEast", "NorthEast", "MidWest", "West", "MISSING")))


#Cohen's Kappa
men <- finalDataPH %>% filter(SEX==1) #arbitrary; could have chosen women
tabledIso <- table(men$actIsoColl, men$spouseIsoColl)
cohen.kappa(tabledIso)

#Getting an ICC with the ICC package and also with lmer running a null model#
ICCbareF(ID13,BSINDEX_3COMP, data=finalDataPH) #ICC 0.4253568
iccmodel<- lmer(BSINDEX_3COMP ~ 1 + (1|ID13), data=finalDataPH, REML=FALSE)
summary(iccmodel)
(0.2670)/(0.2670+0.3599) #=0.4259, virtually identical to the ICC() method.

# Both of these methods indicate that we should in fact be looking at this as a multilevel model, since we have research showing isolation and mortality are related
#and here we see they are highly correlated.

##############################################################################################################################
# Table 1. descriptive characteristics                                                                                      #
##############################################################################################################################
menOnly<- finalDataPH %>% filter(SEX==0)
womenOnly<- finalDataPH %>% filter(SEX==1)

menOnlyNonIso<- finalDataPH %>% filter(SEX==0, actIsoColl==0)
menOnlyIso<- finalDataPH %>% filter(SEX==0, actIsoColl==1)
womenOnlyNonIso<- finalDataPH %>% filter(SEX==1, actIsoColl==0)
womenOnlyIso<- finalDataPH %>% filter(SEX==1, actIsoColl==1)

#age
t.test(menOnlyNonIso$AGE_INT, menOnlyIso$AGE_INT)
sd(menOnlyNonIso$AGE_INT)
sd(menOnlyIso$AGE_INT)

t.test(womenOnlyNonIso$AGE_INT, womenOnlyIso$AGE_INT)
sd(womenOnlyNonIso$AGE_INT)
sd(womenOnlyIso$AGE_INT)

#race

CrossTable(menOnly$RaceNew, menOnly$actIsoColl, prop.r = FALSE, prop.chisq = FALSE, prop.t=FALSE, chisq = TRUE)
CrossTable(womenOnly$RaceNew, womenOnly$actIsoColl, prop.r = FALSE, prop.chisq = FALSE, prop.t=FALSE, chisq = TRUE)

#education (2 or higher )
CrossTable(menOnly$EDUNEW, menOnly$actIsoColl, prop.r = FALSE, prop.chisq = FALSE, prop.t=FALSE, chisq = TRUE)
CrossTable(womenOnly$EDUNEW, womenOnly$actIsoColl, prop.r = FALSE, prop.chisq = FALSE, prop.t=FALSE, chisq = TRUE)

#smoking
CrossTable(menOnly$SMKGRP, menOnly$actIsoColl, prop.r = FALSE, prop.chisq = FALSE, prop.t=FALSE, chisq = TRUE)
CrossTable(womenOnly$SMKGRP, womenOnly$actIsoColl, prop.r = FALSE, prop.chisq = FALSE, prop.t=FALSE, chisq = TRUE)

#diabetes
CrossTable(menOnly$DB82, menOnly$actIsoColl, prop.r = FALSE, prop.chisq = FALSE, prop.t=FALSE, chisq = TRUE)
CrossTable(womenOnly$DB82, womenOnly$actIsoColl, prop.r = FALSE, prop.chisq = FALSE, prop.t=FALSE, chisq = TRUE)

#BMI
CrossTable(menOnly$BMI, menOnly$actIsoColl, prop.r = FALSE, prop.chisq = FALSE, prop.t=FALSE, chisq = TRUE)
CrossTable(womenOnly$BMI, womenOnly$actIsoColl, prop.r = FALSE, prop.chisq = FALSE, prop.t=FALSE, chisq = TRUE)

#interracial couple status
CrossTable(menOnly$interracial, menOnly$actIsoColl, prop.r = FALSE, prop.chisq = FALSE, prop.t=FALSE, chisq = TRUE)
CrossTable(womenOnly$interracial, womenOnly$actIsoColl, prop.r = FALSE, prop.chisq = FALSE, prop.t=FALSE, chisq = TRUE)

#Residency
CrossTable(menOnly$finalGEO, menOnly$actIsoColl, prop.r = FALSE, prop.chisq = FALSE, prop.t=FALSE, chisq = TRUE)
CrossTable(womenOnly$finalGEO, womenOnly$actIsoColl, prop.r = FALSE, prop.chisq = FALSE, prop.t=FALSE, chisq = TRUE)

##############################################################################################################################
# Table 2. isolation distribution                                                                                            #
##############################################################################################################################

CrossTable(finalDataPH$SEX, finalDataPH$BSINDEX_3COMP, chisq = TRUE)
CrossTable(finalDataPH$SEX, finalDataPH$actIsoColl, chisq = TRUE)

###########################################################################################################################################################
# Split the data into the two 18-year time periods. Can round up very small follow up time variables to help with the timefix = FALSE problem, if needed  #
###########################################################################################################################################################

SplitAC <- survSplit(Surv(FAIL_AGELIM, ALLCAUSE)~ actIsoColl + spouseIsoColl + ISOGROUP + SEX + AGE_INT + NewOCC + SMKGRP + DB82 + EDUNEW + BMI + newDyadID + RaceNew + finalGEO + interracial + OCCNOW + OCCLIFE + CensusRegion + ID + SPOUSEID + ISOGROUP, data=finalDataPH, cut=c(18), episode="half")
firstHalfAC <- SplitAC %>% filter(half==1) %>% mutate(FUPTIME = FAIL_AGELIM-tstart)
secondHalfAC <- SplitAC %>% filter(half==2) %>% mutate(FUPTIME = FAIL_AGELIM-tstart)

#Now for Cancer

SplitCAN <- survSplit(Surv(FAIL_AGELIM, ALLCAN)~ actIsoColl + spouseIsoColl + ISOGROUP + SEX + AGE_INT + NewOCC + SMKGRP + DB82 + EDUNEW + BMI + newDyadID + RaceNew + finalGEO + interracial + OCCNOW +OCCLIFE + CensusRegion + ID + SPOUSEID, data=finalDataPH, cut=c(18), episode="half")
firstHalfCAN <- SplitCAN %>% filter(half==1) %>% mutate(FUPTIME = FAIL_AGELIM-tstart)
secondHalfCAN <- SplitCAN %>% filter(half==2) %>% mutate(FUPTIME = FAIL_AGELIM-tstart)


#And CVD

SplitCVD <- survSplit(Surv(FAIL_AGELIM, ALLVAS)~ actIsoColl + spouseIsoColl + ISOGROUP + SEX + AGE_INT + NewOCC + SMKGRP + DB82 + EDUNEW + BMI + newDyadID + RaceNew + finalGEO + interracial + OCCNOW +OCCLIFE  + CensusRegion + ID + SPOUSEID, data=finalDataPH, cut=c(18), episode="half")
firstHalfCVD <- SplitCVD %>% filter(half==1) %>% mutate(FUPTIME = FAIL_AGELIM-tstart)
secondHalfCVD <- SplitCVD %>% filter(half==2) %>% mutate(FUPTIME = FAIL_AGELIM-tstart)

##############################################################################################################################
# Table 3. sandwich models                                                                                                   #
##############################################################################################################################

Model1AC <- coxph(Surv(FUPTIME, ALLCAUSE) ~ spouseIsoColl + actIsoColl + SEX + strata(AGE_INT)  +SMKGRP + DB82 + EDUNEW + BMI + RaceNew + cluster(newDyadID), data=firstHalfAC)
summary(Model1AC)

Model1CAN <- coxph(Surv(FUPTIME, ALLCAN) ~  spouseIsoColl  + actIsoColl + SEX + strata(AGE_INT) + SMKGRP + DB82 + EDUNEW + BMI + RaceNew + cluster(newDyadID), data=firstHalfCAN)
summary(Model1CAN)

Model1CVD <- coxph(Surv(FUPTIME, ALLVAS) ~  spouseIsoColl + actIsoColl + SEX + strata(AGE_INT) + SMKGRP + DB82 + EDUNEW + BMI + RaceNew + cluster(newDyadID), data=firstHalfCVD)
summary(Model1CVD)

#Now second half of follow-up. 

Model2AC <- coxph(Surv(FUPTIME, ALLCAUSE) ~  spouseIsoColl  + actIsoColl + SEX + strata(AGE_INT) + SMKGRP + DB82 + EDUNEW + BMI + RaceNew + cluster(newDyadID), data=secondHalfAC)
summary(Model2AC)


Model2CAN <- coxph(Surv(FUPTIME, ALLCAN) ~ spouseIsoColl  + actIsoColl + SEX + strata(AGE_INT) + SMKGRP + DB82 + EDUNEW + BMI + RaceNew + cluster(newDyadID), data=secondHalfCAN)
summary(Model2CAN)


Model2CVD <- coxph(Surv(FUPTIME, ALLVAS) ~ spouseIsoColl + actIsoColl + SEX + strata(AGE_INT) + SMKGRP + DB82 + EDUNEW + BMI + RaceNew + cluster(newDyadID), data=secondHalfCVD)
summary(Model2CVD)


##############################################################################################################################
# Table 4. stratified isolation group results                                                                                #
##############################################################################################################################

firstmenAC <- firstHalfAC %>% filter(SEX==0)
firstmenCAN <- firstHalfCAN %>% filter(SEX==0)
firstmenCVD <- firstHalfCVD %>% filter(SEX==0)
secondmenAC <- secondHalfAC %>% filter(SEX==0)
secondmenCAN <- secondHalfCAN %>% filter(SEX==0)
secondmenCVD <- secondHalfCVD %>% filter(SEX==0)
firstwomenAC <- firstHalfAC %>% filter(SEX==1)
firstwomenCAN <- firstHalfCAN %>% filter(SEX==1)
firstwomenCVD <- firstHalfCVD %>% filter(SEX==1)
secondwomenAC <- secondHalfAC %>% filter(SEX==1)
secondwomenCAN <- secondHalfCAN %>% filter(SEX==1)
secondwomenCVD <- secondHalfCVD %>% filter(SEX==1)

#men

Model1AC <- coxph(Surv(FUPTIME, ALLCAUSE) ~ ISOGROUP + SEX + strata(AGE_INT)  +SMKGRP + DB82 + EDUNEW + BMI + RaceNew, data=firstmenAC)
summary(Model1AC)

Model1CAN <- coxph(Surv(FUPTIME, ALLCAN) ~  ISOGROUP + SEX + strata(AGE_INT) + SMKGRP + DB82 + EDUNEW + BMI + RaceNew, data=firstmenCAN)
summary(Model1CAN)

Model1CVD <- coxph(Surv(FUPTIME, ALLVAS) ~  ISOGROUP + SEX + strata(AGE_INT) + SMKGRP + DB82 + EDUNEW + BMI + RaceNew, data=firstmenCVD)
summary(Model1CVD)

#Now second half of follow-up. 

Model2AC <- coxph(Surv(FUPTIME, ALLCAUSE) ~  ISOGROUP + SEX + strata(AGE_INT) + SMKGRP + DB82 + EDUNEW + BMI + RaceNew, data=secondmenAC)
summary(Model2AC)

Model2CAN <- coxph(Surv(FUPTIME, ALLCAN) ~ ISOGROUP + SEX + strata(AGE_INT) + SMKGRP + DB82 + EDUNEW + BMI + RaceNew, data=secondmenCAN)
summary(Model2CAN)

Model2CVD <- coxph(Surv(FUPTIME, ALLVAS) ~ ISOGROUP + SEX + strata(AGE_INT) + SMKGRP + DB82 + EDUNEW + BMI + RaceNew, data=secondmenCVD)
summary(Model2CVD)

#women

Model1AC <- coxph(Surv(FUPTIME, ALLCAUSE) ~ ISOGROUP + SEX + strata(AGE_INT)  +SMKGRP + DB82 + EDUNEW + BMI + RaceNew, data=firstwomenAC)
summary(Model1AC)

Model1CAN <- coxph(Surv(FUPTIME, ALLCAN) ~  ISOGROUP + SEX + strata(AGE_INT) + SMKGRP + DB82 + EDUNEW + BMI + RaceNew, data=firstwomenCAN)
summary(Model1CAN)

Model1CVD <- coxph(Surv(FUPTIME, ALLVAS) ~  ISOGROUP + SEX + strata(AGE_INT) + SMKGRP + DB82 + EDUNEW + BMI + RaceNew, data=firstwomenCVD)
summary(Model1CVD)

#Now second half of follow-up. 

Model2AC <- coxph(Surv(FUPTIME, ALLCAUSE) ~  ISOGROUP + SEX + strata(AGE_INT) + SMKGRP + DB82 + EDUNEW + BMI + RaceNew, data=secondwomenAC)
summary(Model2AC)

Model2CAN <- coxph(Surv(FUPTIME, ALLCAN) ~ ISOGROUP + SEX + strata(AGE_INT) + SMKGRP + DB82 + EDUNEW + BMI + RaceNew, data=secondwomenCAN)
summary(Model2CAN)

Model2CVD <- coxph(Surv(FUPTIME, ALLVAS) ~ ISOGROUP + SEX + strata(AGE_INT) + SMKGRP + DB82 + EDUNEW + BMI + RaceNew, data=secondwomenCVD)
summary(Model2CVD)