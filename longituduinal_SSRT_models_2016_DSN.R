# longitudinal SSRT and pars opercularis morphology
setwd("/Users/laurenbutlercurley/Documents/R_codes/")
library(nlme) #LME
library(mgcv) #GAM
rm(list=ls())
dat <- read.csv("allTP_goodMRI_CANTAB_baseline.csv", header = TRUE)

#singleTP <- dat[!duplicated(dat$SubjID),] 
#write.csv(singleTP, file = "firstTP_goodMRI_CANTAB.csv")

#dat <- read.csv("usercache_PLING_lauren.csv", header = TRUE)

##### create bilateral ROIs #####
dat$MRI_cort_area.ctx.bilateral.parsopercularis <- (dat$MRI_cort_area.ctx.lh.parsopercularis + dat$MRI_cort_area.ctx.rh.parsopercularis)
dat$MRI_cort_thick.ctx.bilateral.parsopercularis <- (dat$MRI_cort_thick.ctx.lh.parsopercularis + dat$MRI_cort_thick.ctx.rh.parsopercularis)/2

##### invert SSRT scores #####
dat$SSRT_inverted <- (dat$BEH_CANTAB_SST_SSRT_last_half * (-1))

###########################
###########################
## REGRESSSION ANALYSES  ##
###########################
###########################

# predicting SSRT from age, age2, gender
m.beh = lme(SSRT_inverted ~ Age + I(Age^2) + Gender, random=~1|SubjID, data=dat,na.action="na.omit")
summary(m.beh)

#model for left
m.left = lme(SSRT_inverted ~ MRI_cort_area.ctx.lh.parsopercularis + Age + I(Age^2) + Gender + MRI_cort_area.ctx.total + DeviceSerialNumber, random=~1|SubjID, data=dat,na.action="na.omit")
summary(m.left)
# DeviceSerialNumber

#model for right
m.right = lme(SSRT_inverted ~ MRI_cort_area.ctx.rh.parsopercularis + Age + I(Age^2) + Gender + MRI_cort_area.ctx.total + DeviceSerialNumber, random=~1|SubjID, data=dat,na.action="na.omit")
summary(m.right)
# DeviceSerialNumber

#model for bilateral
m.bilateral = lme(SSRT_inverted ~ MRI_cort_area.ctx.bilateral.parsopercularis  + Age + I(Age^2) + Gender + MRI_cort_area.ctx.total + DeviceSerialNumber, random=~1|SubjID, data=dat,na.action="na.omit")
summary(m.bilateral)
# DeviceSerialNumber

m.total = lme(BEH_CANTAB_SST_SSRT_last_half ~ Age + I(Age^2) + Gender + MRI_cort_area.ctx.total + DeviceSerialNumber, random=~1|SubjID, data=dat,na.action="na.omit")
summary(m.total)

######################################################################
# bilateral pars operc thickness
m.bilthick = lme(SSRT_inverted ~ MRI_cort_thick.ctx.bilateral.parsopercularis  + Age + I(Age^2) + Gender + DeviceSerialNumber, random=~1|SubjID, data=dat,na.action="na.omit")
summary(m.bilthick)
# DeviceSerialNumber

######################################################################
## models predicting ROI with an interaction between age x gender x ROI

#model for left
m.leftint = lme(SSRT_inverted ~  I(Age^2) + Age*Gender*MRI_cort_area.ctx.lh.parsopercularis + MRI_cort_area.ctx.total + DeviceSerialNumber, random=~1|SubjID, data=dat,na.action="na.omit")
summary(m.leftint)

# just age x ROI
m.leftint = lme(SSRT_inverted ~  I(Age^2) + Gender + Age*MRI_cort_area.ctx.lh.parsopercularis + MRI_cort_area.ctx.total + DeviceSerialNumber, random=~1|SubjID, data=dat,na.action="na.omit")
summary(m.leftint)

#model for right
m.rightint = lme(SSRT_inverted ~ I(Age^2) + Age*Gender*MRI_cort_area.ctx.rh.parsopercularis + MRI_cort_area.ctx.total + DeviceSerialNumber, random=~1|SubjID, data=dat,na.action="na.omit")
summary(m.rightint)

#model for bilateral area
m.bilint = lme(SSRT_inverted ~ I(Age^2) + Age*Gender*MRI_cort_area.ctx.bilateral.parsopercularis + MRI_cort_area.ctx.total + DeviceSerialNumber, random=~1|SubjID, data=dat,na.action="na.omit")
summary(m.bilint)

#model for bilateral thickness
m.bilint2 = lme(SSRT_inverted ~ I(Age^2) + Age*Gender*MRI_cort_thick.ctx.bilateral.parsopercularis + DeviceSerialNumber, random=~1|SubjID, data=dat,na.action="na.omit")
summary(m.bilint2)
######################################################################
## predicting ROI from SSRT ## portal analog for vertex-wise maps

#model for left
m.leftroi = lme(MRI_cort_area.ctx.lh.parsopercularis ~ SSRT_inverted + Age + I(Age^2) + Gender + MRI_cort_area.ctx.total + DeviceSerialNumber, random=~1|SubjID, data=dat,na.action="na.omit")
summary(m.leftroi)

#model for right
m.rightroi = lme(MRI_cort_area.ctx.rh.parsopercularis ~ SSRT_inverted + Age + I(Age^2) + Gender + MRI_cort_area.ctx.total + DeviceSerialNumber, random=~1|SubjID, data=dat,na.action="na.omit")
summary(m.rightroi)

#model for bilateral
m.bilateralroi = lme(MRI_cort_area.ctx.bilateral.parsopercularis ~ SSRT_inverted + Age + I(Age^2) + Gender + MRI_cort_area.ctx.total + DeviceSerialNumber, random=~1|SubjID, data=dat,na.action="na.omit")
summary(m.bilateralroi)

############################################################
## GAMM ##
m.bilgam = gam(SSRT_inverted ~ MRI_cort_area.ctx.bilateral.parsopercularis + s(Age) + Gender + MRI_cort_area.ctx.total + DeviceSerialNumber, random=~1|SubjID, data=dat,na.action="na.omit")
summary(m.bilgam)

m.gam1 = gam(SSRT_inverted ~ MRI_cort_area.ctx.bilateral.parsopercularis + s(SubjID, Age, bs="re") + Gender + MRI_cort_area.ctx.total + DeviceSerialNumber, random=~1|SubjID, data=dat,na.action="na.omit")
summary(m.gam1)

# random intercept and slope
m.gam2 = gam(SSRT_inverted ~ MRI_cort_area.ctx.bilateral.parsopercularis + s(SubjID, bs="re") + s(SubjID, Age, bs="re") + Gender + MRI_cort_area.ctx.total + DeviceSerialNumber, random=~1|SubjID, data=dat,na.action="na.omit")
summary(m.gam2)