setwd("/Users/laurenbutlercurley/Documents/R_codes/")
library(nlme) #LME
library(mgcv) #GAM
library(ggplot2) #Plot
rm(list=ls())
dat <- read.csv("allTP_goodMRI_CANTAB.csv", header = TRUE)
##### invert SSRT scores #####
dat$SSRT_inverted <- (dat$BEH_CANTAB_SST_SSRT_last_half * (-1))

# find 1st and 3rd quartile ages
age1stq <- quantile(dat$Age, 0.25)
age3rdq <-quantile(dat$Age, 0.75)

# create two variables centered around the 1st and 3rd quartile - can be amended to be centered on min/max age
dat$age_1st_quart_centered <- dat$Age - age1stq
dat$age_3rd_quart_centered <- dat$Age - age3rdq

################################################################################
# run two models - one using the younger/1st quartile-centered age, the second using the older/3rd-quartile centered # age
m.SSRT1stAge = lme(SSRT_inverted ~ I(age_1st_quart_centered^2) + age_1st_quart_centered*Gender*MRI_cort_area.ctx.lh.parsopercularis 
                   + MRI_cort_area.ctx.total, random=~1|SubjID, data=dat,na.action="na.omit")

m.SSRT3rdAge = lme(SSRT_inverted ~ I(age_3rd_quart_centered^2) + age_3rd_quart_centered*Gender*MRI_cort_area.ctx.lh.parsopercularis 
                   + MRI_cort_area.ctx.total, random=~1|SubjID, data=dat,na.action="na.omit")

# run same two models in reverse, predicting ROI from each scaled age term
m.ROI1stAge = lme(MRI_cort_area.ctx.lh.parsopercularis ~ I(age_1st_quart_centered^2) + age_1st_quart_centered*Gender*SSRT_inverted 
                  + MRI_cort_area.ctx.total, random=~1|SubjID, data=dat,na.action="na.omit")

m.ROI3rdAge = lme(MRI_cort_area.ctx.lh.parsopercularis ~ I(age_3rd_quart_centered^2) + age_3rd_quart_centered*Gender*SSRT_inverted 
                  + MRI_cort_area.ctx.total, random=~1|SubjID, data=dat,na.action="na.omit")

################################################################################
# plot fit lines from each model on the same plot
# i.e. trendline for "younger" centered age term (1st quartile centered) and "older" centered age term (3rd quartile centered) on the same plot with left pars opercularis surface area on the x-axis and SSRT_inverted on the y-axis

pts_ssrt <- ggplot(data = dat, aes(x = MRI_cort_area.ctx.lh.parsopercularis, y = SSRT_inverted, group = SubjID))
ggplot(data = dat, aes(x = MRI_cort_area.ctx.lh.parsopercularis, y = SSRT_inverted, group = SubjID))
################################################################################

# plotting two things together:
plot(6:25,rnorm(20),type="b",xlim=c(1,30),ylim=c(-2.5,2.5),col=2)
par(new=T)
plot(rnorm(30),type="b",axes=F,col=3)
par(new=F)

# if addinig data of the same range:
plot(rnorm(100),type="l",col=2)
lines(rnorm(100),col=3)

# new sample
upvar<-rnorm(10)+seq(1,1.9,by=0.1)
downvar<-rnorm(20)*5+19:10
par(mar=c(5,4,4,4))
plot(6:15,upvar,pch=1,col=3,xlim=c(1,20),xlab="Occasion",ylab="",main="Dual ordinate plot")
mtext("upvar",side=2,line=2,col=3)
abline(lm(upvar~I(1:10)),col=3)

par(new=T)
plot(downvar,axes=F,xlab="",ylab="",pch=2,col=4)
axis(side=4)
abline(lm(downvar~I(1:20)),col=4)
mtext("downvar",side=4,line=2,col=4)

