# SSRT paper figures
library(ggplot2)
library(nlme)
setwd("/Users/laurenbutlercurley/Documents/R_codes")

dat <- read.csv("allTP_goodMRI_SSRT_2016.csv", header = TRUE)
dat$SSRT_inverted = (dat$BEH_CANTAB_SST_SSRT_last_half * (-1))

goodMRI <- dat[!is.na(dat$MRI_cort_thick.ctx.lh.bankssts),]
goodMRI_SST <- goodMRI[!is.na(goodMRI$BEH_CANTAB_SST_SSRT_last_half),]

################################################################################
# create bilateral cortical ROIs
goodMRI_SST$MRI_cort_area.ctx.bilateral.parsopercularis <- (goodMRI_SST$MRI_cort_area.ctx.lh.parsopercularis + goodMRI_SST$MRI_cort_area.ctx.rh.parsopercularis)
goodMRI_SST$MRI_cort_thick.ctx.bilateral.parsopercularis <- (goodMRI_SST$MRI_cort_thick.ctx.lh.parsopercularis + goodMRI_SST$MRI_cort_thick.ctx.rh.parsopercularis)/2

tab =  table(goodMRI_SST$BEH_CANTAB_SST_instance, goodMRI_SST$SubjID)
ids = tab[1,][tab[1,] == 1]
allSST = goodMRI_SST[goodMRI_SST$SubjID %in% names(ids),]

#baseSST <- goodMRI_SST[(goodMRI_SST$BEH_CANTAB_SST_instance==1),]

################################################################################

# plot longitudinal SubjID info
newd=allSST
newd$SubjID <- factor(newd$SubjID, levels=newd$SubjID[order(newd$Age)])
p <- ggplot(newd, aes(x=Age, y=SubjID, color=factor(Gender))) + geom_point() + geom_line(aes(Group=SubjID))
p + theme_bw() + scale_color_discrete(name="Gender") + theme(axis.text.y=element_blank(), axis.ticks=element_blank(), axis.text=element_text(size=16), 
                       axis.title=element_text(size=18), legend.text=element_text(size=16), legend.title=element_text(size=16))

ggsave("SubjID.png", width=12, height=10, dpi=300)
#png("SubjID.png", width=240, height=480, res=120)
#dev.off
######################################################################

# plot behavior
dev.off()
pts_ssrt <- ggplot(data = dat, aes(x = Age, y = SSRT_inverted, group = SubjID))
pts_ssrt + geom_point() + geom_line()
pts_ssrt + geom_point() + geom_line() + stat_smooth(aes(group = 1)) + theme_bw() + theme(axis.text=element_text(size=16), axis.title=element_text(size=18))
ggsave("SSRT_inverted.png", width=14, height=10, dpi=300)

#ggplot(allSST, aes(x=Age, y=SSRT_inverted, color=Gender)) + geom_point(shape=1) + geom_line(aes(Group=SubjID))

# plotting ROI
pts_bil <- ggplot(data = dat, aes(x = Age, y = MRI_cort_area.ctx.bilateral.parsopercularis, group = SubjID))
pts_bil + geom_point() + geom_line() + stat_smooth(aes(group = 1)) + theme_bw() + theme(axis.text=element_text(size=16), axis.title=element_text(size=18)) 
     + labs(y="Bilateral Pars Opercularis Surface Area")
ggsave("Bil_ParsOperc_area.png", width=14, height=10, dpi=300)

################################################################################

data = dat[,c("Age" , "SubjID" , "Gender" , "SSRT_inverted", "MRI_cort_area.ctx.lh.parsopercularis", "MRI_cort_area.ctx.total")]

model = lme(SSRT_inverted ~ I(Age^2) + Age*Gender*MRI_cort_area.ctx.lh.parsopercularis + 
                 MRI_cort_area.ctx.total, random=~1|SubjID, data=data,na.action="na.omit")

data$gender=2*(as.numeric(data$Gender=="F")-.5)
data$y=data$SSRT_inverted
mean_y=mean(data$y)
data$y=data$y-mean(data$y)
data$age=data$Age
data$age=data$age-mean(data$age)
data$x=data$MRI_cort_area.ctx.lh.parsopercularis
data$x=data$x-mean(data$x)

#to residualize y(gadi_sum) ...to data$resid
#base.model=lme(y~I(age^2) + gender*age + MRI_cort_area.ctx.total, random = ~1|SubjID , data=data)
#data$resid=resid(base.model)+mean_y

#model: using y as dependent variable, not resid
model= lme(y ~ I(age^2) + age*gender*x + MRI_cort_area.ctx.total, random=~1|SubjID, data=data,na.action="na.omit")
summary(model)

MRI_cort_area.ctx.lh.parsopercularis=seq(min(data$x),max(data$x),5)
data_pred.hi=cbind.data.frame(
     x=MRI_cort_area.ctx.lh.parsopercularis,
     age=quantile(data$age,.75),
     MRI_cort_area.ctx.total=mean(data$MRI_cort_area.ctx.total),
     gender=0)     

data_pred.lo=cbind.data.frame(
     x=MRI_cort_area.ctx.lh.parsopercularis,
     age=quantile(data$age,.25),
     MRI_cort_area.ctx.total=mean(data$MRI_cort_area.ctx.total),
     gender=0)	
gam_pred.hi=predict(model,newdata=data_pred.hi, level = 0)
gam_pred.lo=predict(model,newdata=data_pred.lo, level = 0)

#par(mfrow=c(1,2))
#non-residualized
plot(data$MRI_cort_area.ctx.lh.parsopercularis , data$SSRT_inverted , pch=15 , cex=.5, type="n",
     xlab="Left Pars Opercularis Surface Area",ylab="SSRT inverted",col=2, ylim=c(-600,0))
lines(data$MRI_cort_area.ctx.lh.parsopercularis[data$age<median(data$age)],
      data$SSRT_inverted[data$age<median(data$age)],col=4,type="p",pch=15,cex=.5, type="n")
lines(MRI_cort_area.ctx.lh.parsopercularis+mean(data$MRI_cort_area.ctx.lh.parsopercularis),gam_pred.hi+mean_y,col=2,lwd=2)
lines(MRI_cort_area.ctx.lh.parsopercularis+mean(data$MRI_cort_area.ctx.lh.parsopercularis),gam_pred.lo+mean_y,col=4,lwd=2)
#legend(2000,-400,legend=c("Older","Younger"),col=c(2,4),lwd=2, cex=0.8)
legend('bottomright', inset=0.2,legend=c("Older","Younger"),col=c(2,4),lwd=2, cex=0.8)
dev.copy(png,'Fig4.png')
dev.off()