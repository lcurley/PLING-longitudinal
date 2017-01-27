## Average covariates across all timepoints to create a portal map

setwd("/Users/laurenbutlercurley/Documents/R_codes/")
rm(list=ls())

dat <- read.csv("allTP_goodMRI_CANTAB.csv", header = TRUE)
dat$MRI_cort_area.ctx.bilateral.parsopercularis <- (dat$MRI_cort_area.ctx.lh.parsopercularis + dat$MRI_cort_area.ctx.rh.parsopercularis)
dat$SSRT_inverted <- (dat$BEH_CANTAB_SST_SSRT_last_half * (-1))
dat$Age2 <- (dat$Age)^2

# get dataframe of SubjID and VisitID
ids <- dat[,2:3]

# get list of unique SubjIDs, no repeats
subs <- dat[!duplicated(dat$SubjID),] 
write.csv(subs, file = "firstTP_SSRT.csv") 

# average scores and demographics by gender
ageG <- aggregate(dat$Age,by=list(name=dat$Gender),data=dat,FUN=mean)
SSRTG <- aggregate(dat$SSRT_inverted,by=list(name=dat$Gender),data=dat,FUN=mean)

### create averages across timepoints for portal maps
avgAge <- aggregate(dat$Age,by=list(name=dat$SubjID),data=dat,FUN=mean)
avgAge2 <- aggregate(dat$Age2,by=list(name=dat$SubjID),data=dat,FUN=mean)
avgBilParsArea <- aggregate(dat$MRI_cort_area.ctx.bilateral.parsopercularis,by=list(name=dat$SubjID),data=dat,FUN=mean)
avgArea <- aggregate(dat$MRI_cort_area.ctx.total,by=list(name=dat$SubjID),data=dat,FUN=mean)
avgSSRT <- aggregate(dat$SSRT_inverted,by=list(name=dat$SubjID),data=dat,FUN=mean)

# merge spreadsheets
ages <- merge(avgAge, avgAge2, by="name")
areas <- merge(avgBilParsArea, avgArea, by="name")
age_areas <- merge(ages, areas, by="name")
final <- merge(age_areas, avgSSRT, by="name")
colnames(final) = c("SubjID", "avgAge", "avgAge2", "avgBilParsArea", "avgArea", "avgSSRTinverted")

# merge final, averaged data with original list of SubjIDs and VisitIDs
# write to file
newport <- merge(ids, final, by="SubjID")
write.csv(newport, file="Average_timepoints_portal_upload_all_VisitIDs.csv")

# fill in column names and write to file
write.csv(final, file="Average_timepoints_portal_upload.csv")