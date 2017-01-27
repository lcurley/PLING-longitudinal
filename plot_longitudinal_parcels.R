######## loop over headers and plot with ggplot ##########
library(ggplot2)
df <- read.csv("usercache_PLING_lauren.csv", header=TRUE)

# find everything with "MRI_cort_area.ctx.*" in the header
area <- df[,grepl("MRI_cort_area.ctx.",names(df))]
thick <- df[,grepl("MRI_cort_thick.ctx.",names(df))]

#add back Age, Gender, VisitID, SubjID
dem <- subset(df, select=c("Age", "Age.at.Imaging", "Gender", "VisitID", "SubjID"))
dem_area <- cbind(dem, area)
dem_thick <- cbind(dem, thick)
# maybe remove all Age > 13?
dem_area_u13 <- dem_area[!(dem_area$Age.at.Imaging>13),]
dem_thick_u13 <- dem_thick[!(dem_thick$Age.at.Imaging>13),]

names <- colnames(thick)
lenn <- length(names)

setwd("/Users/laurenbutlercurley/Documents/R_codes/PLING_MRI_plots_longitudinal")
# plot using age.at.imaging

# ages 3-13
for (n in names){
     ggplot(data=dem_thick_u13, aes_string(x = "Age.at.Imaging", y = n, group = "SubjID")) + geom_point() + geom_line() + stat_smooth(aes(group = 1))
     ggsave(filename=paste("PLING_u13_",n,".png",sep=""))
}

# all ages 3-20
for (n in names){
     ggplot(data=dem_thick, aes_string(x = "Age.at.Imaging", y = n, group = "SubjID")) + geom_point() + geom_line() + stat_smooth(aes(group = 1))
     ggsave(filename=paste("PLING_",n,".png",sep=""))
}

############## plot DTI parcels ##############

dti_FA <- df[,grepl("DTI_fiber_FA",names(df))]
dti_MD <- df[,grepl("DTI_fiber_MD",names(df))]

dem <- subset(df, select=c("Age", "Age.at.Imaging", "Gender", "VisitID", "SubjID"))
dem_dtiFA <- cbind(dem, dti_FA)
dem_dtiMD <- cbind(dem, dti_MD)
# maybe remove all Age > 13?
dem_dtiFA_u13 <- dem_dtiFA[!(dem_dtiFA$Age.at.Imaging>13),]
dem_dtiMD_u13 <- dem_dtiMD[!(dem_dtiMD$Age.at.Imaging>13),]

names <- colnames(dti_MD)
lenn <- length(names)

setwd("/Users/laurenbutlercurley/Documents/R_codes/PLING_DTI_plots_longitudinal")
# plot using age.at.imaging

# ages 3-13
for (n in names){
     ggplot(data=dem_dtiMD_u13, aes_string(x = "Age.at.Imaging", y = n, group = "SubjID")) + geom_point() + geom_line() + stat_smooth(aes(group = 1))
     ggsave(filename=paste("PLING_u13_",n,".png",sep=""))
}

# all ages 3-20
for (n in names){
     ggplot(data=dem_dtiMD, aes_string(x = "Age.at.Imaging", y = n, group = "SubjID")) + geom_point() + geom_line() + stat_smooth(aes(group = 1))
     ggsave(filename=paste("PLING_",n,".png",sep=""))
}


#plot average size for cortical regions over time on the same plot
#(not total)
#maybe use % change from baseline to standardize differences in size? or proportion of total or both?

# loop and calculate % of total for each Desikan parcel
