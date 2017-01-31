setwd("/Users/chasereuter/Downloads")
library("spls")

################## DATA ##################
dat <- read.csv("allTP_goodMRI_CANTAB_baseline.csv", header = TRUE)
dat$SSRT_inverted <- (dat$BEH_CANTAB_SST_SSRT_last_half * (-1))

dat$nvisit=NA
dat$maxvisit=NA
for(i in 1:length(unique(dat$SubjID))){
  dat_i=dat[dat$SubjID==unique(dat$SubjID)[i],]
  dat_i$nvisit=1:dim(dat_i)[1]
  dat_i$maxvisit=dim(dat_i)[1]
  dat[dat$SubjID==unique(dat$SubjID)[i],]=dat_i
}
#maxvisit with 2 or more visits
dat = dat[dat$maxvisit>1,]
maxvisit = dat$maxvisit
###############################################################
area <- dat[,grepl("MRI_cort_area.ctx",names(dat))]
area = area[,!grepl("fuzzy",names(area)) & !grepl("total",names(area))]
l <- area[,grepl("ctx.lh.",names(area))]
r <- area[,grepl("ctx.rh.",names(area))]
tot = l + r
split = strsplit(names(tot), "[.]")
names(tot) = unlist(lapply(split, function(x){paste(x[c(1,2,4)],collapse = '.')}))

#as a ratio of total
tot = tot/dat$MRI_cort_area.ctx.total
#include total & SSRT to find intercepts
tot$MRI_cort_area.ctx.total = dat$MRI_cort_area.ctx.total
tot$SSRT_inverted = dat$SSRT_inverted

length.subj = length(unique(dat$SubjID))
Age = dat$Age
#calculate intercepts for each person for each variable
ints = matrix(nrow = length.subj , ncol = ncol(tot)) #rows = number of subjects, columns = 1 b0 for every variable
for(jj in 1:length(tot)){ #each variable
  for(i in 1:length.subj){ #each subject
    which = dat$SubjID==unique(dat$SubjID)[i]
    y = tot[which,jj]
    age = Age[which]
    age.c = age-mean(age)
    tmp.lm = lm(y ~ age.c)
    ints[i,jj] = tmp.lm$coef[1]
  }
}
#calculate average age per person
age.mean = c()
for(i in 1:length.subj){
  which = dat$SubjID==unique(dat$SubjID)[i]
  age=Age[which]
  age.mean[i] = mean(age)
}
ints = as.data.frame(ints)
names(ints) = paste0(names(tot),".b0")
#take out SSRT_inverted.b0 from predictors matrix
SSRT_inverted.b0 = ints$SSRT_inverted.b0 
ints$SSRT_inverted.b0 = NULL
#include age.mean into predictors
ints$age.mean = age.mean

######################################
######         SPLS #1          ######
######################################
nfold = length.subj-1

cv = cv.spls(ints, SSRT_inverted.b0 , fold = nfold , eta = seq(0.20,0.99,0.01), K = c(1:10) )

f <- spls( ints, SSRT_inverted.b0, eta = cv$eta.opt, K = cv$K.opt )
print(f)
coef.f = coef(f)
ci.f <- ci.spls( f, plot.it=TRUE, plot.fix='y')#, plot.var=10 )
cis <- ci.f$cibeta
cis =cis[[1]]

### for R^2 ###
vars = rownames(cis)
lm.dat = cbind(SSRT_inverted.b0, ints[,vars])
m = lm(SSRT_inverted.b0 ~., data = lm.dat)
summary(m)

##########################################################
# #2: slopes as predictors + age + interactions with age #
##########################################################
#calculate slopes
slopes = matrix(nrow = length.subj , ncol = ncol(tot)) #rows = number of subjects, columns = 1 b0 for every variable
for(jj in 1:length(tot)){ #each variable
  for(i in 1:length.subj){ #each subject
    which = dat$SubjID==unique(dat$SubjID)[i]
    y = tot[which,jj]
    age = Age[which]
    age.c = age-mean(age)
    tmp.lm = lm(y ~ age.c)
    slopes[i,jj] = tmp.lm$coef[2]
  }
}
slopes = as.data.frame(slopes)
names(slopes) = paste0(names(tot),".b1")
#take out SSRT_inverted.b1 from predictors matrix
SSRT_inverted.b1 = slopes$SSRT_inverted.b1 
slopes$SSRT_inverted.b1 = NULL
#interactions
slopes.x.age = slopes * age.mean
names(slopes.x.age) = paste0(names(slopes),".x.age")
#include age.mean into predictors
slopes$age.mean = age.mean
#combine slopes and slopes.x.age
slopes.all = cbind(slopes,slopes.x.age)

####################################
######         SPLS #2        ######
####################################
nfold = length.subj-1
cv2 = cv.spls(slopes.all, SSRT_inverted.b1 , fold = nfold , eta = seq(0.20,0.99,0.01), K = c(1:10) )
f2 <- spls( slopes.all , SSRT_inverted.b1, eta = cv2$eta.opt, K = cv2$K.opt )
print(f2)
coef.f2 = coef(f2)
ci.f2 <- ci.spls( f2, plot.it=TRUE, plot.fix='y')
cis2 <- ci.f2$cibeta
cis2 =cis2[[1]]

### for R^2 ###
vars = rownames(cis2)
lm.dat = cbind(SSRT_inverted.b1, slopes.all[,vars])
m = lm(SSRT_inverted.b1 ~., data = lm.dat)
summary(m)


####################################
######         SPLS #3        ######
####################################

### ints, slopes, age, ints x age, slopes x age

names(ints)
names(slopes.all)
ints$age.mean = NULL
ints.x.age = ints*age.mean
names(ints.x.age) = paste0(names(ints),".x.age")
ints.all = cbind(ints, ints.x.age)

ints.slopes.all = cbind(ints.all,slopes.all)

cv = cv.spls(ints.slopes.all, SSRT_inverted.b0 , fold = nfold , eta = seq(0.20,0.99,0.01), K = c(1:10) )
f <- spls( ints.slopes.all , SSRT_inverted.b0, eta = cv$eta.opt, K = cv$K.opt )
print(f)
coef.f = coef(f)
ci.f <- ci.spls( f, plot.it=TRUE, plot.fix='y')
cis <- ci.f$cibeta
cis =cis[[1]]

### for R^2 ###
vars = rownames(cis)
lm.dat = cbind(SSRT_inverted.b0, ints.slopes.all[,vars])
m = lm(SSRT_inverted.b0 ~., data = lm.dat)
summary(m)



############################################################################################################


##SPLS #4
names(ints)
names(slopes.all)
slopes.no.age.int = slopes.all[,!grepl("age",names(slopes.all))]
ints.slopes.no.ageinteraction = cbind(ints,slopes.no.age.int)
names(ints.slopes.no.ageinteraction)

cv4 = cv.spls(ints.slopes.no.ageinteraction, SSRT_inverted.b0 , fold = nfold , eta = seq(0.20,0.99,0.01), K = c(1:10) )
f4 <- spls( ints.slopes.no.ageinteraction , SSRT_inverted.b0, eta = cv4$eta.opt, K = cv4$K.opt )
print(f4)
coef.f4 = coef(f4)
ci.f4 <- ci.spls( f4, plot.it=TRUE, plot.fix='y')
cis4 <- ci.f4$cibeta
cis4 =cis4[[1]]

### for R^2 ###
vars = rownames(cis4)
lm.dat = cbind(SSRT_inverted.b0, ints.slopes.no.ageinteraction[,vars])
m = lm(SSRT_inverted.b0 ~., data = lm.dat)
summary(m)


