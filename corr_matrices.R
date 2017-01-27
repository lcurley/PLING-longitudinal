# Correlation matrices and plotting

setwd("/Users/laurenbutlercurley/Documents/R_codes/first TP/")

library(lattice)
dat <- read.csv("firstTP_DTI_FA.csv", header = TRUE)
dat <- read.csv("firstTP_MRI_area_ChiHua.csv", header = TRUE)
# select only right hemisphere
df1 <- dat[ , grepl("rh",names(dat))] # get right hemisphere ROIS
df2 <- df1[ , grepl("12",names(df1))] # keep only the fuzzy 12 clusters
correls <- cor(df2)
levelplot(correls)

library(corrplot)
corrplot(correls, method='circle') # number, shade color

#corrgram(x, order= TRUE, panel=, lower.panel=panel.shade, upper.panel=panel.ellipse, 
#          text.panel=, diag.panel=, mean = "Title of my chart")
          # x is a dataframe with 1 obs per row

library(ggplot2)
library(reshape)
correls<- cor(dat)
correls.m <- melt(correls) 
ggplot(correls.m, aes(X1, X2, fill = value)) + geom_tile() + 
    scale_fill_gradient2(low = "blue",  high = "yellow")


library(plotrix)
library(seriation)
library(MASS)
plotcor(cor(dat), mar=c(0.1, 4, 4, 0.1))