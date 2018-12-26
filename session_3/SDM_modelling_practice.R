###################################
#MODELLING PRACTICE
###################################

#loading input data
load("input_data.Rdata")

#objects in there:
#presence - presences of Abies alba in Europe
#presence.pseudoabsence
#presence.background
#predictors - predictors in Europe
#predictors.past.21k - climate during the last glacial maximum (colder than today)
#predictors.past.130k - climate during the last interglacial (warmer than today)
#predictors.future.2080 - future climate change scenario

#load functions
source("SDM_functions.R")

#loading libraries (COMMENT THE ONES YOU DON'T NEED)
library(tree)
library(partykit)
library(randomForest)
library(kernlab)
library(nnet)
library(NeuralNetTools)
library(mgcv)
library(plotmo)
library(earth)
library(maxnet)
library(dismo)
library(raster)
library(vegan)
library(lm.beta)
library(glmulti)
library(gbm)
library(devtools)
library(RColorBrewer)
library(HH)



################################################################
################################################################
#REDUCING COLLINEARITY OF PREDICTORS
################################################################
################################################################

#to data frame
predictors.df<-as.data.frame(predictors)

#remove NA
predictors.df<-na.omit(predictors.df)

#correlation matrix
predictors.correlation<-cor(predictors.df)

#correlation matrix as distance matrix
predictors.dist<-as.dist(abs(predictors.correlation))

#cluster
predictors.cluster<-hclust(1-predictors.dist)

#plot
plot(predictors.cluster)

#select one predictor per cluster below the 0.5 line
selected.predictors<-c("bio5", "bio6", "bio14", "ndvi_range")

#df with those predictors
predictors.df2<-predictors.df[ , selected.predictors]

#VARIANCE INFLATION FACTOR (library HH)
vif.result<-vif(predictors.df2)
vif.result

#example on how to remove predictors. Rinse and repeat
predictors.df2$"VARIABLE TO REMOVE HERE"<-NULL
vif.result<-vif(predictors.df2)
vif.result

#removing predictor

#removing predictor

#....

#new list of selected predictors
selected.predictors<-names(vif.result)

#selected predictors to brick
predictors<-predictors[[selected.predictors]]

###############################################################################
###############################################################################
#LOOKING AT THE DATA
###############################################################################
###############################################################################
#TAKE THE CODE YOU NEED FROM SDM_example_code.R

###############################################################################
###############################################################################
#MODELLING
###############################################################################
###############################################################################
#TAKE THE CODE YOU NEED FROM SDM_example_code.R