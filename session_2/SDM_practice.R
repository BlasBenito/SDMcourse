#################
#CALIBRATING SDMs
#################

setwd("PATH")

#we need a few new packages
install.packages("tree", dep=TRUE) #regression trees
install.packages("rJava", dep=TRUE) #run "sudo R CMD javareconf" in linux
install.packages("partykit", dep=TRUE) #conditional inference trees (requires rJava)
install.packages("randomForest", dep=TRUE) #random forest
install.packages("kernlab", dep=TRUE) #support vector machines
install.packages("nnet", dep=TRUE) #artificial neural networks
install.packages("mgcv", dep=TRUE) #GAM
install.packages("plotmo", dep=TRUE) #response curves
install.packages("earth", dep=TRUE) #MARS
install.packages("maxnet", dep=TRUE) #MAXENT
install.packages("vegan", dep=TRUE) #NMDS
install.packages("lm.beta", dep=TRUE) #standardized coefficients for GLM
install.packages("glmulti", dep=TRUE) #variable selection for glm
install.packages("gbm", dep=TRUE) #gradient boosting modelling (boosted regression trees)
install.packages("devtools", dep=TRUE) #to access functions stored in Github

#loading libraries
library(tree)
library(partykit)
library(randomForest)
library(kernlab)
library(nnet)
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


#load some functions
source("SDM_functions.R")

#load presences and predictors
load("./data/predictors_and_presence.Rdata")


################################################
#PRELIMINARY EXPLORATION OF PRESENCE DATA
################################################

#saving plot
#pdf("./results/exploratory_plots.pdf", width=15, height=6, pointsize=30)

#iterating through predictors
for (predictor in names(predictors)){

par(mfrow=c(1,3))

#límites de la predictor en los gráficos
predictor.min<-min(presence.background[, predictor])
predictor.max<-max(presence.background[, predictor])

#density plot
#------------
#presences vs. background
d0<-density(presence.background[presence.background$presence==0, predictor], from=predictor.min, to=predictor.max)
d1<-density(presence.background[presence.background$presence==1, predictor], from=predictor.min, to=predictor.max)
plot(d0, col="red", xlim=c(predictor.min,predictor.max), lwd=3, main=paste(predictor, " vs. background", sep=""), xlab=predictor, ylim=c(0, max(d0$y, d1$y)))
lines(d1, col="forestgreen", lwd=3)

#presences vs absences
d0<-density(presence.absence[presence.absence$presence==0, predictor], from=predictor.min, to=predictor.max)
d1<-density(presence.absence[presence.absence$presence==1, predictor], from=predictor.min, to=predictor.max)
plot(d0, col="red", xlim=c(predictor.min,predictor.max), lwd=3, main=paste(predictor, " vs. absences", sep=""), xlab=predictor, ylim=c(0, max(d0$y, d1$y)))
lines(d1, col="forestgreen", lwd=3)

#presences vs pseudoabsences
d0<-density(presence.pseudoabsence[presence.pseudoabsence$presence==0, predictor], from=predictor.min, to=predictor.max)
d1<-density(presence.pseudoabsence[presence.pseudoabsence$presence==1, predictor], from=predictor.min, to=predictor.max)
plot(d0, col="red", xlim=c(predictor.min,predictor.max), lwd=3, main=paste(predictor, " vs. pseudoabsences", sep=""), xlab=predictor, ylim=c(0, max(d0$y, d1$y)))
lines(d1, col="forestgreen", lwd=3)
}

#closing pdf
#dev.off()


####################################
#2D PLOTS
####################################
#creating a temporary table
#NOTE: WE USE BACKGROUND DATA BECAUSE IT REPRESENTS THE WHOLE RANGE OF AVAILABLE CONDITIONS
table.visualization<-presence.background
#reordering the table to have presences at the end (so they will be plotted on top)
table.visualization<-table.visualization[order(table.visualization$presence, decreasing=FALSE), ]
#colors
table.visualization$color<-"gray60"
table.visualization[which(table.visualization$presence==1), "color"]<-"red4" 
#different size for points
table.visualization$cex<-0.6
table.visualization[which(table.visualization$presence==1), "cex"]<-0.8

#let's see the data in 2D with scatterplots
par(mfrow=c(2,2), mar=c(4,4,2,2))
plot(table.visualization$bio5 ~ table.visualization$bio6, col=table.visualization$color, cex=table.visualization$cex, ylab="bio5", xlab="bio6", pch=20)
plot(table.visualization$bio5 ~ table.visualization$bio14, col=table.visualization$color, cex=table.visualization$cex, ylab="bio5", xlab="bio14", pch=20)
plot(table.visualization$landcover_veg_tree ~ table.visualization$ndvi_range, col=table.visualization$color, cex=table.visualization$cex, ylab="landcover_veg_tree", xlab="ndvi_range", pch=20)
plot(table.visualization$human_footprint ~ table.visualization$topo_diversity, col=table.visualization$color, cex=table.visualization$cex, ylab="human_footprint", xlab="topo_diversity", pch=20)


##########################################
#NON METRIC MULTIDIMENSIONAL SCALING
##########################################
#WE USE PSEUDOABSENCES BECAUSE NMDS TAKE FOREVER WHEN DONE WITH THOUSANDS OF POINTS
table.visualization<-presence.pseudoabsence
#reordering the table to have presences at the end (so they will be plotted on top)
table.visualization<-table.visualization[order(table.visualization$presence, decreasing=FALSE), ]
#colors
table.visualization$color<-"gray60"
table.visualization[which(table.visualization$presence==1), "color"]<-"red4" 
#different size for points
table.visualization$cex<-1
table.visualization[which(table.visualization$presence==1), "cex"]<-0.8

#NMDS
nmds<-metaMDS(table.visualization[,names(predictors)], distance="euclidean")
#rules to interpret stress:
#stress < 0.05: excellent representation
#stress < 0.1: good representation
#stress < 0.2: acceptable representation,
#stress > 0.3: unsatisfactory representation

#PLOT
# pdf("./results/nmds.pdf", width=10, height=7, pointsize=15)
x11()
plot(nmds, type="n", main="NMDS")
points(nmds$points, col=table.visualization$color, cex=table.visualization$cex, pch=19)
ordisurf(nmds, table.visualization$bio6, add=T, col="blue", lwd=2, cex=3)
ordisurf(nmds, table.visualization$bio5, add=T, col="forestgreen", lwd=2, cex=3)
legend("bottomright", c("pseudoabsence","presence"), pch=19, col=c("gray60","red4"))
legend("topright", title="Lines", c("Bio6", "Bio5"), lty=1, lwd=2, col=c("blue", "forestgreen"))
# dev.off()

#REMOVING OBJECTS WE DON'T NEED ANY LONGER
rm(d0, d1, nmds,table.visualization, predictor.max, predictor, predictor.min)
gc()



#################################################################
#################################################################
#CALIBRATING MODELS
#################################################################
#################################################################

#################################################################
#BIOCLIM
#################################################################
#fitting model (ONLY PRESENCE)
bioclim.model=bioclim(predictors, presence[, c("x", "y")])

#plot model (use names(predictors) as reference for the numbers in "a" and "b")
x11()
par(mfrow=c(2,2), mar=c(4,4,3,3))
plot(bioclim.model, a=1, b=2)
plot(bioclim.model, a=1, b=3)
plot(bioclim.model, a=1, b=4)
plot(bioclim.model, a=1, b=5)

#predict map and plot
bioclim.model.map=predict(predictors, bioclim.model)
x11()
plot(bioclim.model.map)


#################################################################
#DOMAIN
#################################################################
#fitting model (ONLY PRESENCE)
domain.model=domain(predictors, presence[, c("x", "y")])

#predict map and plot
domain.model.map=predict(x=predictors, domain.model) #FREAKING SLOW
x11()
plot(domain.model.map)

 

#################################################################
#GENERALIZED LINEAR MODELS (GLM)
#################################################################

#TWO PREDICTORS
#--------------
glm.model<-glm(presence ~ bio5 + ndvi_range, family=binomial(link=logit), data=presence.pseudoabsence)
#family=binomial(link=logit) -> we are doing LOGISTIC REGRESSION HERE

#summary
summary(glm.model)
#NOTE: estimates cannot be interpreted as effect size if predictors are not in the same scale

#we can compute standardized coefficients with lm.beta
lm.beta(glm.model)

#D-squared (explained deviance)
Dsquared(glm.model)

#response curve
plotmo(glm.model, level=0.68)

#prediction
glm.map<-predict(predictors, glm.model, type="response")

#plot
x11()
plot(glm.map)
points(presence[, c("x","y")], cex=0.1)


#AHORA CON DOS predictors QUE INTERACCIONAN
#-----------------------------------------
glm.model<-glm(presence ~ bio5 * ndvi_range, family=binomial(link=logit), data=presence.pseudoabsence)

summary(glm.model)
lm.beta(glm.model)
Dsquared(glm.model)
plotmo(glm.model, level=0.68, all2=TRUE)
glm.map<-predict(predictors, glm.model, type="response")
x11()
plot(glm.map)
points(presence[, c("x","y")], cex=0.1)


#TWO INTERACTING POLYNOMIAL PREDICTORS
#----------------------------------------------------
glm.model<-glm(presence ~ poly(bio5, 2) * poly(ndvi_range, 2), family=binomial(link=logit), data=presence.pseudoabsence)

summary(glm.model)
lm.beta(glm.model)
Dsquared(glm.model)
plotmo(glm.model, level=0.68, all2=TRUE)
glm.map<-predict(predictors, glm.model, type="response")
x11
plot(glm.map)
points(presence[, c("x","y")], cex=0.1)


#TWO INTERACTING POLYNOMIAL PREDICTORS (FOURTH DEGREE POLYNOMIAL)
#----------------------------------------------------
glm.model<-glm(presence ~ poly(bio5, 4) * poly(ndvi_range, 4), family=binomial(link=logit), data=presence.pseudoabsence)

summary(glm.model)
lm.beta(glm.model)
Dsquared(glm.model)
plotmo(glm.model, level=0.68, all2=TRUE)
glm.map<-predict(predictors, glm.model, type="response")
plot(glm.map)
points(presence[, c("x","y")], cex=0.1)



#WORKING WITH ALL PREDICTORS, SECOND DEGREE POLYNOMIAL, NO INTERACTIONS
#----------------------------------------
#formula
formula.regression.poly<-as.formula(paste("presence ~ poly(", paste(names(predictors), collapse=", 2) + poly("), ", 2)", collapse=""))
formula.regression.poly

#NOTE: we'll calibrate the models with every dataset (background, weighted background, absence and pseudoabsence)

#BACKGROUND
glm.model.background<-glm(formula.regression.poly, family=binomial(link=logit), data=presence.background)

#WEIGHTED BACKGROUND (check "weights")
glm.model.backgroundw<-glm(formula.regression.poly, family=quasibinomial(link=logit), data=presence.background, weights=WeightPresenceBackground(presence.background[ , "presence"]))

#ABSENCE
glm.model.absence<-glm(formula.regression.poly, family=binomial(link=logit), data=presence.absence)

#PSEUDOabsence
glm.model.pseudoabsence<-glm(formula.regression.poly, family=binomial(link=logit), data=presence.pseudoabsence)


#response curves
x11()
plotmo(glm.model.background, level=0.68)
x11()
plotmo(glm.model.backgroundw, level=0.68)
x11()
plotmo(glm.model.absence, level=0.68)
x11()
plotmo(glm.model.pseudoabsence, level=0.68)

#summary
summary(glm.model.background)
summary(glm.model.backgroundw)
summary(glm.model.absence)
summary(glm.model.pseudoabsence)

#lm.beta
lm.beta(glm.model.background)
lm.beta(glm.model.backgroundw)
lm.beta(glm.model.absence)
lm.beta(glm.model.pseudoabsence)

#explained deviance
Dsquared(glm.model.background)
Dsquared(glm.model.backgroundw)
Dsquared(glm.model.absence)
Dsquared(glm.model.pseudoabsence)

#map
glm.model.background.map<-predict(predictors, glm.model.background, type="response")
glm.model.backgroundw.map<-predict(predictors, glm.model.backgroundw, type="response")
glm.model.absence.map<-predict(predictors, glm.model.absence, type="response")
glm.model.pseudoabsence.map<-predict(predictors, glm.model.pseudoabsence, type="response")

#maps
x11()
par(mfrow=c(2,2))
plot(glm.model.background.map, main="background")
plot(glm.model.backgroundw.map, main="background ponderado")
plot(glm.model.absence.map, main="absence")
plot(glm.model.pseudoabsence.map, main="pseudo-absence")


#variable selection with glmulti
#-------------------------------
#first we fit a model
glm.model<-glm(formula.regression.poly, family=binomial(link="logit"), data=presence.pseudoabsence)

#model selection (many options here, read the help file!)
glm.model.selection<-glmulti(glm.model, level=1, method="g")

#result
summary(glm.model.selection)

#best model
formula.best.model=as.formula(summary(glm.model.selection)$bestmodel)
formula.best.model
#and from here you can just fit the model again with the new formula.


#################################################################
#GENERALIZED ADDITIVE MODELS
#################################################################
#formula
formula.gam<-as.formula(paste("presence ~ s(", paste(names(predictors), collapse=", k=4) + s("), ", k=4)", collapse=""))
formula.gam

#background
gam.model.background<-gam(formula.gam, family=binomial(link=logit), data=presence.background)

#weighted background (sloooooow)
gam.model.backgroundw<-gam(formula.gam, family=quasibinomial(link=logit), data=presence.background, weights=WeightPresenceBackground(presence.background[ , "presence"]))

#absence
gam.model.absence<-gam(formula.gam, family=binomial(link=logit), data=presence.absence)

#pseudoabsence
gam.model.pseudoabsence<-gam(formula.gam, family=binomial(link=logit), data=presence.pseudoabsence)

#response curves
x11()
plotmo(gam.model.background, level=0.68, type="response")
x11()
plotmo(gam.model.backgroundw, level=0.68, type="response")
x11()
plotmo(gam.model.absence, level=0.68, type="response")
x11()
plotmo(gam.model.pseudoabsence, level=0.68, type="response")

#summary
summary(gam.model.background)
summary(gam.model.backgroundw)
summary(gam.model.absence)
summary(gam.model.pseudoabsence)

#explained deviance
Dsquared(gam.model.background)
Dsquared(gam.model.backgroundw)
Dsquared(gam.model.absence)
Dsquared(gam.model.pseudoabsence)

#maps
gam.model.background.map<-predict(predictors, gam.model.background, type="response")
gam.model.backgroundw.map<-predict(predictors, gam.model.backgroundw, type="response")
gam.model.absence.map<-predict(predictors, gam.model.absence, type="response")
gam.model.pseudoabsence.map<-predict(predictors, gam.model.pseudoabsence, type="response")

#maps
x11()
par(mfrow=c(2,2))
plot(gam.model.background.map, main="background")
plot(gam.model.backgroundw.map, main="background ponderado")
plot(gam.model.absence.map, main="absence")
plot(gam.model.pseudoabsence.map, main="pseudo-absence")



#################################################################
#MULTIVARIATE ADAPTIVE REGRESSION SPLINES (MARS)
#################################################################
#we cannot use a formula here..
#-------------------------
formula.mars<-as.formula(paste("presence ~ ", paste(names(predictors), collapse="+")))
formula.mars

#background
mars.model.background<-earth(formula.mars, glm = list(family=binomial(link="logit")), data=presence.background)

#weighted background SLOOOOOOOOW
# mars.model.backgroundw<-earth(formula.mars, glm = list(family=quasibinomial(link="logit")), data=presence.background, weights=WeightPresenceBackground(presence.background[ , "presence"]))

#absence
mars.model.absence<-earth(formula.mars, glm = list(family=binomial(link="logit")), data=presence.absence)

#pseudoabsence
mars.model.pseudoabsence<-earth(formula.mars, glm = list(family=binomial(link="logit")), data=presence.pseudoabsence)

#response curves
x11()
plotmo(mars.model.background, type="response")
# x11()
# plotmo(mars.model.backgroundw, type="response")
x11()
plotmo(mars.model.absence, type="response")
x11()
plotmo(mars.model.pseudoabsence, letype="response")

#summary
summary(mars.model.background)
# summary(mars.model.backgroundw)
summary(mars.model.absence)
summary(mars.model.pseudoabsence)

#variable importance
mars.model.background.evimp<-evimp(mars.model.background)
# summary(mars.model.backgroundw)
mars.model.absence.evimp<-evimp(mars.model.absence)
mars.model.pseudoabsence.evimp<-evimp(mars.model.pseudoabsence)


#maps
mars.model.background.map<-predict(predictors, mars.model.background, type="response")
# mars.model.backgroundw.map<-predict(predictors, mars.model.backgroundw, type="response")
mars.model.absence.map<-predict(predictors, mars.model.absence, type="response")
mars.model.pseudoabsence.map<-predict(predictors, mars.model.pseudoabsence, type="response")

#maps
x11()
par(mfrow=c(2,2))
plot(mars.model.background.map, main="background")
# plot(mars.model.backgroundw.map, main="background ponderado")
plot(mars.model.absence.map, main="absence")
plot(mars.model.pseudoabsence.map, main="pseudo-absence")


#################################################################
#MAXENT
#################################################################

#fit maxent model (PRESENCE AND BACKGROUND ONLY)
maxent.model=maxnet(p=presence.background$presence, data=presence.background[, names(predictors)], regmult=6)

#what's inside the model?
str(maxent.model)

#coefficients of predictors, linear combinations of predictors, and particular ranges of predictors
sort(abs(maxent.model$betas), decreasing=TRUE)

#response curves
x11()
plot(maxent.model, type="logistic")

#map
maxent.model.map<-predict(predictors, maxent.model, type="logistic")
x11()
plot(maxent.model.map)


#################################################################
#REGRESSION TREES
#################################################################
#we need a formula like this
formula.mars

#to avoid confussions...
formula.tree=formula.mars

#fitting regression trees with the "tree" package
tree.model.background<-tree(formula.tree, data=presence.background)
tree.model.pseudoabsence<-tree(formula.tree, data=presence.pseudoabsence)
tree.model.absence<-tree(formula.tree, data=presence.absence)

#plotting trees
x11()
par(mfrow=c(1,3))
plot(tree.model.background)
text(tree.model.background, all=TRUE, cex=1)
plot(tree.model.pseudoabsence)
text(tree.model.pseudoabsence, all=TRUE, cex=1)
plot(tree.model.absence)
text(tree.model.absence, all=TRUE, cex=1)

#prediction to map
tree.model.background.map<-predict(predictors, tree.model.background)
tree.model.pseudoabsence.map<-predict(predictors, tree.model.pseudoabsence)
tree.model.absence.map<-predict(predictors, tree.model.absence)

#maps
x11()
par(mfrow=c(1,3), oma=c(0,0,0,0))
plot(tree.model.background.map)
plot(tree.model.pseudoabsence.map)
plot(tree.model.absence.map)


#################################################################
#CONDITIONAL INFERENCE TREES
#################################################################
#more robust than "tree"
#we can use the same formula (formula.tree)

#fitting trees
ctree.model.background<-ctree(formula.tree, data=presence.background)
ctree.model.pseudoabsence<-ctree(formula.tree, data=presence.pseudoabsence)
ctree.model.absence<-ctree(formula.tree, data=presence.absence)

#plots (they can get quite complex with many variables)
plot(ctree.model.background)
x11()
plot(ctree.model.pseudoabsence)
x11()
plot(ctree.model.absence)

#predict
ctree.model.background.map<-predict(predictors, ctree.model.background, type="response")
ctree.model.pseudoabsence.map<-predict(predictors, ctree.model.pseudoabsence, type="response")
ctree.model.absence.map<-predict(predictors, ctree.model.absence, type="response")

#maps
x11()
par(mfrow=c(1,3))
plot(ctree.model.background.map)
plot(ctree.model.pseudoabsence.map)
plot(ctree.model.absence.map)


#################################################################
#RANDOM FOREST
#################################################################
#we can use the same "formula.tree"

#BACKGROUND (a bit slow)
rf.model.background<-randomForest(formula.tree, data=presence.background, importance=TRUE, ntree=500, mtry=3, nodesize=10)

#WEIGHTED BACKGROUND (CANNOT BE DONE IN RANDOM FOREST)

#absence
rf.model.absence<-randomForest(formula.tree, data=presence.absence, importance=TRUE, ntree=500, mtry=3, nodesize=10)

#PSEUDO-absence
rf.model.pseudoabsence<-randomForest(formula.tree, data=presence.pseudoabsence, importance=TRUE, ntree=500, mtry=3, nodesize=10)

#ERROR PLOT
x11()
par(mfrow=c(3, 1))
plot(rf.model.background)
plot(rf.model.absence)
plot(rf.model.pseudoabsence)

#model summary
print(rf.model.background)
print(rf.model.absence)
print(rf.model.pseudoabsence)

#response curves
x11()
plotmo(rf.model.background)
x11()
plotmo(rf.model.absence)
x11()
plotmo(rf.model.pseudoabsence)

#variable importance
varImpPlot(rf.model.background)
varImpPlot(rf.model.absence)
varImpPlot(rf.model.pseudoabsence)

#as a table
importance(rf.model.background)
importance(rf.model.absence)
importance(rf.model.pseudoabsence)

#to map
rf.model.background.map<-predict(predictors, rf.model.background, type="response")
rf.model.absence.map<-predict(predictors, rf.model.absence, type="response")
rf.model.pseudoabsence.map<-predict(predictors, rf.model.pseudoabsence, type="response")

#map
x11()
par(mfrow=c(1,3))
plot(rf.model.background.map, main="background")
plot(rf.model.absence.map, main="absence")
plot(rf.model.pseudoabsence.map, main="pseudoabsence")


#################################################################
#BOOSTED REGRESSION TREES
#################################################################
#CHECK https://cran.r-project.org/web/packages/dismo/vignettes/brt.pdf

#fitting models
#a realistic learning.rate would be 0.005, but that's quite slow

#background
gbm.model.background=gbm.step(data=presence.background, gbm.x = names(predictors), gbm.y = "presence", family = "bernoulli", tree.complexity = 5, learning.rate = 0.01, bag.fraction = 0.5)

#weighted background
gbm.model.backgroundw=gbm.step(data=presence.background, gbm.x = names(predictors), gbm.y = "presence", family = "bernoulli", tree.complexity = 5, learning.rate = 0.01, bag.fraction = 0.5, site.weights=WeightPresenceBackground(presence.background[ , "presence"]))

#pseudoabsencce
gbm.model.pseudoabsence=gbm.step(data=presence.pseudoabsence, gbm.x = names(predictors), gbm.y = "presence", family = "bernoulli", tree.complexity = 5, learning.rate = 0.01, bag.fraction = 0.5)

#absence
gbm.model.absence=gbm.step(data=presence.absence, gbm.x = names(predictors), gbm.y = "presence", family = "bernoulli", tree.complexity = 5, learning.rate = 0.01, bag.fraction = 0.5)

#plot error
plot(gbm.model.background)
plot(gbm.model.pseudoabsence)
plot(gbm.model.absence)

#variable importance
summary(gbm.model.background)
summary(gbm.model.pseudoabsence)
summary(gbm.model.absence)

#plot
gbm.plot(gbm.model.background)
gbm.plot(gbm.model.backgroundw)
gbm.plot(gbm.model.pseudoabsence)
gbm.plot(gbm.model.absence)

#also with plotmo
plotmo(gbm.model.background)
plotmo(gbm.model.backgroundw)
plotmo(gbm.model.pseudoabsence)
plotmo(gbm.model.absence)

#prediction (check the ntrees argument)
gbm.model.background.map<-predict(predictors, gbm.model.background, n.trees=gbm.model.background$gbm.call$best.trees, type="response")

gbm.model.backgroundw.map<-predict(predictors, gbm.model.backgroundw, n.trees=gbm.model.backgroundw$gbm.call$best.trees, type="response")

gbm.model.absence.map<-predict(predictors, gbm.model.absence, n.trees=gbm.model.absence$gbm.call$best.trees, type="response")

gbm.model.pseudoabsence.map<-predict(predictors, gbm.model.pseudoabsence, n.trees=gbm.model.pseudoabsence$gbm.call$best.trees, type="response")

#map
x11()
par(mfrow=c(2,2))
plot(gbm.model.background.map, main="background")
plot(gbm.model.backgroundw.map, main="backgroundw")
plot(gbm.model.absence.map, main="absence")
plot(gbm.model.pseudoabsence.map, main="pseudoabsence")

#CHECKING INTERACTIONS
#analyzing interactions
gbm.model.pseudoabsence.interactions <- gbm.interactions(gbm.model.pseudoabsence)
#interaction matrix
gbm.model.pseudoabsence.interactions$interactions
#relevant interactions
gbm.model.pseudoabsence.interactions$rank.list

#plotting them specifically
par(mfrow=c(1,3))
gbm.perspec(gbm.model.absence, 1, 8)
gbm.perspec(gbm.model.pseudoabsence, 1, 6)
gbm.perspec(gbm.model.pseudoabsence, 1, 5)



#################################################################
#ARTIFICIAL NEURAL NETWORKS
#################################################################

#recommended links:
#https://stevenmiller888.github.io/mind-how-to-build-a-neural-network/
#http://stevenmiller888.github.io/mind-how-to-build-a-neural-network-part-2/
#http://neuralnetworksanddeeplearning.com/chap4.html


#number of hidden neurons (no sharp rule for this, but nnet requires more hidden neurons than input predictors to work well)
hidden.neurons=20

#background (slower, reduce maxit)
nnet.model.background<-nnet(formula.tree, size=hidden.neurons, maxit=100, data=presence.background)

#weighted background (slower, reduce maxit)
nnet.model.backgroundw<-nnet(formula.tree, size=hidden.neurons, maxit=100, data=presence.background, weights=WeightPresenceBackground(presence.background[ , "presence"]))

#pseudoabsence
nnet.model.pseudoabsence<-nnet(formula.tree, size=hidden.neurons, maxit=10000, data=presence.pseudoabsence)

#absence
nnet.model.absence<-nnet(formula.tree, size=hidden.neurons, maxit=10000, data=presence.absence)

#installing a new package
install.packages("NeuralNetTools", dep=TRUE)
library(NeuralNetTools)

#plotting neural networks
plotnet(nnet.model.background)
plotnet(nnet.model.backgroundw)
plotnet(nnet.model.pseudoabsence)
plotnet(nnet.model.absence)

#variable importance
garson(nnet.model.background)
garson(nnet.model.backgroundw)
garson(nnet.model.pseudoabsence)
garson(nnet.model.absence)

#response curves
plotmo(nnet.model.background, all2=TRUE)
plotmo(nnet.model.backgroundw)
plotmo(nnet.model.pseudoabsence, all2=TRUE)
plotmo(nnet.model.absence)

#to maps
nnet.model.background.map<-predict(predictors, nnet.model.background)
nnet.model.backgroundw.map<-predict(predictors, nnet.model.backgroundw)
nnet.model.pseudoabsence.map<-predict(predictors, nnet.model.pseudoabsence)
nnet.model.absence.map<-predict(predictors, nnet.model.absence)

x11()
par(mfrow=c(2,2), oma=c(0,0,0,0))
plot(nnet.model.background.map, main="background")
plot(nnet.model.backgroundw.map, main="weighted background")
plot(nnet.model.pseudoabsence.map, main="pseudoabsence")
plot(nnet.model.absence.map, main="absence")

################################################################
#SUPPORT VECTOR MACHINES
#################################################################
help(ksvm)
# kernlab
#check out the "kernel" argument

#lets try different kernels with the same dataset
#kernel lineal
svm.model.pseudoabsence.vanilladot<-ksvm(formula.tree, data=presence.pseudoabsence, kernel="vanilladot")
svm.model.pseudoabsence.vanilladot

#kernel radial
svm.model.pseudoabsence.rbfdot<-ksvm(formula.tree, data=presence.pseudoabsence, kernel="rbfdot", kpar=list(sigma=1))
svm.model.pseudoabsence.rbfdot

#kernel laplace
svm.model.pseudoabsence.laplace<-ksvm(formula.tree, data=presence.pseudoabsence, kernel="laplacedot", kpar=list(sigma=1))
svm.model.pseudoabsence.laplace

#predict
svm.model.pseudoabsence.vanilladot.map<-predict(predictors, svm.model.pseudoabsence.vanilladot, type="response")
svm.model.pseudoabsence.rbfdot.map<-predict(predictors, svm.model.pseudoabsence.rbfdot, type="response")
svm.model.pseudoabsence.laplace.map<-predict(predictors, svm.model.pseudoabsence.laplace, type="response")

#plot
x11()
par(mfrow=c(1,3))
plot(svm.model.pseudoabsence.vanilladot.map, main="vanilladot")
plot(svm.model.pseudoabsence.rbfdot.map, main="rbfdot")
plot(svm.model.pseudoabsence.laplace.map, main="laplacedot")


#LET'S CHECK THE EFFECT OF THE DIFFERENT PRESENCE DATASETS
#background
svm.model.background<-ksvm(formula.tree, data=presence.background, kernel="laplacedot", kpar=list(sigma=1))
svm.model.background #check the training error!

#weighted background
svm.model.backgroundw<-ksvm(formula.tree, data=presence.background, kernel="laplacedot", kpar=list(sigma=1), class.weigths=WeightPresenceBackground(presence.background[ , "presence"]))
svm.model.backgroundw #check the training error!

#pseudoabsence
svm.model.pseudoabsence<-ksvm(formula.tree, data=presence.pseudoabsence, kernel="laplacedot", kpar=list(sigma=1))
svm.model.pseudoabsence

#absence
svm.model.absence<-ksvm(formula.tree, data=presence.absence, kernel="laplacedot", kpar=list(sigma=1))
svm.model.absence


#PREDICT
svm.model.background.map<-predict(predictors, svm.model.background, type="response")
svm.model.backgroundw.map<-predict(predictors, svm.model.backgroundw, type="response")
svm.model.pseudoabsence.map<-predict(predictors, svm.model.pseudoabsence, type="response")
svm.model.absence.map<-predict(predictors, svm.model.absence, type="response")

#map
x11()
par(mfrow=c(2,2))
plot(svm.model.background.map, main="background")
plot(svm.model.backgroundw.map, main="weighted background")
plot(svm.model.absence.map, main="absence")
plot(svm.model.pseudoabsence.map, main="pseudoabsence")



################################################################
################################################################
#ENSEMBLE MODEL FORECASTING
################################################################
################################################################

#BRICKS WITH ALL MODELS (by presence class + maxent)
models.background<-brick(rf.model.background.map, gbm.model.background.map, nnet.model.background.map, glm.model.background.map, mars.model.background.map, maxent.model.map, svm.model.background.map)
names(models.background)<-c("rf", "gbm", "nnet", "glm", "mars", "maxent", "svm")

models.backgroundw<-brick(gbm.model.backgroundw.map, nnet.model.backgroundw.map, glm.model.backgroundw.map, maxent.model.map, svm.model.backgroundw.map)
names(models.backgroundw)<-c("gbm", "nnet", "glm", "maxent", "svm")

models.pseudoabsence<-brick(gam.model.pseudoabsence.map, rf.model.pseudoabsence.map, gbm.model.pseudoabsence.map, nnet.model.pseudoabsence.map, glm.model.pseudoabsence.map, mars.model.pseudoabsence.map, maxent.model.map, svm.model.pseudoabsence.map)
names(models.pseudoabsence)<-c("gam","rf", "gbm", "nnet", "glm", "mars", "maxent", "svm")

models.absence<-brick(rf.model.absence.map, gbm.model.absence.map, rf.model.absence.map, nnet.model.absence.map, glm.model.absence.map, mars.model.absence.map, maxent.model.map, svm.model.absence.map)
names(models.absence)<-c("rf", "gbm", "rf", "nnet", "glm", "mars", "maxent", "svm")


#MEAN
ensemble.background.mean<-calc(models.background, mean)
ensemble.backgroundw.mean<-calc(models.backgroundw, mean)
ensemble.pseudoabsence.mean<-calc(models.pseudoabsence, mean)
ensemble.absence.mean<-calc(models.absence, mean)

#STANDARD DEVIATION
ensemble.background.sd<-calc(models.background, sd)
ensemble.backgroundw.sd<-calc(models.backgroundw, sd)
ensemble.pseudoabsence.sd<-calc(models.pseudoabsence, sd)
ensemble.absence.sd<-calc(models.absence, sd)

#PLOT
par(mfrow=c(2,4))
plot(ensemble.background.mean, main="background")
plot(ensemble.backgroundw.mean, main="backgroundw")
plot(ensemble.pseudoabsence.mean, main="pseudoabsence")
plot(ensemble.absence.mean, main="absence")
plot(ensemble.background.sd)
plot(ensemble.backgroundw.sd)
plot(ensemble.pseudoabsence.sd)
plot(ensemble.absence.sd)

#
x11()
par(mfrow=c(1,2))
plot(ensemble.pseudoabsence.mean, main="mean")
plot(ensemble.pseudoabsence.sd, main="sd")



#ANALIZING AN ENSEMBLE IN DETAIL
########################################
#brick to dataframe
models.pseudoabsence.df<-na.omit(as.data.frame(models.pseudoabsence))
#correlation matrix
models.pseudoabsence.cor<-cor(models.pseudoabsence.df)
#cor to dist
models.pseudoabsence.dis<-abs(as.dist(models.pseudoabsence.cor))
#plot
x11()
plot(hclust(1-models.pseudoabsence.dis))


#borra objetos que no vamos a usar más
rm(modelos, modelos.df, modelos.cor, modelos.dis)
gc()



##################################
##################################
#MODEL EVALUATION
##################################
##################################

#LET'S EVALUATE A SINGLE MODEL FIRST, JUST TO SEE HOW IT WORKS
#wE WILL DO WHAT IS CALLED "INTRINSIC EVALUATION", OR "EVALUATING WITH THE SAME PRESENCES WE USED TO TRAIN THE MODEL". THIS IS NOT OPTIMAL, BUT A GOOD STARTING POINT TO UNDERSTAND THE EVALUATION PROCESS

#let's choose a model and its map
model.to.evaluate<-rf.model.pseudoabsence
model.to.evaluate.map<-rf.model.pseudoabsence.map

#we need the values of the presences and pseudoabsences over the model
#it can be done in two different ways

#1. we extract the suitability values of the presences and pseudoabsences from the model's map
suitability.values.1<-extract(model.to.evaluate.map, presence.pseudoabsence[, c("x","y")])

#2. we use the model itself to predict suitability values for the presences and pseudoabsences
suitability.values.2<-predict(model.to.evaluate, newdata=presence.pseudoabsence, type="response")

#they are the same
x11()
plot(suitability.values.1, suitability.values.2)

#we separate the suitability values of the presences from the suitability values of the pseudoabsences
values.presences<-suitability.values.1[presence.pseudoabsence$presence==1]
values.pseudoabsences<-suitability.values.1[presence.pseudoabsence$presence==0]

#aplicamos la función evaluate de la librería dismo
evaluation<-evaluate(p=values.presences, a=values.pseudoabsences)
evaluation

#plots
x11()
par(mfrow=c(1,3), mar=c(2,2,4,2), oma=c(3,3,5,3))
density(evaluation)
boxplot(evaluation, col=c("blue", "red"))
plot(evaluation, "ROC")

#structure
str(evaluation)

#getting the auc value only
auc<-evaluation@auc
auc


#LET'S DO THIS FOR ALL MODELS
#############################
#we have a brick with every pseudoabsence model named "models.pseudoabsence"

auc.values<-vector()

#let's start the loop
for (model.name in names(models.pseudoabsence)){
  
  #select the current model
  model.map=models.pseudoabsence[[model.name]]
  
  #extract values of presences
  presence.values=extract(model.map, presence.pseudoabsence[presence.pseudoabsence$presence==1, c("x","y")])
  
  #extract absence values
  absence.values=extract(model.map, presence.pseudoabsence[presence.pseudoabsence$presence==0, c("x","y")])
  
  #evaluate and save result
  auc.values=c(auc.values, evaluate(p=presence.values, a=absence.values)@auc)
  
}

#putting the results together
auc.values<-data.frame(model=names(models.pseudoabsence), auc=auc.values)
auc.values<-auc.values[order(auc.values$auc, decreasing=TRUE), ]
auc.values


#USING KFOLD TO EVALUATE A MODEL
#################################################################################
#selecting the presence data
presence.data<-presence.pseudoabsence

#generating 10 groups
k<-10
groups<-kfold(presence.pseudoabsence, k)

#objects to store results
auc.kfold<-vector()
maps.kfold<-brick()

#model formula (change it depending on the model)
model.formula<-as.formula(paste("presence ~ ", paste(names(predictors), collapse="+")))

#for every group
for (group in 1:k) {
  
  #select all groups but the current one as training data
  training.data<-presence.data[groups != group,]
  
  #select the group data to test the model
  testing.data<-presence.data[groups == group,]
  
  #fit model here (any model can be here)
  trained.model<-randomForest(formula.tree, data=training.data, ntree=500, mtry=3, nodesize=10)
  
  #predict to map and save in the maps brick
  maps.kfold<-addLayer(maps.kfold, predict(predictors, trained.model))
  
  #evaluation
  auc.kfold[group]<-evaluate(p=predict(trained.model, testing.data[testing.data$presence==1, ]), a=predict(trained.model, testing.data[testing.data$presence==0, ]))@auc
  
} #end of loop

#summarizing results
auc.kfold.mean<-mean(auc.kfold)
auc.kfold.mean
auc.kfold.sd<-sd(auc.kfold)
auc.kfold.sd
maps.kfold.mean<-calc(maps.kfold, mean)
maps.kfold.sd<-calc(maps.kfold, sd)

#plotting maps
x11()
par(mfrow=c(1,2))
plot(maps.kfold.mean, main="Mean")
plot(maps.kfold.sd, main="Standard deviation")
 

trained.model<-randomForest(formula.tree, data=presence.data, ntree=500, mtry=3, nodesize=10)

#USING BOOTSTRAP TO EVALUATE A MODEL
####################################################################################

#computing number of rows to get as evaluation data
presence.data<-presence.pseudoabsence
presence.data.nrow=nrow(presence.data)
evaluation.data.percentage<-40
evaluation.data.rows<-round((evaluation.data.percentage*presence.data.nrow)/100, 0)

#number of iterations
iterations=10

#object to save auc values
auc.bootstrap<-vector()
maps.bootstrap<-brick()

#iteramos por los grupos
for (iteration in 1:iterations){
  
  print(paste("Iteration:", iteration))
  
  #random sampling of the presence data
  sampling<-sample(x=presence.data.nrow, size=evaluation.data.rows, replace=TRUE)
  
  #select all groups but the current one as training data
  training.data<-presence.data[-sampling, ]
  
  #select the group data to test the model
  testing.data<-presence.data[sampling, ]
  
  #fit model here (any model can be here)
  trained.model<-randomForest(formula.tree, data=training.data, ntree=500, mtry=3, nodesize=10)
  
  #predict to map and save in the maps brick
  maps.bootstrap<-addLayer(maps.bootstrap, predict(predictors, trained.model))
  
  #evaluation
  auc.bootstrap[iteration]<-evaluate(p=predict(trained.model, testing.data[testing.data$presence==1, ]), a=predict(trained.model, testing.data[testing.data$presence==0, ]))@auc
  
} #end of loop

#summarizing results
auc.bootstrap.mean<-mean(auc.bootstrap)
auc.bootstrap.mean
auc.bootstrap.sd<-sd(auc.bootstrap)
auc.bootstrap.sd
maps.bootstrap.mean<-calc(maps.bootstrap, mean)
maps.bootstrap.sd<-calc(maps.bootstrap, sd)

#plotting maps
par(mfrow=c(1,2))
plot(maps.bootstrap.mean, main="Mean")
plot(maps.bootstrap.sd, main="Standard deviation")

trained.model<-randomForest(formula.tree, data=presence.data, ntree=500, mtry=3, nodesize=10)

####################################################################################
####################################################################################
#THRESHOLDS
####################################################################################
####################################################################################

#we use the results of "evaluate" to compute thresholds
model.evaluation<-evaluate(p=predict(rf.model.pseudoabsence, presence.pseudoabsence[presence.pseudoabsence$presence==1, ]), a=predict(rf.model.pseudoabsence, presence.pseudoabsence[presence.pseudoabsence$presence==0, ]))
model.evaluation
str(model.evaluation)
#max TPR+TNR at : that is already a threshold maximizing the sum of correct classifications for presences and absences
#TPR: TRUE POSITIVE RATE
#TNR: TRUE NEGATIVE RATE

#we can compute it ourselves
threshold.TPR.TNR<-model.evaluation@t[which.max(model.evaluation@TPR + model.evaluation@TNR)]
threshold.TPR.TNR


#we can maximize other measures like "kappa"
threshold.MAXKAPPA<-model.evaluation@t[which.max(model.evaluation@kappa)]
threshold.MAXKAPPA

#or the CORRECT CLASSIFICATION RATE
threshold.CCR<-model.evaluation@t[which.max(model.evaluation@CCR)]
threshold.CCR

#FOR MANY MODELS ALL THESE THRESHOLDS WILL USUALLY BE THE SAME

#PLOT
x11()
par(mfrow=c(1,3))
plot(rf.model.pseudoabsence.map > threshold.TPR.TNR, main="TPR + TNR")
plot(rf.model.pseudoabsence.map > threshold.CCR, main="CCR")
plot(rf.model.pseudoabsence.map > threshold.MAXKAPPA, main="Maximum kappa")


######################################################################################
######################################################################################
#PROJECTION IN SPACE
######################################################################################
######################################################################################

#CHANGE NAME OF THE PREDICTORS BRICK
predictors.eu<-predictors

#IMPORT AND PREPARE PREDICTORS FROM NORTH AMERICA
list.predictors.na<-list.files(path="./data/spatial_projection/",pattern='*.asc', full.names=TRUE)
#brick
predictors.na<-brick(stack(list.predictors.na))
#keeping same predictors as in europe
predictors.na<-predictors.na[[names(predictors.eu)]]
#plot
plot(predictors.na)

#MESS ANALYSIS
mess.na<-mess(x=predictors.na, v=presence.pseudoabsence[, names(predictors.eu)], full=TRUE)
#naming and plotting result
names(mess.na)=c(names(predictors.na), "mess")
x11()
plot(mess.na)


#MAX DEVIATED PREDICTOR PER CELL
mess.na.max<-which.max(mess.na)
plot(mess.na.max)

#making it a bit more easy to interpret
library(RColorBrewer)
paleta.color<-brewer.pal(n=8, name="Set1")
#plot
x11()
plot(mess.na.max, col=paleta.color, legend=FALSE)
legend("bottomleft", c(names(predictors.na)), pch=15, col=paleta.color, cex=1.6)


################################
#MAKING A PROJECTION
################################
projection.eu<-predict(predictors.eu, rf.model.pseudoabsence, type="response")
projection.na<-predict(predictors.na, rf.model.pseudoabsence, type="response")

#plot
x11()
par(mfrow=c(1,2))
plot(projection.eu)
plot(projection.na)


#importing presence data from North America
presence.na<-read.table("./data/presence/occurrence.txt",header=T, sep='\t', fill=TRUE, check.names=TRUE, stringsAsFactors=FALSE)
str(presence.na)

#cropping presences to the study area (North America this time)
presences.over.predictors<-data.frame(extract(x=predictors.na, y=presence.na[ , c("longitude","latitude")]))
presence.na<-cbind(presence.na, presences.over.predictors)
presence.na<-presence.na[!(is.na(presence.na[ , names(predictors.na)])), ]

#removing duplicated coordinates
duplicated.xy<-duplicated(presence.na[ , c("latitude", "longitude")])
presence.na<-presence.na[!duplicated.xy, ]

#getting only coordinates
presence.na<-presence.na[ , c("longitude","latitude")]
names(presence.na)<-c("x", "y")

#plot
x11()
plot(projection.na)
points(presence.na, cex=0.1)


