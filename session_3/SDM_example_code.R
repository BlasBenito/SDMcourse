################################################
#PRELIMINARY EXPLORATION OF PRESENCE DATA
################################################

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
#GENERALIZED LINEAR MODELS (GLM)
#################################################################

#FORMULA SECOND DEGREE POLYNOMIAL, NO INTERACTIONS
#----------------------------------------
#formula
formula.regression.poly<-as.formula(paste("presence ~ poly(", paste(names(predictors), collapse=", 2) + poly("), ", 2)", collapse=""))
formula.regression.poly


#WEIGHTED BACKGROUND
glm.model.backgroundw<-glm(formula.regression.poly, family=quasibinomial(link=logit), data=presence.background, weights=WeightPresenceBackground(presence.background[ , "presence"]))

#PSEUDOabsence
glm.model.pseudoabsence<-glm(formula.regression.poly, family=binomial(link=logit), data=presence.pseudoabsence)

#response curves
plotmo(glm.model.backgroundw, level=0.68)
plotmo(glm.model.pseudoabsence, level=0.68)

#summary
summary(glm.model.backgroundw)
summary(glm.model.pseudoabsence)

#lm.beta
lm.beta(glm.model.backgroundw)
lm.beta(glm.model.pseudoabsence)

#explained deviance
Dsquared(glm.model.backgroundw)
Dsquared(glm.model.pseudoabsence)

#map
glm.model.backgroundw.map<-predict(predictors, glm.model.backgroundw, type="response")
glm.model.pseudoabsence.map<-predict(predictors, glm.model.pseudoabsence, type="response")

#maps
plot(glm.model.backgroundw.map, main="weighted background")
plot(glm.model.pseudoabsence.map, main="pseudo-absence")


#variable selection with glmulti
#-------------------------------
#model selection (many options here, read the help file!)
glm.model.selection<-glmulti(glm.model.pseudoabsence, level=1, method="g")

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


#weighted background (sloooooow)
gam.model.backgroundw<-gam(formula.gam, family=quasibinomial(link=logit), data=presence.background, weights=WeightPresenceBackground(presence.background[ , "presence"]))

#pseudoabsence
gam.model.pseudoabsence<-gam(formula.gam, family=binomial(link=logit), data=presence.pseudoabsence)

#response curves
plotmo(gam.model.backgroundw, level=0.68, type="response")
plotmo(gam.model.pseudoabsence, level=0.68, type="response")

#summary
summary(gam.model.backgroundw)
summary(gam.model.pseudoabsence)

#explained deviance
Dsquared(gam.model.backgroundw)
Dsquared(gam.model.pseudoabsence)

#maps
gam.model.backgroundw.map<-predict(predictors, gam.model.backgroundw, type="response")
gam.model.pseudoabsence.map<-predict(predictors, gam.model.pseudoabsence, type="response")

#maps
plot(gam.model.backgroundw.map, main="background ponderado")
plot(gam.model.pseudoabsence.map, main="pseudo-absence")



#################################################################
#MULTIVARIATE ADAPTIVE REGRESSION SPLINES (MARS)
#################################################################
formula.mars<-as.formula(paste("presence ~ ", paste(names(predictors), collapse="+")))
formula.mars

#weighted background SLOOOOOOOOW
mars.model.backgroundw<-earth(formula.mars, glm = list(family=quasibinomial(link="logit")), data=presence.background, weights=WeightPresenceBackground(presence.background[ , "presence"]))

#pseudoabsence
mars.model.pseudoabsence<-earth(formula.mars, glm = list(family=binomial(link="logit")), data=presence.pseudoabsence)

#response curves)
plotmo(mars.model.backgroundw, type="response")
plotmo(mars.model.pseudoabsence, letype="response")

#summary
summary(mars.model.backgroundw)
summary(mars.model.pseudoabsence)

#variable importance
mars.model.backgroundw.evimp<-evimp(mars.model.backgroundw)
mars.model.pseudoabsence.evimp<-evimp(mars.model.pseudoabsence)


#maps
mars.model.backgroundw.map<-predict(predictors, mars.model.backgroundw, type="response")
mars.model.pseudoabsence.map<-predict(predictors, mars.model.pseudoabsence, type="response")

#maps
plot(mars.model.backgroundw.map, main="background ponderado")
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
plot(maxent.model, type="logistic")

#map
maxent.model.map<-predict(predictors, maxent.model, type="logistic")
plot(maxent.model.map)


#################################################################
#CONDITIONAL INFERENCE TREES
#################################################################
formula.tree<-formula.mars

#fitting tree
ctree.model.pseudoabsence<-ctree(formula.tree, data=presence.pseudoabsence)

#plot tree
plot(ctree.model.pseudoabsence)

#predict
ctree.model.pseudoabsence.map<-predict(predictors, ctree.model.pseudoabsence, type="response")

#maps
plot(ctree.model.pseudoabsence.map)


#################################################################
#RANDOM FOREST
#################################################################
formula.tree<-formula.mars

#pseudoabsence
rf.model.pseudoabsence<-randomForest(formula.tree, data=presence.pseudoabsence, importance=TRUE, ntree=500, mtry=3, nodesize=10)

#ERROR PLOT
plot(rf.model.pseudoabsence)

#model summary
print(rf.model.pseudoabsence)

#response curves
plotmo(rf.model.pseudoabsence)

#variable importance
varImpPlot(rf.model.background)
varImpPlot(rf.model.absence)
varImpPlot(rf.model.pseudoabsence)

#as a table
importance(rf.model.pseudoabsence)

#to map
rf.model.pseudoabsence.map<-predict(predictors, rf.model.pseudoabsence, type="response")

#map
plot(rf.model.pseudoabsence.map, main="pseudoabsence")


#################################################################
#BOOSTED REGRESSION TREES
#################################################################

#weighted background
gbm.model.backgroundw=gbm.step(data=presence.background, gbm.x = names(predictors), gbm.y = "presence", family = "bernoulli", tree.complexity = 5, learning.rate = 0.01, bag.fraction = 0.5, site.weights=WeightPresenceBackground(presence.background[ , "presence"]))

#pseudoabsencce
gbm.model.pseudoabsence=gbm.step(data=presence.pseudoabsence, gbm.x = names(predictors), gbm.y = "presence", family = "bernoulli", tree.complexity = 5, learning.rate = 0.01, bag.fraction = 0.5)


#plot error
plot(gbm.model.backgroundw)
plot(gbm.model.pseudoabsence)

#variable importance
summary(gbm.model.backgroundw)
summary(gbm.model.pseudoabsence)

#plot
gbm.plot(gbm.model.backgroundw)
gbm.plot(gbm.model.pseudoabsence)

#also with plotmo
plotmo(gbm.model.backgroundw)
plotmo(gbm.model.pseudoabsence)

#prediction (check the ntrees argument)
gbm.model.backgroundw.map<-predict(predictors, gbm.model.backgroundw, n.trees=gbm.model.backgroundw$gbm.call$best.trees, type="response")

gbm.model.pseudoabsence.map<-predict(predictors, gbm.model.pseudoabsence, n.trees=gbm.model.pseudoabsence$gbm.call$best.trees, type="response")

#map

plot(gbm.model.backgroundw.map, main="backgroundw")
plot(gbm.model.pseudoabsence.map, main="pseudoabsence")

#CHECKING INTERACTIONS
gbm.model.pseudoabsence.interactions <- gbm.interactions(gbm.model.pseudoabsence)
#interaction matrix
gbm.model.pseudoabsence.interactions$interactions
#relevant interactions
gbm.model.pseudoabsence.interactions$rank.list

#plotting them specifically (change variable numbers as you need)
gbm.perspec(gbm.model.pseudoabsence, 1, 6)
gbm.perspec(gbm.model.pseudoabsence, 1, 5)



#################################################################
#ARTIFICIAL NEURAL NETWORKS
#################################################################

#number of hidden neurons
hidden.neurons=20

#weighted background (slower, reduce maxit)
nnet.model.backgroundw<-nnet(formula.tree, size=hidden.neurons, maxit=10000, data=presence.background, weights=WeightPresenceBackground(presence.background[ , "presence"]))

#pseudoabsence
nnet.model.pseudoabsence<-nnet(formula.tree, size=hidden.neurons, maxit=10000, data=presence.pseudoabsence)

#plotting neural networks
plotnet(nnet.model.backgroundw)
plotnet(nnet.model.pseudoabsence)

#variable importance
garson(nnet.model.backgroundw)
garson(nnet.model.pseudoabsence)

#response curves
plotmo(nnet.model.backgroundw)
plotmo(nnet.model.pseudoabsence)

#to maps
nnet.model.backgroundw.map<-predict(predictors, nnet.model.backgroundw)
nnet.model.pseudoabsence.map<-predict(predictors, nnet.model.pseudoabsence)

#maps
plot(nnet.model.backgroundw.map, main="weighted background")
plot(nnet.model.pseudoabsence.map, main="pseudoabsence")


################################################################
#SUPPORT VECTOR MACHINES
#################################################################
#weighted background
svm.model.backgroundw<-ksvm(formula.tree, data=presence.background, kernel="laplacedot", kpar=list(sigma=1), class.weigths=WeightPresenceBackground(presence.background[ , "presence"]))
svm.model.backgroundw #check the training error!

#pseudoabsence
svm.model.pseudoabsence<-ksvm(formula.tree, data=presence.pseudoabsence, kernel="laplacedot", kpar=list(sigma=1))
svm.model.pseudoabsence

#PREDICT
svm.model.backgroundw.map<-predict(predictors, svm.model.backgroundw, type="response")
svm.model.pseudoabsence.map<-predict(predictors, svm.model.pseudoabsence, type="response")

#map
plot(svm.model.backgroundw.map, main="weighted background")
plot(svm.model.pseudoabsence.map, main="pseudoabsence")



################################################################
################################################################
#ENSEMBLE MODEL FORECASTING
################################################################
################################################################

#BRICKS WITH ALL MODELS (by presence class + maxent)
models.brick<-brick("YOUR MDOELS HERE")
names(models.brick)<-c("YOUR MODEL NAMES HERE")

#MEAN
ensemble.mean<-calc(models.brick, mean)

#STANDARD DEVIATION
ensemble.sd<-calc(models.brick, sd)

#PLOT
plot(ensemble.mean, main="Mean")
plot(ensemble.sd, main="Standard deviation")


#ANALIZING THE ENSEMBLE
########################################
#brick to dataframe
models.df<-na.omit(as.data.frame(models.brick))
#correlation matrix
models.cor<-cor(models.df)
#cor to dist
modelsmodels.dis<-abs(as.dist(models.models.cor))
#plot
plot(hclust(1-models.dis))



##################################
##################################
#MODEL EVALUATION
##################################
##################################

#################################################################################
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
par(mfrow=c(1,2))
plot(maps.kfold.mean, main="Mean")
plot(maps.kfold.sd, main="Standard deviation")
 

####################################################################################
#USING BOOTSTRAP TO EVALUATE A MODEL
####################################################################################

#computing number of rows to get as evaluation data
presence.data<-presence.pseudoabsence
presence.data.nrow=nrow(presence.data)
evaluation.data.percentage<-40
evaluation.data.rows<-round((evaluation.data.percentage*presence.data.nrow)/100, 0)

#number of iterations
iterations=30

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



######################################################################################
######################################################################################
#PROJECTION (MESS ANALYSIS)
######################################################################################
######################################################################################

#MESS ANALYSIS
#WARNING: replace PROJECTION BRICK with whatever brick you want to use to project your model
mess.result<-mess(x="PROJECTION BRICK", v=presence[ , names(PROJECTION BRICK)], full=TRUE)
#naming and plotting result
names(mess.result)=c(names(PROJECTION BRICK), "mess")
plot(mess.result)


#MAX DEVIATED PREDICTOR PER CELL
mess.result.max<-which.max(mess.result)
plot(mess.result.max)

#making it a bit more easy to interpret
colors.map<-brewer.pal(n=8, name="Set1")
#plot
plot(mess.result.max, col=colors.map, legend=FALSE)
legend("bottomleft", c(names(PROJECTION BRICK)), pch=15, col=paleta.color, cex=1.6)

