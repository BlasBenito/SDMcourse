#LOADING REQUIRED LIBRARY
library(raster)

#LOADING INPUT DATA
####################
load("inputData.RData")
source("migration_function.R")

#presence data
str(presence)

#sdms
names(sdm)
par(mfrow=c(1,3))
plot(sdm[[1]], main="2000")
plot(sdm[[50]], main="2050")
plot(sdm[[100]], main="2100")

#cropping raster
sdm=crop(sdm, c(463000, 475000, 4102000, 4113000))
plot(sdm[[1]])
points(presence$x, presence$y)

#running the simulation
run.simulation(population=presence, sdm=sdm, max.dispersal.distance=2000, max.age=150, reproductive.age=20, max.fecundity=5, constant.climate=FALSE)
