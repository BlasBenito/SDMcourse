############
#INTRO TO R
############

#DEFINE WORKING DIRECTORY
setwd("path to the folder here")

#INSTALL SOME BASIC PACKAGES
####################
install.packages(c("car","plotmo","raster","corrgram", "dismo"), dep=TRUE)

#LOADING PACKAGES
library(car)
library(plotmo)
library(raster)
library(corrgram)
library(dismo)


######
#HELP!
######
#web page
help.start() 

#search for related functions
help.search("regression") 

#look up in the R site
RSiteSearch("regression") 

#manual of a particular function
help(lm) 

#manual of a library
help(package = "stats") 

#exampel of a particular function (if it exists)
example(lm) 

#demos
demo(persp)
demo(graphics)



###############
#DATA TYPES
###############
#NUMERIC
x<-10.5
x
is.numeric(x)
is.character(x)

#removing objects
rm(x)


#CHARACTER
y<-"dog"
y
is.numeric(y)
is.character(y)

rm(y)


#MISSING VALUES
x<-c(1, 2, 3, 4, 5, 6, 7, NA) #NA -> Not Available
is.na(x)
mean(x)
help(mean)
mean(x, na.rm=T)

rm(x)


#rounding
round(5.9, digits=0)

#
abs(-5.9)


#####################
#DATA STRUCTURES
#####################

#VECTOR
#######
z<-c(1, 2, 3, 4, 5, 6)
z
length(z)
is.vector(z)
z<-1:6
z
z<-seq(1, 6, 1)
z

#some functions
max(z)
min(z)
sum(z)
mean(z)
median(z)
range(z)
sort(z)

rm(z)


#DATA FRAME
###########
data(cars)
cars
is.data.frame(cars)
str(cars) 

#USING INDEXES
cars$speed
cars$dist
cars[, "speed"]
cars[, "dist"]
cars[10,]
cars[5:15,]
cars[5:15, "speed"]

#CONDITIONAL SUBSETTING
cars2<-cars[cars$dist < 60, ]
cars2

rm(cars)
rm(cars2)

#CREATING A DATAFRAME
x<-c(10, 20, 40, 80)
y<-c("s","s","n","n")
xy<-data.frame(x, y)
names(xy)<-c("value","group")

#creating an empty dataframe
empty.data.frame<-dataframe(especie=as.character(), latitude=numeric(), longitude=numeric(), stringsAsFactors=FALSE)

rm(x, y, xy, d1, empty.data.frame)


#import data
climate.data<-read.table("./data/climate.csv", header=TRUE, sep=",")
str(climate.data)



###############################
#RASTER DATA (TO MANAGE PREDICTORS!)
###############################
#IMPORTING RASTER MAPS
#lista de variables
list.predictors <- list.files(path="./data/maps", pattern='bio', full.names=TRUE)
list.predictors
#stack (predictors NOT in memory)
predictors.stack <- stack(list.predictors)

#plotting the stack
plot(predictors.stack)

#names
names(predictors.stack)

#plotting a single map
plot(predictors.stack, "bio5")
plot(predictors.stack[["bio5"]])

#brick (data IN memory, faster, but you need enough RAM memory)
predictors.brick <- brick(predictors.stack)
plot(predictors.brick)

#features
predictors.brick

#slots
slotNames(predictors.brick)

#checking one slot
slot(predictors.brick, "extent")

#other way to access an slot
predictors.brick@extent

#another way to check the extent
extension<-extent(predictors.brick)
extension

#coordinate system, check http://www.spatialreference.org/
reference.system<-CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")
projection(predictors.brick)<-reference.system
predictors.brick

#resolution
xres(predictors.brick)
yres(predictors.brick)
#in km? (fast approximation)
xres(predictors.brick)*111.19

#map algebra
temperature.range<-predictors.brick[["bio5"]] - predictors.brick[["bio6"]]
names(temperature.range)<-"temperature.range"

#adding it to the brick
predictors.brick<-addLayer(temperature.range, predictors.brick)
plot(predictors.brick)

#continuous to binary
plot(predictors.brick[["temperature.range"]] > 250)

#brick to dataframe
predictors.df<-as.data.frame(predictors.brick)
str(predictors.df)
#removing NA values
predictors.df<-na.omit(predictors.df)
str(predictors.df)
#NO BORRAMOS predictors.df, VAMOS A USARLO LUEGO

#getting one raster from the brick
bio5<-predictors.brick[["bio5"]]
plot(bio5)

#decreasing resolution (aggregation)
bio5.lowres<-aggregate(bio5, fact=5, fun=mean)
plot(bio5.lowres)

#map cropping
#xmin, xmax, ymin, ymax
bio5.lowres.iberian<-crop(bio5.lowres, c(-10, 5, 35, 45))
plot(bio5.lowres.iberian)

#change projection to ETRS89 / UTM zone 30N at 50KM resolution
bio5.lowres.iberian.utm<-projectRaster(bio5.lowres.iberian, crs="+proj=utm +zone=30 +ellps=GRS80 +units=m +no_defs ", res=c(50000, 50000))

#vemos la resolución
yres(bio5.lowres.iberian.utm)
xres(bio5.lowres.iberian.utm)

#ploteamos los dos
par(mfrow=c(1,2))
plot(bio5.lowres.iberian)
plot(bio5.lowres.iberian.utm)

#saving maps
writeRaster(bio5.lowres.iberian.utm, filename="bio5_lowres_iberian.utm.asc", format="ascii", overwrite=TRUE)



#######################
#CONTROL STRUCTURES
#######################

#FOR LOOP
###########

#numeric sequence
for (year in 1980:2000){
  print(year) #you can do whatever you want here using the "year" variable.
}

rm(year)

#over character vector
genuses<-c("Abies","Fagus","Pinus","Quercus")
for (genus in genuses){
  print(genus)
}

rm(genuses)

#plotting maps
pdf("predictors.pdf", width=20, height=15, pointsize=20)
for (predictor in names(predictors.brick)){
  plot(predictors.brick, predictor)
}
dev.off()


#stopping a loop
for (genus in genuses){
  print(genus)
  if(genus=="Pinus"){break} #you can use whatever condition you need here
}


#DATA ANALYSIS
########################

#sampling rows from a big table
n.rows=nrow(predictors.df)
n.rows
predictors.df.small<-predictors.df[sample(n.rows, 1000), ]
nrow(predictors.df.small)

#saving the new dataframe (we don't need to do it, just an example)
write.table(predictors.df.small, "table.csv", sep=",", row.names=FALSE, col.names=TRUE, quote=FALSE)

#CORRELOGRAM
corrgram(predictors.df.small, order=TRUE, lower.panel=panel.ellipse, upper.panel=panel.pts, text.panel=panel.txt)

#CORRELATION MATRIX
variables.cor<-cor(predictors.df.small)
variables.cor

#CORRELATION MATRIX AS DENDROGRAM
variables.dist=as.dist(1-abs(variables.cor))
plot(hclust(variables.dist))


#LINEAR MODEL
#question: can we predict bio5 (summer temperature) using bio6 (winter temperature) and bio12 (rainfall) as predictors?
#model fit
model<-lm(bio5 ~ bio6 + bio12, data=predictors.df.small)
summary(model)
#response curves
plotmo(model, level=0.95, all2=TRUE)
#lets make a prediction to a raster map
help(predict) #there are a few predict functions
#we wan to use the version of the raster package
prediction<-raster::predict(object=predictors.brick, model=model)
plot(prediction)
#plotting errors
plot(prediction - predictors.brick[["bio5"]])

##################
#CREATING FUNCTIONS
##################
#let's create some data
x<-rnorm(n=1000, mean=20, sd=5)

#histogram
hist(x)

#density
plot(density(x))

#sum values
x.sum<-sum(x)

#number of cases
n<-length(x)

#average
mean<-x.sum/n

#a function to compute the average of a vector
average<-function(input.data){
  result=sum(input.data)/length(input.data)
  return(result)
}

#using it
average(x)

#it already exists, but you got the point
mean(x)

##############
#HOUSE KEEPING
##############
#objects created during the session
objects()
#listed in the Environment tab in Rstudio

#saving work environment (Rstudio does it automatically)
save(list=ls(all=TRUE), file="workspace.Rdata")

#remove everything from the workspace
rm(list=ls())

#load workspace
load("workspace.Rdata")

#look up libraries and dataframes in the working space
search()

