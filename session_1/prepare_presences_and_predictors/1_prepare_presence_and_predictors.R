###################################
#PREPARE PRESENCES AND PREDICTORS
###################################

setwd("path")

#CARGAMOS FUNCIONES (iremos dando un vistazo a las funciones cuando las usemos)
source("SDM_functions.R")

#install some new packages
install.packages(c("rgeos", "HH"), dep=TRUE)

#loading libraries
library(raster)
library(HH)
library(rgeos)
library(dismo)



##########################################################################
#SAME EXTENSION AND RESOLUTION FOR RASTER MAPS
##########################################################################
#this is just an example, the predictors we'll use are already prepared.

#importamos mapas
elev<-raster("./mask_and_region/elevation.asc")
hfp<-raster("./mask_and_region/hfp.asc")
ndvi<-raster("./mask_and_region/ndvi.asc")

#plot
par(mfrow=c(1,3))
plot(elev)
plot(hfp)
plot(ndvi)

#extension
extent(elev)
extent(hfp)
extent(ndvi)

#resolution
xres(elev)
xres(hfp)
xres(ndvi)


#1 - EQUAL RESOLUTION
help(resample)
elev2<-resample(x=elev, y=ndvi, method="bilinear")
#comparing extension
extent(ndvi)
extent(elev2)
#resolution
xres(ndvi)
xres(elev2)
#plot
par(mfrow=c(1,2))
plot(elev)
plot(elev2)

#hfp
hfp2<-resample(x=hfp, y=ndvi, method="bilinear")

par(mfrow=c(1, 3))
plot(ndvi)
plot(elev2)
plot(hfp2)


#2 - PREPARING A MASK
#hint: we use multiplication to propagate null cells
null.cells=ndvi*elev2*hfp2 #la multiplicación de todos los mapas propaga los valores nulos
plot(null.cells) #this map has all the null cells of every map


#3 - APPLYING A MASK
#all predictors to a brick
predictors.brick<-brick(ndvi, elev2, hfp2)
names(predictors.brick)<-c("ndvi", "elev", "hfp")
#applying mask
predictors.brick<-mask(predictors.brick, null.cells)
plot(predictors.brick) #so far so good


#4 - CROPPING
#to the map with the smaller extension
extent(elev)
extent(hfp) #this one
predictors.brick<-crop(x=predictors.brick, y=extent(hfp))
plot(predictors.brick)

#5 - EXPORTING MAPS (don't run this)
writeRaster(predictors.brick[["elev"]], filename="elev_final.asc", format="ascii", overwrite=TRUE)
writeRaster(predictors.brick[["hfp"]], filename="hfp_final.asc", format="ascii", overwrite=TRUE)
writeRaster(predictors.brick[["ndvi"]], filename="ndvi_final.asc", format="ascii", overwrite=TRUE)

rm(elev, elev2, hfp, hfp2, ndvi, null.cells, predictors.brick)


################################################################
################################################################
#REDUCING COLLINEARITY OF PREDICTORS
################################################################
################################################################

#list of files to import
list.predictors <- list.files(path="./predictors",pattern='*.asc', full.names=TRUE)

#creating an stack
predictors <- stack(list.predictors)
names(predictors)

#resolution
res.degrees<-xres(predictors)
res.km<-res.degrees*111.19
res.km

#plot
plot(predictors, maxnl=length(names(predictors)))

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
x11()
plot(predictors.cluster)

#select one predictor per cluster below the 0.5 line
selected.predictors<-c("landcover_veg_bare", "bio5", "bio15", "bio14", "ndvi_average", "landcover_veg_herb", "landcover_veg_tree", "topo_diversity", "bio12", "bio6", "human_footprint", "bio7", "ndvi_range")

#df with those predictors
predictors.df2<-predictors.df[ , selected.predictors]

#VARIANCE INFLATION FACTOR (library HH)
vif.result<-vif(predictors.df2)
vif.result

#remember, the maximum recommended value is 5
#I decide to remove bio7
predictors.df2$bio7<-NULL
vif.result<-vif(predictors.df2)
vif.result

#now landcover_veg_bare
predictors.df2$landcover_veg_bare<-NULL
vif.result<-vif(predictors.df2)
vif.result

#bio15
predictors.df2$bio15<-NULL
vif.result<-vif(predictors.df2)
vif.result

#landcover_veg_herb (because I want to keep landcover_veg_tree)
predictors.df2$landcover_veg_herb<-NULL
vif.result<-vif(predictors.df2)
vif.result

#bio12 (I want to keep bio6)
predictors.df2$bio12<-NULL
vif.result<-vif(predictors.df2)
vif.result

#new list of selected predictors
selected.predictors<-names(vif.result)

#selected predictors to brick
predictors<-brick(predictors[[selected.predictors]])


#removing stuff we don't need anymore
rm(predictors.correlation, predictors.df, predictors.df2, list.predictors, res.degrees, vif.result, predictors.cluster, predictors.dist, selected.predictors, res.km)

#cleaning the ram memory
gc()


################################################################
################################################################
#PREPARING PRESENCE DATA
################################################################
################################################################

#importing data
presence.complete<-read.table("./presence/occurrence.txt",header=T, sep='\t', fill=TRUE, check.names=TRUE, stringsAsFactors=FALSE)

#copy (if we make a mistake, we start here again)
presence<-presence.complete
str(presence)

#CLEANING THE DATA
#####################

#REMOVING PRESENCES OUTSIDE THE STUDY REGION
#--------------------------------------
#extracting values of the presences over the predictors
presence.over.predictors<-data.frame(extract(x=predictors, y=presence[ , c("longitude","latitude")]))

#joining it with he presence data
presence<-data.frame(presence, presence.over.predictors)
rm(presence.over.predictors)

#removing records with NA values for any of the predictors (all predictors have the same null cells)
presence<-presence[!(is.na(presence$bio5)), ]

#plot
plot(predictors[[1]], col="gray80")
points(presence$longitude,presence$latitude)


#REMOVING OTHER SPECIES (CAN BE COMMON IN GBIF DATASETS)
#------------------------------------
unique(presence$specific_epithet)
unique(presence$species_id)
presence<-presence[presence$species_id==2433433, ]
presence$specific_epithet<-"Ursus arctos"


#REMOVING FOSSIL RECORDS IF ANY
#-------------------------------------------
barplot(table(presence$basis_of_record))
presence<-presence[presence$basis_of_record != "FOSSIL_SPECIMEN", ]
barplot(table(presence$basis_of_record))


#REMOVING OLD RECORDS
#---------------------------
hist(presence$year)
#how many records have no year?
sum(is.na(presence$year))
#where are these records?
plot(predictors[[1]], col="gray80")
points(presence$longitude,presence$latitude)
points(presence[is.na(presence$year), ]$longitude, presence[is.na(presence$year), ]$latitude, col="red", cex=1.2)
#data without year is quite redundant, we get rid of it
presence<-presence[!is.na(presence$year), ]
#check
sum(is.na(presence$year))
#we keep data with year > 1970 (congruent with most predictors)
presence<-presence[presence$year >= 1970, ]


#CHECK COORDINATES PRECISION
#--------------------------
#checking the precision field
str(presence$coordinate_precision)
#character type for numeric variable!
presence$coordinate_precision<-as.numeric(presence$coordinate_precision)
#what values do we have?
unique(presence$coordinate_precision)
barplot(table(presence$coordinate_precision))
#how many records with no data?
sum(is.na(presence$coordinate_precision))
#where?
plot(predictors[[1]], col="gray80")
points(presence$longitude,presence$latitude)
points(presence[is.na(presence$coordinate_precision), ]$longitude, presence[is.na(presence$coordinate_precision), ]$latitude, col="red", cex=1.2)
#we are in trouble, all Central European data has no coordinates precision...
#I decided to keep it.


#INCREASING DISTANCE AMONG CONSECUTIVE PRESENCES
#-----------------------------------------------------
#to reduce autocorrelation

#predictor resolution
xres(predictors)
yres(predictors)
res.degrees<-xres(predictors)
#empty cells between consecutive presence records
empty.cells<-1
#minimum distance between consecutive presence records
minimum.distance<-res.degrees*empty.cells
#in km
minimum.distance*111.19

#let's see how the function works by opening a plot
x11()
plot(predictors[[1]], col="gray80", ext=c(-8, -3, 42, 44))
points(presence$longitude,presence$latitude)

#functon to reduce autocorrelation
presence=ReduceSpatialClustering(data=presence, minimum.distance=minimum.distance)

#ending the plot
points(points(presence$longitude,presence$latitude, col="red", cex=1.2, lwd=2))


#REMOVING RECORDS WITH DUPLICATED COORDINATES
#-----------------------------------------
#finding duplicates
duplicated.records<-duplicated(presence[ , c("latitude", "longitude")])
#how many?
length(duplicated.records[duplicated.records==TRUE])
#getting non-duplicates only
presence<-presence[!duplicated.records, ]


#FINAL PRESENCES
plot(predictors[[1]], main="Presencias final", col="gray80")
points(presence$longitude,presence$latitude, cex=0.3)

#how many presence records?
nrow(presence)

#removing and renaming some columns
#cambiamos latitude y longitude por y y x (más corto)
colnames(presence)[5]<-"y"
colnames(presence)[6]<-"x"
presence$specific_epithet<-NULL
presence$basis_of_record<-NULL
presence$species_id<-NULL
presence$country_code<-NULL
presence$year<-NULL
presence$coordinate_precision<-NULL
#adding a presence column
presence$presence<-1


#########################
#########################
#PREPARING "ABSENCE" DATA
#########################
#########################

#BACKGROUND DATA
###################################
#getting the coordinates of non-null cells
background <- as.data.frame(predictors[["bio5"]], xy=TRUE)
background <- na.omit(background)
background$bio5 <- NULL
#extracting values of predictors
background.predictors<-data.frame(extract(predictors, background))
#merging with the coordinates fields
background<-cbind(background, background.predictors)
#adding the presence column with 0 values
background$presence<-0
#adding it to the presence data
presence.background<-rbind(presence, background)

#plot
x11()
plot(predictors[[1]], main="Background", col="gray80")
points(presence.background[presence.background$presence==1, ]$x, presence.background[presence.background$presence==1, ]$y, cex=0.3, col="red")
points(presence.background[presence.background$presence==0, ]$x, presence.background[presence.background$presence==0, ]$y, pch=20, cex=0.01, col="gray40")


#PSEUDOABSENCES
##############################
#Pseudoabsences shouldn't overlap with presences, therefore we will use "presence.complete" as a template, since we have removed presences from there.
#we want as many pseudoabsences as presences we have in the modelling dataset
pseudoabsence <- data.frame(randomPoints(mask=predictors, n=nrow(presence), p=presence.complete[ , c("latitude","longitude")], excludep=TRUE))
str(pseudoabsence)

#rest of the stuff we did with the background data
pseudoabsence.predictors<-data.frame(extract(predictors, pseudoabsence))
pseudoabsence<-cbind(pseudoabsence, pseudoabsence.predictors)
pseudoabsence$presence<-0
presence.pseudoabsence<-rbind(presence, pseudoabsence)

#plot
plot(predictors[[1]], main="Pseudoausencia", col="gray80")
points(presence.pseudoabsence[presence.pseudoabsence$presence==1, ]$x, presence.pseudoabsence[presence.pseudoabsence$presence==1, ]$y, cex=0.3, col="red")
points(presence.pseudoabsence[presence.pseudoabsence$presence==0, ]$x, presence.pseudoabsence[presence.pseudoabsence$presence==0, ]$y, pch=20, cex=0.3, col="gray40")


#ABSENCES
#########
#We don't have true absences, but we will generate points resembling true absences just for pedagogical purposes (aka "don't do this at home").

#DEFINING A RADIUS AROUND KNOWN PRESENCES
radius=100000 #en metros, 100km

#we convert the presence data into an "sp" (spatial) object
presence.sp<-presence[ , c("x", "y")]
coordinates(presence.sp)<-c("x", "y")

#generating a polygon around presences using that radius
buffer<-circles(presence.sp, d=radius, lonlat=TRUE)
buffer.dissolve<-gUnaryUnion(buffer@polygons)
plot(buffer.dissolve)

#extracting the IDs of the cells inside of the buffer
buffer.cells<-unlist(cellFromPolygon(predictors, p=buffer.dissolve))

#we create a template for a mask from a raster file
mask.temp<-raster(predictors)

#we create a vector with as many empty values as cells has the mask template
mask.values<-rep(NaN, ncell(mask.temp))

#giving value 1 to all cells inside the buffer
mask.values[buffer.cells]<-1

#transferring these values to the mask template
mask.temp<-setValues(mask.temp, values=mask.values)
plot(mask.temp)

#I don't even remember what the heck I did here
#actually I do, I intersected the mask of the predictors (predictors[[1]] > 0)) with the mask created from the buffer
mask<-mask.temp*(predictors[[1]] > 0)
plot(mask)

#now we generate as many absences as presences inside of the buffer, without overlapping them
absence <- data.frame(randomPoints(mask=mask, n=nrow(presence), p=presence.complete[ , c("latitude","longitude")], excludep=TRUE))
absence.predictors<-data.frame(extract(predictors, absence))
absence<-cbind(absence, absence.predictors)
absence$presence<-0
presence.absence<-rbind(presence, absence)

#plot
plot(predictors[[1]], main="Ausencia", col="gray80")
points(presence.absence[presence.absence$presence==1, ]$x, presence.absence[presence.absence$presence==1, ]$y, cex=0.3, col="red")
points(presence.absence[presence.absence$presence==0, ]$x, presence.absence[presence.absence$presence==0, ]$y, pch=20, cex=0.3, col="gray40")


#saving all the important data in an Rdata file for the next session
save(presence, presence.absence, presence.background, presence.pseudoabsence, predictors, file="predictors_and_presence.Rdata")
