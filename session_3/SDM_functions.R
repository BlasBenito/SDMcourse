


#############################################################################
#ReduceSpatialClustering
#This function reduces the spatial clustering of a set of presence records. It is intended to reduce spatial autocorrelation, and reduce sampling bias, specially at larger geographical scales.

#It requires two different arguments:
#data.table: a table with two fields representing latitude (named 'y') and longitude (named 'x')

#a minimum.distance value, provided in the same units as the coordinates, that will define the search radius when looking for pairs of coordinates within search distance to get rid of. Hint: the minimum distance can be extracted from the resolution of a raster containint the environmental factors, like "minimum.distance<-xres(v.brick.20km)"
ReduceSpatialClustering = function(data, minimum.distance){

#count rows
row<-1


#repite la operación hasta que se cumple la condición de salida
repeat{
  
  #contenido de la fila (para no tirar de toda la tabla en todas las operaciones)
  f<-data[row, ]
  
  #genera los límites de la cuadrícula de búsqueda
  ymax<-f$latitude + minimum.distance
  ymin<-f$latitude - minimum.distance
  xmax<-f$longitude + minimum.distance
  xmin<-f$longitude - minimum.distance
  
  #selecciona de la tabla los datos con coordenadas dentro del rectángulo que no tienen las mismas coordenadas que la fila con la que estamos trabajando, y las elimina de la tabla
  data<-data[!((data$latitude <= ymax) & (data$latitude >= ymin) & (data$longitude <= xmax) & (data$longitude >= xmin) & (data$latitude != f$latitude | data$longitude != f$longitude)), ]
  
  #estima de filas por procesar
  print(paste("Processed rows: ", row, " out of ", nrow(data), sep=""))
  
  #suma 1 al contador de la fila
  row<-row+1
  
  #condición de salida cuando llega a la última fila
  if(row>=nrow(data))break
}

return(data)

}






#############################################################################
#WEIGHT PRESENCE/BACKGROUND DATA
WeightPresenceBackground=function(presence.column){
  
  #computing weight for presences
  n.presences=sum(presence.column)
  print(paste("Presence points = ", n.presences, sep=""))
  weight.presences=1/n.presences
  print(paste("Weight for presences = ", weight.presences, sep=""))
  
  n.background=length(presence.column)-n.presences
  print(paste("Background points = ", n.background, sep=""))
  weight.background=1/n.background
  print(paste("Weight for background = ", weight.background, sep=""))
  
  #generamos un vector con los los pesos
  weights<-c(rep(weight.presences, n.presences), rep(weight.background, n.background))
  #return(weights)
}


#############################################################################
#http://modtools.wordpress.com/2013/08/14/dsquared/ 
# Linear models come with an R-squared value that measures the proportion of variation that the model accounts for. The R-squared is provided with summary(model) in R. For generalized linear models (GLMs), the equivalent is the amount of deviance accounted for (D-squared; Guisan & Zimmermann 2000), but this value is not normally provided with the model summary. The Dsquared function, now included in the modEvA package (Barbosa et al. 2014), calculates it. There is also an option to calculate the adjusted D-squared, which takes into account the number of observations and the number of predictors, thus allowing direct comparison among different models (Weisberg 1980, Guisan & Zimmermann 2000).
Dsquared <- function(model, adjust = TRUE) {
  # version 1.1 (13 Aug 2013)
  # calculates the explained deviance of a GLM
  # model: a model object of class "glm"
  # adjust: logical, whether or not to use the adjusted deviance taking into acount the nr of observations and parameters (Weisberg 1980; Guisan & Zimmermann 2000)
  d2 <- (model$null.deviance - model$deviance) / model$null.deviance
  if (adjust) {
    n <- length(model$fitted.values)
    p <- length(model$coefficients)
    d2 <- 1 - ((n - 1) / (n - p)) * (1 - d2)
  }
  return(d2)
}  # end Dsquared function



