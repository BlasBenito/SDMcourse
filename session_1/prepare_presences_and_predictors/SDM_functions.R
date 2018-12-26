

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


