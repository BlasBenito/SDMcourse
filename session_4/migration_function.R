# population=presence
# sdm=sdm
# max.dispersal.distance=500
# max.age=150
# reproductive.age=20
# max.fecundity=10

run.simulation=function(population, sdm, max.dispersal.distance, max.age, reproductive.age, max.fecundity, constant.climate){
  
  require(raster)
  
  x11()

  #SOME GRAPHICAL PARAMETERS
  ##########################
  #colors
  sdm.colors=gray.colors(99)

  #PARAMETERS
  ###########
  #TIME
  years=2000:2100
  years.length=length(years)
  
  #AGE (Random start)
  presence$age=sample(1:max.age, size=nrow(presence), replace=TRUE)
  
  #POP SIZE PER YEAR
  pop.size=vector()
  mature.size=vector()
  years.vector=vector()
  
  #GENERATING DISTRIBUTION FOR DISPERSAL SIMULATION
  dispersal.distances <- ceiling(rnorm(10000, mean=0, sd=(max.dispersal.distance/3)))

  #SIMULATION
  #############################################################################
  #LOOPING THROUGH YEARS
  #start with burn-in
  for (year.index in c(rep(1,50), 1:years.length)){
  # for (year.index in 1:years.length){
    
    #getting the current year
    year=years[year.index]
    cat(paste("Simulating year", year, sep=" "), sep="\n")
    
    #filling the year vector
    years.vector[year.index]=year
    
    #INCREASE AGE AND SIMULATING MORTALITY BY OLD AGE
    ##################################################################################
    presence$age=presence$age+1
    presence=presence[which(presence$age <= max.age), ]

    #HABITAT SUITABILITY CHANGE
    ##################################################################################
    if (constant.climate==FALSE){
      presence$suitability=as.vector(extract(sdm[[year.index]], presence[, c("x", "y")]))
    }
    
    if (constant.climate==TRUE){
      presence$suitability=as.vector(extract(sdm[[1]], presence[, c("x", "y")]))
    }
    
    #removing presence with no environmental data (seed could go out of bounds)
    presence=na.omit(presence)
    
    #mortality by low suitability
    presence=presence[which((presence$suitability + (presence$age / 2)) > sample(0:100, size=nrow(presence), replace=TRUE)), ]
    
    #DENSITY DEPENDENT PROCESS
    ###############################################################################
    presence = ReduceSpatialClustering(presence)
    
    #REPRODUCTION AND DISPERSAL
    ##################################################################################
    #selecting reproductives
    reproductives=presence[which(presence$age >= reproductive.age), ]
    
    #producing seeds
    if (nrow(reproductives) > 0){
    
      #computing fecundity based on habitat suitability and max.fecundity
      reproductives$fecundity=ceiling((reproductives$suitability*max.fecundity)/100)
      
      #producing seeds
      seeds=reproductives[rep(row.names(reproductives), reproductives$fecundity), ]
      seeds$age=0
      seeds$fecundity=NULL
      
      #dispersal
      seeds$x=seeds$x + sample(dispersal.distances, size=nrow(seeds), replace=TRUE)
      seeds$y=seeds$y + sample(dispersal.distances, size=nrow(seeds), replace=TRUE)
      
      #PLOTTING
      #-----------------------------------------------------------------------------
      par(mfrow=c(1,2))
      
      if (constant.climate==TRUE){
      plot(sdm[[1]], col=sdm.colors, breaks=seq(min(minValue(sdm[[1]])), max(maxValue(sdm[[1]])), length.out=100), legend=FALSE, main=year)
      }
      
      if (constant.climate==FALSE){
        plot(sdm[[year.index]], col=sdm.colors, breaks=seq(min(minValue(sdm[[year.index]])), max(maxValue(sdm[[year.index]])), length.out=100), legend=FALSE, main=year)
      }
      
      #presence
      points(presence[which(presence$age < reproductive.age), "x"], presence[which(presence$age < reproductive.age), "y"], cex=0.8, pch=20, col="green")
      points(presence[which(presence$age >= reproductive.age), "x"], presence[which(presence$age >= reproductive.age), "y"], cex=1.5, pch=20, col="red")
      #seeds
      points(seeds$x, seeds$y, cex=0.3, pch=20, col="yellow")
      #pop.size
      pop.size[year.index]<-nrow(presence)
      plot(years.vector, pop.size, xlab="Year", ylab="Individuals", type="l", main="Population", ylim=c(0, max(pop.size)))
      #mature individuals
      mature.size[year.index]<-nrow(presence[presence$age > reproductive.age, ])
      lines(years.vector, mature.size, xlab="", ylab="", type="l", col="red")
      #-----------------------------------------------------------------------------
      
      #joining seeds to the population
      presence=rbind(presence, seeds)
      
    }

    #CHECKING IF THE POPULATION CRASHED
    if (nrow(presence)==0) {
      message("Population crashed, end of simulation.")
      break #end of simulation
    }
    
    #PREVENTING THE POPULATION TO BECOME TOO BIG
    if (nrow(presence)>5000) {
      message("Population beyond safe performance bounds, simulation aborted.")
      break #end of simulation
    }
    
    #we give it time to plot (but this slows down the simulation, only for demonstraton purposes)
    # Sys.sleep(0.1)
    
  }#end of loop

}#end of function



ReduceSpatialClustering = function(data){
  
  # data<-data[order(data$age, decreasing=TRUE), ]
  
  #count rows
  row<-1
  
  
  #repite la operación hasta que se cumple la condición de salida
  repeat{
    
    #contenido de la fila (para no tirar de toda la tabla en todas las operaciones)
    f<-data[row, ]
    
    #minimum distance according age
    minimum.distance = f$age/2
    
    #genera los límites de la cuadrícula de búsqueda
    ymax<-f$y + minimum.distance
    ymin<-f$y - minimum.distance
    xmax<-f$x + minimum.distance
    xmin<-f$x - minimum.distance
    
    #selects the neighbors
    data<-data[!((data$y <= ymax) & (data$y >= ymin) & (data$x <= xmax) & (data$x >= xmin) & (data$y != f$y | data$x != f$x)), ]
    
    #suma 1 al contador de la fila
    row<-row+1
    
    #condición de salida cuando llega a la última fila
    if(row>=nrow(data))break
  }
  
  return(data)
  
}
