
## ===============================================================================
## R-Script - Usefull functions to help in repetitive tasks of Remote Sensing data
## author: Javier Lopatin
## mail:javierlopatin@gmail.com & javier.lopatin@kit.edu
## last changes: 12/02/2016
## ===============================================================================


## -----------------------------------------------------##
## List the names of the rasters in a folder            ##
## -----------------------------------------------------##

rasterListNames <- function(fileExtantion, folder){
  # make a list of all fileExtantion files
  rast_list = list.files(folder, pattern = fileExtantion)
  # delete the ".dat" from the name
  rast_list = gsub('.{4}$', '', rast_list)
  return(rast_list)
}

## -----------------------------------------------------##
## List and load the rasters contained in a folder      ##
## -----------------------------------------------------##

rasterList <- function(fileExtantion, folder, rasterNames=NULL){
  # make a list of all fileExtantion files
  rast_list = list.files(folder, pattern = fileExtantion)
  # import rasters
  setwd(file.path(home, folder))
  rasterlist <- list()
  for(i in 1:length(rast_list)){
    rast <- stack(rast_list[i])
    if (is.null(rasterNames)){
      names(rast) <- rasterNames
      }
    rasterlist[[i]] <- rast
  }
  setwd(home)
  return(rasterlist)
}


## -------------------------------------------------##
## Plot all spectra signatures                      ##
## x = the spectral bands ordered in columns        ##
## y = the band wavelength, e.g. seq(475, 875, 10)  ##
## -------------------------------------------------##

plot.spectra <- function(x, y, ...){
  # set parameters
  N2 = length(x[,1]) # observations
  refl_2<-t(x)       # transpose data

  plot(y, refl_2[,1],type="l", axes=T, ylim=c(min(x), max(x)), ...)
  for(i in 1:N2-1){  lines(y, refl_2[, i+1])  }
}


## -----------------------------------------------------##
## Plot all first derivative signatures                 ##
## x = derivative results from the first.deriv function ##
## -----------------------------------------------------##

plot.deriv <- function(x, ...){
  plot(x[,1], x[,2],type="l", axes=T, ylim=c(min(apply(x[, 2:length(x[1, ])],1,min)), max(apply(x[, 2:length(x[1, ])],1,max))), ...)
  for(i in 1:(length(x[1,])-1)){  lines(x[,1], x[, i+1])}
}
