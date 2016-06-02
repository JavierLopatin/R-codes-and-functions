## ------------------------------------------------------------------------------------------------------##
## Description: This function allows to convert a raster layer into a polygon using directly the         ##
## GDAL python script gdal_polygonize.py. This will outperform in speed the rasterToPolygon function      ##   
## of the raster package.                                                                                ##  
## Requirements: GDAL install in your computer                                                           ##
## ------------------------------------------------------------------------------------------------------##

####### Version for Unix-like OS #######

# Source: http://johnbaumgartner.wordpress.com/2012/07/26/getting-rasters-into-shape-from-r/
# Based on: Lyndon Estes code, http://r-sig-geo.2731867.n2.nabble.com/Memory-management-with-rasterToPolygons-raster-and-rgeos-td7153049.html
polygonizeR <- function(x, outshape=NULL, attname='DN', gdalformat = 'ESRI Shapefile', quiet=TRUE) {
  py.c <- Sys.which('gdal_polygonize.py')
  if (!length(py.c)) stop("Can't find gdal_polygonize.py on your system.")
  if (!is.null(outshape)) {
    outshape <- sub('\\.shp$', '', outshape)
    f.exists <- file.exists(paste(outshape, c('shp', 'shx', 'dbg'), sep='.'))
    if (any(f.exists)) stop(sprintf('File already exists: %s', 
      toString(paste(outshape, c('shp', 'shx', 'dbg'), sep='.')[f.exists])), call.=FALSE)
  } else outshape <- tempfile()
  if (is(x, 'Raster')) {
    require(raster)
    writeRaster(x, {f <- tempfile(fileext='.asc')})
    rast.nm <- normalizePath(f)
  } else if (is.character(x)) {
    rast.nm <- normalizePath(x) 
  } else stop('x must be either a file path (as a character string), or a Raster object.')
  full.c <- sprintf("%1$s %2$s -f '%3$s' %4$s.shp %4$s %5$s", py.c, rast.nm, gdalformat, outshape, attname)
  system(full.c)
  shp <- readOGR(dirname(outshape), layer = basename(outshape), verbose=!quiet)
  return(shp)
}

# #### Example ####
# # Import raster
# download.file('http://dl.dropbox.com/u/1058849/blog/NEcountries.asc.zip',
#   destfile={f <- tempfile()}, quiet=TRUE, cacheOK=FALSE)
# unzip(f, exdir={d <- tempdir()})
# library(rasterVis)
# r <- raster(file.path(d, 'NEcountries.asc'), crs=CRS('+proj=longlat'))
# levelplot(r, margin=FALSE, col.regions=rainbow)
# 
# # polygonize and plot
# p <- rasterToPolygons(r, dissolve=TRUE)
# spplot(p, col.regions=rainbow(200))

####### Version for Windows #######

win.polygonizeR <- function(x, outshape=NULL, gdalformat = 'ESRI Shapefile', 
                        pypath=NULL, readpoly=TRUE, quietish=TRUE) {
  # x: an R Raster layer, or the file path to a raster file recognised by GDAL
  # outshape: the path to the output shapefile (if NULL, a temporary file will be created)
  # gdalformat: the desired OGR vector format
  # pypath: the path to gdal_polygonize.py (if NULL, an attempt will be made to determine the location
  # readpoly: should the polygon shapefile be read back into R, and returned by this function? (logical)
  # quietish: should (some) messages be suppressed? (logical)
  if (isTRUE(readpoly)) require(rgdal)
  if (is.null(pypath)) {
    pypath <- Sys.which('gdal_polygonize.py')
  }
  ## The line below has been commented:
  # if (!file.exists(pypath)) stop("Can't find gdal_polygonize.py on your system.") 
  owd <- getwd()
  on.exit(setwd(owd))
  setwd(dirname(pypath))
  if (!is.null(outshape)) {
    outshape <- sub('\\.shp$', '', outshape)
    f.exists <- file.exists(paste(outshape, c('shp', 'shx', 'dbf'), sep='.'))
    if (any(f.exists)) 
      stop(sprintf('File already exists: %s', 
                   toString(paste(outshape, c('shp', 'shx', 'dbf'), 
                                  sep='.')[f.exists])), call.=FALSE)
  } else outshape <- tempfile()
  if (is(x, 'Raster')) {
    require(raster)
    writeRaster(x, {f <- tempfile(fileext='.asc')})
    rastpath <- normalizePath(f)
  } else if (is.character(x)) {
    rastpath <- normalizePath(x)
  } else stop('x must be a file path (character string), or a Raster object.')
  
  ## Now 'python' has to be substituted by OSGeo4W
  #system2('python',
  system2('C:\\OSGeo4W64\\OSGeo4W.bat',
    args=(sprintf('"%1$s" "%2$s" -f "%3$s" "%4$s.shp"', 
    pypath, rastpath, gdalformat, outshape)))
  if (isTRUE(readpoly)) {
    shp <- readOGR(dirname(outshape), layer = basename(outshape), verbose=!quietish)
    return(shp) 
  }
  return(NULL)
}
