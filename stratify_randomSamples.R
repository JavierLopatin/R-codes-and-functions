################################################################
#
# Generate stratitied random samples inside a polygons
# according to their FID (cluster) Field
#
# @Author: Javier Lopatin | javier.lopatin@kit.edu
# usage:
#    shp = input shapefile
#    sample = number of samples to generate
#    minDistance = minimun distance betweem samples (default=60m)
#
# Example:
#    shapefile = readOGR('.', 'shape')
#    p1 = generate_samples(shapefile, 60)
#    p1 = SpatialPointsDataFrame(p1, data.frame(p1))
#    writeOGR(p1, ".", "samplePlots", driver="ESRI Shapefile"); rm(p1)
#
################################################################


# Function to generate stratify random samples
generate_samples <- function(shp, sample, minDistance=60){
  # shp = input shapefile
  # sample = number of samples to generate
  # minDistance = minimun distance betweem samples (default=60m)

  # load packages
  require(raster) # to use the function 'area'
  require(rgdal) # load and write shapefiles
  require(sp) # generate samples
  require(rgeos)

  # calculate area of polygons
  shp$area_ha = area(shp)/10000
  # proportion of size according to number of samples
  prop = ( shp$area_ha/sum(shp$area_ha) ) * sample
  prop = as.integer(prop)
  # if any number is zero, change it to one
  if (any(prop == 0)){
    idx = grep(0, prop)
    prop[idx] = 1
  }

  #### generate stratified random samples
  # regular grid of 60 m distance
  points = spsample(shp[1,], type='regular', cellsize=minDistance)
  # random selection of points
  points = sample(points, prop[1])
  strata = points; rm(points)

  # other clusters
  for (i in 2:length(shp$FID)){
    # regular grid of e.g. 60 m distance
    points = spsample(shp[i,], type='regular', cellsize=minDistance)
    # random selection of points
    points = sample(points, prop[i])
    # store
    strata = rbind(strata, points); rm(points)
  }
  return(strata)
}
