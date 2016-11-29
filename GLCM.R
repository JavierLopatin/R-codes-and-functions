##------------------------------------------------------------------------------##
##                                                                              ##
## GLCM: Apply Gray Level Covariance Matrix to an image (not used in the paper) ##
##                                                                              ##
## Arguments:                                                                   ##
## -  img       input image                                                     ##
##                                                                              ##
##------------------------------------------------------------------------------##

GLCM <- function(img){
  # Estimate the GLCM values of one band
  texture <- glcm(img, n_grey = 32, window = c(3, 3), shift = c(1, 1), 
                  statistics =c("homogeneity", "contrast", "dissimilarity", "entropy","second_moment", "correlation"),
                  min_x=NULL, max_x=NULL, na_opt="any",na_val=NA, scale_factor=1, asinteger=FALSE)
  stats <- c("homogeneity", "contrast", "dissimilarity", "entropy","second_moment", "correlation")
  Names <- paste(names(img), "_", stats, sep="")
  names(texture) <- Names
  return (texture)
}
