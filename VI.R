##################################################
## Vegetation indices for Multispectral data    
##################################################

## ------------------------------------------------------------------------------------------------------##
## Description: The script collect several vegetation indices for the analysis of multispectral data     ##
## B, G, R, NIR = blue, green, red and infra red channels                                                ##
## waveB, waveG, waveR, waveNIR = wavelength of the channels                                             ##
## ------------------------------------------------------------------------------------------------------##

msVI <- function (B, G, R, NIR, waveB, waveG, waveR, waveNIR){
  # Difference vegetation index
  DVI <- NIR - R

  # Simple ratio
  SR <- NIR/R

  # Modified simple ratio
  MSE <- ((NIR/R) - 1)/sqrt(1 + (NIR/R))

  # Normalized difference vegetation index
  NDVI <- ((NIR-R)/(NIR+R))

  # Green normalized difference vegetation index (GNDVI)
  GNDVI <- ((NIR-G)/(NIR+G))

  # Transformed vegetation Index
  TVI <- sqrt((NIR-R)/(NIR+R)) + 0.5

  # Soil adjusted vegetation index
  SAVI <- 1 + 0.5 * ((NIR-R)/(NIR+R+0.5))

  # Modified soil adjusted vegetation index 2 (MSAVI2)
  MSAVI2 <- (2 * NIR + 1 - sqrt(((2 * NIR + 1)^2) - 8 * (NIR * R)))/2

  # GEMI (Global Environment Monitoring Index)
  n <- (2*(NIR^2- R^2)+1.5*NIR+0.5*R)/(NIR+R+0.5)
  GEMI <- n*(1-0.25*n)*(R-0.125)/(1-R)

  # OSAVI (Optimized Soil Adjusted Vegetation Index)
  OSAVI <- (NIR - R)/(NIR + R + 0.16)

  # SRxNDVI (Simple Ratio x Normalized Difference Vegetation Index
  SRxNDVI <- (NIR^2 - R)/(NIR + R^2)

  # Renormalized difference vegetation index (RDVI)
  RDVI <- (NIR-R)/sqrt(NIR+R)

  # Angular Vegetation Index B-G-R
  v_a <- sqrt((waveB - waveG)^2+(B - G)^2)
  v_b <- sqrt((waveG - waveR)^2+(G - R)^2)
  v_c <- sqrt((waveB - waveR)^2+(B - R)^2)

  AVI_RGB1 <- acos((-v_c^2+v_a^2+v_b^2)/(2*v_a*v_b))
  AVI_RGB2 <- acos((-v_b^2+v_a^2+v_c^2)/(2*v_a*v_c))
  AVI_RGB3 <- acos((-v_a^2+v_b^2+v_c^2)/(2*v_b*v_c))

  AVI_RGB1[is.na(AVI_RGB1)] <- 0
  AVI_RGB2[is.na(AVI_RGB2)] <- 0
  AVI_RGB3[is.na(AVI_RGB3)] <- 0

  # Angular Vegetation Index G-R-NIR
  v_a <- sqrt((waveG - waveR)^2+(G - R)^2)
  v_b <- sqrt((waveR - waveNIR)^2+(R - NIR)^2)
  v_c <- sqrt((waveG - waveNIR)^2+(G - NIR)^2)

  AVI_GRNIR1 <- acos((-v_c^2+v_a^2+v_b^2)/(2*v_a*v_b))
  AVI_GRNIR2 <- acos((-v_b^2+v_a^2+v_c^2)/(2*v_a*v_c))
  AVI_GRNIR3 <- acos((-v_a^2+v_b^2+v_c^2)/(2*v_b*v_c))

  AVI_GRNIR1[is.na(AVI_GRNIR1)] <- 0
  AVI_GRNIR2[is.na(AVI_GRNIR2)] <- 0
  AVI_GRNIR3[is.na(AVI_GRNIR3)] <- 0

  # store results
  output <- data.frame(DVI, SR, MSE, NDVI, GNDVI, TVI, SAVI, MSAVI2, GEMI,  OSAVI, SRxNDVI,
                       RDVI, AVI_RGB1, AVI_RGB2, AVI_RGB3, AVI_GRNIR1, AVI_GRNIR2, AVI_GRNIR3)

  return(output)

}
