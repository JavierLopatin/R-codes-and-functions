## R-Script - Usefull functions for Remote Sensing data
## author: Javier Lopatin
## mail:javierlopatin@gmail.com & javier.lopatin@kit.edu
## last changes: 12/02/2016

##################################################################################
##                                                                              ##
## Multi-method ensemble selection of spectral bands                            ##
##                                                                              ##
## This function performs a band selection based on a multi-method ensemble     ##
## assessment of the variable importance and classification coefficients of     ##
## three different model types: Partial Least Squares Discriminant Analysis,    ##
## Random Forest and Support Vector Machine classifications                     ## 
##                                                                              ## 
## Arguments:                                                                   ##
## - x        Numeric matrix containing the spectra (samples as rows)           ##
## - y        Numeric vector containing the response variable                   ##
## - wl       Numeric vector containing the wavelength information of the bands ##
##                                                                              ##
## function based on the paper:                                                 ##
## Feilhauer, H., Asner, G. P., & Martin, R. E. (2015). Multi-method ensemble   ##
## selection of spectral bands related to leaf biochemistry. Remote Sensing of  ## 
## Environment, 164(November), 57-65. http://doi.org/10.1016/j.rse.2015.03.033  ##                                           ##
##                                                                              ##
##################################################################################

classificationEnsemble <- function(classes, spec, wl=NA){
  
  ## load required libraries
  library (caret)
  library(e1071)
  library(doParallel)
  
  # set data
  data <- data.frame(classes, spec)
  data <- na.omit(data)
  
  # Set the random number seed so we can reproduce the results
  set.seed(123)
  # Split data in training and test
  forTraining <- createDataPartition(data$classes, p = 0.6, list=F)
  train <- data [ forTraining,]
  test<- data [-forTraining,]
  
  # Each model used 5 repeated 10-fold cross-validation. Use AUC to pick the best model
  controlObject <- trainControl(method = "cv", number = 10,  repeats=10, classProbs=TRUE, allowParallel = TRUE)
  
  # initialize parallel processing
  cl <- makeCluster(detectCores())
  registerDoParallel(cl)
  
  #############################
  ### PLS-DA classification ###
  #############################
  
  # apply classification
  set.seed(123)
  plsClas <- train(x=train[, 2:length(train)], y=train$classes, method = "pls", tuneLength=20, 
                   preProc = c("center", "scale"), trControl = controlObject)
  
  # predict
  pls.pred <- predict(plsClas, test[, 2:length(train)])
  
  # confusion matix
  preMatrix <- function(pred, test){ # functionn to prevent caret error for different length
    u = union(pred, test)
    t = table(factor(pred, u), factor(test, u))
    return(t)
  }  
  conf.pls <- confusionMatrix(preMatrix(pls.pred, test$classes))
  
  # get accuracies
  PA.pls    <- conf.pls$byClass[,3] 
  UA.pls    <- conf.pls$byClass[,4]
  OA.pls    <- conf.pls$overall["Accuracy"]
  kappa.pls <- conf.pls$overall["Kappa"]
  
  ### variable importance
  plscf <- as.vector(rowMeans(plsClas$finalModel$coefficients)) ## extract coeff.
  plscf <- plscf / sd (plscf) ## scale regression coefficients
  
  #########################
  ### RF classification ###
  #########################
  
  set.seed(123)
  rfClas <- train(x=train[, 2:length(train)], y=train$classes, method = "rf", tuneLength=15, trControl = controlObject)
  
  # predict
  rf.pred <- predict(rfClas, test[,2:length(train)])
  
  # confusion matix
  conf.rf <- confusionMatrix(preMatrix(rf.pred, test$classes))
  
  # get accuracies
  PA.rf    <- conf.rf$byClass[,3] 
  UA.rf    <- conf.rf$byClass[,4]
  OA.rf    <- conf.rf$overall["Accuracy"]
  kappa.rf <- conf.rf$overall["Kappa"]
  
  ### variable importance
  rfcf <- varImp(rfClas$finalModel)[[1]]
  rfcf <- as.vector (rfcf / sd(rfcf))
  
  ##########################
  ### SVM classification ###
  ##########################
  
  set.seed(123)
  svmClas <- train(x=train[, 2:length(train)], y=train$classes, method = "svmLinear2", tuneLength=15, 
                   preProc = c("center", "scale"), trControl = controlObject)
  
  # predict
  svm.pred <- predict(svmClas, test[,2:length(train)])
  
  # confusion matix
  conf.svm <- confusionMatrix(preMatrix(svm.pred, test$classes))
  
  # get accuracies
  PA.svm    <- conf.svm$byClass[,3] 
  UA.svm    <- conf.svm$byClass[,4]
  OA.svm    <- conf.svm$overall["Accuracy"]
  kappa.svm <- conf.svm$overall["Kappa"]
  
  ### variable importance
  vr.alpha <- t(svmClas$finalModel$coefs) ## extract alpha vector
  svr.index <-  svmClas$finalModel$index ## extract alpha index
  ## calculate pseudo-regression coefficients from the alpha vector
  svrcf <- numeric (ncol (spec))
  for(i in 1:ncol(spec)){
    svrcf[i] <- svr.alpha %*% spec[svr.index, i]
  } 
  svrcf <- svrcf / sd (svrcf) ## scale pseudo-coefficients
  
  #####################################################################    
  ### get ensemble from all models and identify important variables ###
  #####################################################################
  
  ## get ensemble from all models and identify important variables
  ensemblecf <- abs(plscf) * OA.pls + abs(rfcf) * OA.rf + abs(svrcf) * OA.svm 
  th <- mean(ensemblecf) + sd(ensemblecf) ## calculate threshold
  r <- ensemblecf > th ## apply threshold
  
  # stop parallel process
  stopCluster(cl)
  
  ######################
  ### prepare output ###
  ######################
  
  cf <- rbind (wl, plscf, rfcf, svrcf, ensemblecf, selbands)
  colnames(cf) <- colnames(spec)
  
  fit <- c (OA.pls, OA.rf, OA.svm)
  names (fit) <- c ("PLS-DA OA", "RF OA", "SVR OA")
  output <- list (cf, fit, th, plsClas, rfClas, svmClas)
  names (output) <- c ("selection", "fits", "threshold", "PLS", "RF", "SVM")
  class (output) <- "ensemble"
  output
  
}

################################################################################
##                                                                            ##
## plot.ensemble: visualization of ensemble objects                           ##
##                                                                            ##
## Arguments:                                                                 ##
## - spec     spectral information. Used to create the quantiles of spectra   ##
## - en       classificationEnsemble object                                   ##
##                                                                            ##
################################################################################

plot.classificationEnsemble <- function (spec, en, label=TRUE, ...) {
  # extract the data from the classification Ensamble function
  wl <- en[[1]][1,]
  cf <- en[[1]][2:4,]
  fit <- en[[2]]
  cf <- cf * fit
  fit <- round (fit, 2)
  encf <- en[[1]][5,]
  z1 <- matrix (rep (encf, 100), ncol=100)
  sel <- as.logical (en[[1]][6,])
  z2 <- matrix (rep (sel, 100), ncol=100)
  z2[,11:100] <- NA
  z2[z2==0] <- NA
  # obtain the quantiles of the spectras
  quant <- apply(spec, 2, quantile, probs =c(0.05, 0.25, 0.5, 0.75, 0.95))
  # Plot the spectra
  par (mar=c (5, 4, 4, 7) + 0.1, xpd=NA)
  plot(wl, quant[1,], type="l", ylim = c(0,0.8), axes=F, ylab = NA, xlab=NA, las=1)
  lines(wl, quant[2,], type="l")
  lines(wl, quant[3,], type="l")
  lines(wl, quant[4,], type="l")
  lines(wl, quant[5,], type="l")
  polygon(c(wl, rev(wl)), c(quant[2,], rev(quant[1,])), col = "grey70")
  polygon(c(wl, rev(wl)), c(quant[3,], rev(quant[2,])), col = "grey50")
  polygon(c(wl, rev(wl)), c(quant[4,], rev(quant[3,])), col = "grey50")
  polygon(c(wl, rev(wl)), c(quant[5,], rev(quant[4,])), col = "grey70")
  axis(side = 4, las=1)
  mtext(side = 4, line = 3, 'Reflectance')
  # add coefficients
  par(new = T)
  plot(wl,  cf[1,], type = "l", col=2, xlab=expression(lambda(nm)), ylab="Weighted coefficients", las=1,
       ylim=c(min(cf), max(cf)))
  lines(wl,  cf[2,], type = "l", col=3)
  lines(wl,  cf[3,], type = "l", col=4) 
  # add spectral bands selected by the ensemble method
  image (wl, seq(min (cf) * 1.1, max (cf) * 1.1, length.out=100), z2, col=7, 
         xlab="", ylab="", add=T)
  # Labels
  if(label==TRUE){
    labels <- c (paste (c ("PLS ", "RF ", "SVM "), "OA", "=", fit, sep=""), NA, 
                 "Ensemble selection")
    legend ("topright", bty="n", col=c (2, 3, 4, NA, NA, rep (1, 3)), 
            pt.bg=c(rep (NA, 4), 7), lwd=c(rep (2, 3), rep (NA, 2)),
            pch=c (rep (NA, 4), 22), cex=0.7, pt.cex=1, legend=labels)
   }
 }


################################################################################
##                                                                            ##
## GLCM: Apply Gray Level Covariance Matrix to an image                       ##
##                                                                            ##
## Arguments:                                                                 ##
## -  img       input image                                                   ##
##                                                                            ##
################################################################################

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

################################################################################
##                                                                            ##
## Small functions to help in repetitive tasks                                ##
##                                                                            ##
################################################################################

## List the names of the rasters in a folder
rasterListNames <- function(fileExtantion, folder){
  # make a list of all fileExtantion files
  rast_list = list.files(folder, pattern = fileExtantion)
  # delete the ".dat" from the name
  rast_list = gsub('.{4}$', '', rast_list)
  return(rast_list)
}

## List and load the rasters contained in a folder
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


##################################################
## 1- Vegetation indices for Multispectral data ##
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

###################################################
## 2- Obtain the first derivative of the spectra ##
###################################################

## -------------------------------------------------##
## x = the spectral bands ordered in columns        ##   
## y = the band wavelength, e.g. seq(475, 875, 10)  ##   
## -------------------------------------------------##

first.deriv <- function(x, y){
  refl_2<-t(x) # transpose data
  refl_2<-data.frame(y, refl_2)
  N  = length(x)     # N° bands
  N2 = length(x[,1]) # observations
  B  = y[2] - y[1]   # band with
  der<-matrix(,N,N2) # sotored matrix
  
  for(i in 1:N){ 
    for(i2 in 2:N2){
      t1 <- refl_2[,i2]
      der[i,i2]<- (((t1[i+1]) - (t1[i])) / B) 
    }
  }
  der<-der[,-1] # delate column 1: only NA
  der_1<-data.frame(y,der)
  der_1<-na.omit(der_1)
}

####################################
## 3- plot all spectra signatures ##
####################################

## -------------------------------------------------##
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

#############################################
## 4- plot all first derivative signatures ##
#############################################

## -----------------------------------------------------##
## x = derivative results from the first.deriv function ##   
## -----------------------------------------------------##

plot.deriv <- function(x, ...){
  plot(x[,1], x[,2],type="l", axes=T, ylim=c(min(apply(x[, 2:length(x[1, ])],1,min)), max(apply(x[, 2:length(x[1, ])],1,max))), ...)
  for(i in 1:(length(x[1,])-1)){  lines(x[,1], x[, i+1])}
}

## ----------------------------- ##
##   Savitzky-Golay Algorithm
## ------------------------------##

# ------------------------------------------------------------------------------#
# Obtained from:https://stat.ethz.ch/pipermail/r-help/2004-February/045568.html #                                                                           #
# T2 <- sav.gol(T, fl, forder=4, dorder=0);                                     #
#                                                                               #
# Polynomial filtering method of Savitzky and Golay                             #
# See Numerical Recipes, 1992, Chapter 14.8, for details.                       #
#                                                                               #
# T      = vector of signals to be filtered                                     #
#          (the derivative is calculated for each ROW)                          #
# fl     = filter length ( must be odd. for instance fl = 51..151)              #
# forder = filter order (2 = quadratic filter, 4= quartic)                      #
# dorder = derivative order (0 = smoothing, 1 = first derivative, etc.)         #
# ------------------------------------------------------------------------------#

sav.gol <- function(T, fl, forder=4, dorder=0)
{
    m <- length(T)
    dorder <- dorder + 1

    # -- calculate filter coefficients --
    fc <- (fl-1)/2                          # index: window left and right
    X  <- outer(-fc:fc, 0:forder, FUN="^")  # polynomial terms and 
coefficients
    Y  <- pinv(X);                          # pseudoinverse

    # -- filter via convolution and take care of the end points --
    T2 <- convolve(T, rev(Y[dorder,]), type="o")    # convolve(...)
    T2 <- T2[(fc+1):(length(T2)-fc)]
}
#-----------------------------------------------------------------------
#   *** PseudoInvers of a Matrix ***
#   using singular value decomposition
#
pinv <- function (A)
{
    s <- svd(A)
    # D <- diag(s$d); Dinv <- diag(1/s$d)
    # U <- s$u; V <- s$v
    # A = U D V'
    # X = V Dinv U'
    s$v %*% diag(1/s$d) %*% t(s$u)
}
#-----------------------------------------------------------------------
