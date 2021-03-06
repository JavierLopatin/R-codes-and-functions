
################################################################################
##
## R-Script - Obtaning the continuous wavelet analysis variables from
## hyperspectral data for the reflectance and first derivative of the data
## author: Javier Lopatin
## mail: javierlopatin@gmail.com
## last changes: 18-12-2015

# data needed
# spectra = the spectral bands ordered in columns
# wavel = the band wavelength, e.g. seq(475, 875, 10)

##################
##### RESULTS ####
##################
# the results obtained with the script are:
#
# WAfit:          find the best family using autopls package. Input variable is the 
#                 observed variable only. e.g. WA_fit(Biomass)
# r2_all_ref:     list of all r2 of the reflectance data performed in WAfit
# rmse_all_ref:   list of all RMSE of the reflectance data performed in WAfit
# r2_all_deriv:   list of all r2 of the derivative data performed in WAfit
# rmse_all_deriv: list of all RMSE of the derivative data performed in WAfit
# bestWA:         see wich WA variables had the best fit. The input data are:
#                 WA_fit (r2_all_ref, rmse_all_ref, r2_all_deriv, rmse_all_deriv)
# plotWA:         plot the WA_fit results. The input data are the same that bestWA
#
################################################################################
# load required packages
library(Rwave)     #  For "Cauchy's" function
library(biwavelet) # For "morlet", "paul" and "dog" functions

# prepare paramiters
N  = length(spectra)     # N° band
N2 = length(spectra[,1]) # observations
B  = wavel[2] - wavel[1] # band with
refl_2<-t(spectra)       # transpose data
refl_2<-data.frame(wavel, refl_2)

## first derivative of the reflectance data
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

# Obtain the first derivative
der_1 <- first.deriv(spectra, wavel)

############################################
### Calculate wavelet for spectral input ###
############################################

## create empty matrix in which the transformations can be stored
# Morlet
morlet_refl_1<-data.frame(matrix(,N,N2))
morlet_refl_2<-data.frame(matrix(,N,N2))
morlet_refl_4<-data.frame(matrix(,N,N2))
morlet_refl_8<-data.frame(matrix(,N,N2))
morlet_refl_16<-data.frame(matrix(,N,N2))
morlet_refl_32<-data.frame(matrix(,N,N2))
# paul
paul_refl_1<-data.frame(matrix(,N,N2))
paul_refl_2<-data.frame(matrix(,N,N2))
paul_refl_4<-data.frame(matrix(,N,N2))
paul_refl_8<-data.frame(matrix(,N,N2))
paul_refl_16<-data.frame(matrix(,N,N2))
paul_refl_32<-data.frame(matrix(,N,N2))
# dog
dog_refl_1<-data.frame(matrix(,N,N2))
dog_refl_2<-data.frame(matrix(,N,N2))
dog_refl_4<-data.frame(matrix(,N,N2))
dog_refl_8<-data.frame(matrix(,N,N2))
dog_refl_16<-data.frame(matrix(,N,N2))
dog_refl_32<-data.frame(matrix(,N,N2))
# cauchy
cauchy_refl_1<-data.frame(matrix(,N,N2))
cauchy_refl_2<-data.frame(matrix(,N,N2))
cauchy_refl_4<-data.frame(matrix(,N,N2))
cauchy_refl_8<-data.frame(matrix(,N,N2))
cauchy_refl_16<-data.frame(matrix(,N,N2))
cauchy_refl_32<-data.frame(matrix(,N,N2))
# squiz
squiz_refl_1<-data.frame(matrix(,N,N2))
squiz_refl_2<-data.frame(matrix(,N,N2))
squiz_refl_4<-data.frame(matrix(,N,N2))
squiz_refl_8<-data.frame(matrix(,N,N2))
squiz_refl_16<-data.frame(matrix(,N,N2))
squiz_refl_32<-data.frame(matrix(,N,N2))
# gabor
gabor_refl_1<-data.frame(matrix(,N,N2))
gabor_refl_2<-data.frame(matrix(,N,N2))
gabor_refl_4<-data.frame(matrix(,N,N2))
gabor_refl_8<-data.frame(matrix(,N,N2))
gabor_refl_16<-data.frame(matrix(,N,N2))
gabor_refl_32<-data.frame(matrix(,N,N2))

for(i in 1:N2){
  # Morlet
  morlet<-wt(refl_2[,c(1, i+1)], mother="morlet")
  morlet_refl_1[[i]]<-  Re(morlet$wave[1,])  # Sclae 1
  morlet_refl_2[[i]]<-  Re(morlet$wave[2,])  # Scale 2
  morlet_refl_4[[i]]<-  Re(morlet$wave[4,])  # Scale 4
  morlet_refl_8[[i]]<-  Re(morlet$wave[8,])  # Scale 8
  morlet_refl_16[[i]]<- Re(morlet$wave[16,]) # Scale 16
  morlet_refl_32[[i]]<- Re(morlet$wave[32,]) # Scale 32
  # paul
  paul<-wt(refl_2[,c(1, i+1)], mother="paul")
  paul_refl_1[[i]]<-  Re(paul$wave[1,])
  paul_refl_2[[i]]<-  Re(paul$wave[2,])
  paul_refl_4[[i]]<-  Re(paul$wave[4,])
  paul_refl_8[[i]]<-  Re(paul$wave[8,])
  paul_refl_16[[i]]<- Re(paul$wave[16,])
  paul_refl_32[[i]]<- Re(paul$wave[32,])
  # dog
  dog<-wt(refl_2[,c(1, i+1)], mother="dog")
  dog_refl_1[[i]]<-  Re(dog$wave[1,])
  dog_refl_2[[i]]<-  Re(dog$wave[2,])
  dog_refl_4[[i]]<-  Re(dog$wave[4,])
  dog_refl_8[[i]]<-  Re(dog$wave[8,])
  dog_refl_16[[i]]<- Re(dog$wave[16,])
  dog_refl_32[[i]]<- Re(dog$wave[32,])
  # cauchy
  cauchy <-cwtTh(refl_2[,i+1],noctave=5, nvoice=12, moments=20)
  cauchy_refl_1[[i]]<-  Re(cauchy[,1])
  cauchy_refl_2[[i]]<-  Re(cauchy[,2])
  cauchy_refl_4[[i]]<-  Re(cauchy[,4])
  cauchy_refl_8[[i]]<-  Re(cauchy[,8])
  cauchy_refl_16[[i]]<- Re(cauchy[,16])
  cauchy_refl_32[[i]]<- Re(cauchy[,32])
  # squiz
  squiz<- cwtsquiz(refl_2[,i+1],noctave=5, nvoice=12)
  squiz_refl_1[[i]]<-  Re(squiz[,1])
  squiz_refl_2[[i]]<-  Re(squiz[,2])
  squiz_refl_4[[i]]<-  Re(squiz[,4])
  squiz_refl_8[[i]]<-  Re(squiz[,8])
  squiz_refl_16[[i]]<- Re(squiz[,16])
  squiz_refl_32[[i]]<- Re(squiz[,32])
  # gabor
  gabor<- cgt(refl_2[,i+1], nvoice=32)
  gabor_refl_1[[i]]<-  Re(gabor[,1])
  gabor_refl_2[[i]]<-  Re(gabor[,2])
  gabor_refl_4[[i]]<-  Re(gabor[,4])
  gabor_refl_8[[i]]<-  Re(gabor[,8])
  gabor_refl_16[[i]]<- Re(gabor[,16])
  gabor_refl_32[[i]]<- Re(gabor[,32])
}

#### Calculate wavelet for derivative input
# Morlet
morlet_deriv1_1<-data.frame(matrix(,N-1,N2-1))
morlet_deriv1_2<-data.frame(matrix(,N-1,N2-1))
morlet_deriv1_4<-data.frame(matrix(,N-1,N2-1))
morlet_deriv1_8<-data.frame(matrix(,N-1,N2-1))
morlet_deriv1_16<-data.frame(matrix(,N-1,N2-1))
morlet_deriv1_32<-data.frame(matrix(,N-1,N2-1))
# paul
paul_deriv1_1<-data.frame(matrix(,N-1,N2-1))
paul_deriv1_2<-data.frame(matrix(,N-1,N2-1))
paul_deriv1_4<-data.frame(matrix(,N-1,N2-1))
paul_deriv1_8<-data.frame(matrix(,N-1,N2-1))
paul_deriv1_16<-data.frame(matrix(,N-1,N2-1))
paul_deriv1_32<-data.frame(matrix(,N-1,N2-1))
# dog
dog_deriv1_1<-data.frame(matrix(,N-1,N2-1))
dog_deriv1_2<-data.frame(matrix(,N-1,N2-1))
dog_deriv1_4<-data.frame(matrix(,N-1,N2-1))
dog_deriv1_8<-data.frame(matrix(,N-1,N2-1))
dog_deriv1_16<-data.frame(matrix(,N-1,N2-1))
dog_deriv1_32<-data.frame(matrix(,N-1,N2-1))
# cauchy
cauchy_deriv1_1<-data.frame(matrix(,N-1,N2-1))
cauchy_deriv1_2<-data.frame(matrix(,N-1,N2-1))
cauchy_deriv1_4<-data.frame(matrix(,N-1,N2-1))
cauchy_deriv1_8<-data.frame(matrix(,N-1,N2-1))
cauchy_deriv1_16<-data.frame(matrix(,N-1,N2-1))
cauchy_deriv1_32<-data.frame(matrix(,N-1,N2-1))
# squiz
squiz_deriv1_1<-data.frame(matrix(,N-1,N2-1))
squiz_deriv1_2<-data.frame(matrix(,N-1,N2-1))
squiz_deriv1_4<-data.frame(matrix(,N-1,N2-1))
squiz_deriv1_8<-data.frame(matrix(,N-1,N2-1))
squiz_deriv1_16<-data.frame(matrix(,N-1,N2-1))
squiz_deriv1_32<-data.frame(matrix(,N-1,N2-1))
# gabor
gabor_deriv1_1<-data.frame(matrix(,N-1,N2-1))
gabor_deriv1_2<-data.frame(matrix(,N-1,N2-1))
gabor_deriv1_4<-data.frame(matrix(,N-1,N2-1))
gabor_deriv1_8<-data.frame(matrix(,N-1,N2-1))
gabor_deriv1_16<-data.frame(matrix(,N-1,N2-1))
gabor_deriv1_32<-data.frame(matrix(,N-1,N2-1))

for(i in 1:(N2-1)){
  # Morlet
  morlet<-wt(der_1[, c(1, i+1)], mother="morlet")
  morlet_deriv1_1[[i]]<-  Re( morlet$wave[1,])
  morlet_deriv1_2[[i]]<-  Re(morlet$wave[2,])
  morlet_deriv1_4[[i]]<-  Re(morlet$wave[4,])
  morlet_deriv1_8[[i]]<-  Re(morlet$wave[8,])
  morlet_deriv1_16[[i]]<- Re(morlet$wave[16,])
  morlet_deriv1_32[[i]]<- Re(morlet$wave[32,])
  # paul
  paul<-wt(der_1[,c(1, i+1)], mother="paul")
  paul_deriv1_1[[i]]<-  Re(paul$wave[1,])
  paul_deriv1_2[[i]]<-  Re(paul$wave[2,])
  paul_deriv1_4[[i]]<-  Re(paul$wave[4,])
  paul_deriv1_8[[i]]<-  Re(paul$wave[8,])
  paul_deriv1_16[[i]]<- Re(paul$wave[16,])
  paul_deriv1_32[[i]]<- Re(paul$wave[32,])
  # dog
  dog<-wt(der_1[,c(1, i+1)], mother="dog")
  dog_deriv1_1[[i]]<-  Re(dog$wave[1,])
  dog_deriv1_2[[i]]<-  Re(dog$wave[2,])
  dog_deriv1_4[[i]]<-  Re(dog$wave[4,])
  dog_deriv1_8[[i]]<-  Re(dog$wave[8,])
  dog_deriv1_16[[i]]<- Re(dog$wave[16,])
  dog_deriv1_32[[i]]<- Re(dog$wave[32,])
  # cauchy
  cauchy <-cwtTh(der_1[,i+1],noctave=5, nvoice=12, moments=20)
  cauchy_deriv1_1[[i]]<-  Re(cauchy[,1])
  cauchy_deriv1_2[[i]]<-  Re(cauchy[,2])
  cauchy_deriv1_4[[i]]<-  Re(cauchy[,4])
  cauchy_deriv1_8[[i]]<-  Re(cauchy[,8])
  cauchy_deriv1_16[[i]]<- Re(cauchy[,16])
  cauchy_deriv1_32[[i]]<- Re(cauchy[,32])
  # squiz
  squiz<- cwtsquiz(der_1[,i+1],noctave=5, nvoice=12)
  squiz_deriv1_1[[i]]<-  Re(squiz[,1])
  squiz_deriv1_2[[i]]<-  Re(squiz[,2])
  squiz_deriv1_4[[i]]<-  Re(squiz[,4])
  squiz_deriv1_8[[i]]<-  Re(squiz[,8])
  squiz_deriv1_16[[i]]<- Re(squiz[,16])
  squiz_deriv1_32[[i]]<- Re(squiz[,32])
  # gabor
  gabor<- cgt(der_1[,i+1], nvoice=32)
  gabor_deriv1_1[[i]]<-  Re(gabor[,1])
  gabor_deriv1_2[[i]]<-  Re(gabor[,2])
  gabor_deriv1_4[[i]]<-  Re(gabor[,4])
  gabor_deriv1_8[[i]]<-  Re(gabor[,8])
  gabor_deriv1_16[[i]]<- Re(gabor[,16])
  gabor_deriv1_32[[i]]<- Re(gabor[,32])

}

# the list of WA outputs
WA_names <- c("morlet_refl_1", "morlet_refl_2", "morlet_refl_4", "morlet_refl_8", "morlet_refl_16", "morlet_refl_32",
             "morlet_deriv1_1", "morlet_deriv1_2", "morlet_deriv1_4", "morlet_deriv1_8", "morlet_deriv1_16", "morlet_deriv1_32",
             "paul_refl_1", "paul_refl_2", "paul_refl_4", "paul_refl_8", "paul_refl_16", "paul_refl_32",
             "paul_deriv1_1", "paul_deriv1_2", "paul_deriv1_4", "paul_deriv1_8", "paul_deriv1_16", "paul_deriv1_32",
             "dog_refl_1", "dog_refl_2", "dog_refl_4", "dog_refl_8", "dog_refl_16", "dog_refl_32",
             "dog_deriv1_1", "dog_deriv1_2", "dog_deriv1_4", "dog_deriv1_8", "dog_deriv1_16", "dog_deriv1_32",
             "cauchy_refl_1", "cauchy_refl_2", "cauchy_refl_4", "cauchy_refl_8", "cauchy_refl_16", "cauchy_refl_32",
             "cauchy_deriv1_1", "cauchy_deriv1_2", "cauchy_deriv1_4", "cauchy_deriv1_8", "cauchy_deriv1_16", "cauchy_deriv1_32",
             "squiz_refl_1", "squiz_refl_2", "squiz_refl_4", "squiz_refl_8", "squiz_refl_16", "squiz_refl_32",
             "squiz_deriv1_1", "squiz_deriv1_2", "squiz_deriv1_4", "squiz_deriv1_8", "squiz_deriv1_16", "squiz_deriv1_32",
             "gabor_refl_1", "gabor_refl_2", "gabor_refl_4", "gabor_refl_8", "gabor_refl_16", "gabor_refl_32",
             "gabor_deriv1_1", "gabor_deriv1_2", "gabor_deriv1_4", "gabor_deriv1_8", "gabor_deriv1_16", "gabor_deriv1_32")

WA_values <- list(morlet_refl_1, morlet_refl_2, morlet_refl_4, morlet_refl_8, morlet_refl_16, morlet_refl_32,
              morlet_deriv1_1, morlet_deriv1_2, morlet_deriv1_4, morlet_deriv1_8, morlet_deriv1_16, morlet_deriv1_32,
              paul_refl_1, paul_refl_2, paul_refl_4, paul_refl_8, paul_refl_16, paul_refl_32,
              paul_deriv1_1, paul_deriv1_2, paul_deriv1_4, paul_deriv1_8, paul_deriv1_16, paul_deriv1_32,
              dog_refl_1, dog_refl_2, dog_refl_4, dog_refl_8, dog_refl_16, dog_refl_32,
              dog_deriv1_1, dog_deriv1_2, dog_deriv1_4, dog_deriv1_8, dog_deriv1_16, dog_deriv1_32,
              cauchy_refl_1, cauchy_refl_2, cauchy_refl_4, cauchy_refl_8, cauchy_refl_16, cauchy_refl_32,
              cauchy_deriv1_1, cauchy_deriv1_2, cauchy_deriv1_4, cauchy_deriv1_8, cauchy_deriv1_16, cauchy_deriv1_32,
              squiz_refl_1, squiz_refl_2, squiz_refl_4, squiz_refl_8, squiz_refl_16, squiz_refl_32,
              squiz_deriv1_1, squiz_deriv1_2, squiz_deriv1_4, squiz_deriv1_8, squiz_deriv1_16, squiz_deriv1_32,
              gabor_refl_1, gabor_refl_2, gabor_refl_4, gabor_refl_8, gabor_refl_16, gabor_refl_32,
              gabor_deriv1_1, gabor_deriv1_2, gabor_deriv1_4, gabor_deriv1_8, gabor_deriv1_16, gabor_deriv1_32)

# trampose data frames
for (i in 1:length(WA_names)){
  r <- WA_values[[i]];   r <- as.data.frame(t(r))
  assign(WA_names[[i]], r)
}
rm(r) # delate the list of data frames


###############################
### Find the best WA family ###
###############################

# function where x is the predictor variable
WAfit <- function(x){

library(autopls)

### Reflectance input ###
#
###### Morlet
#
# scale 1
morlet_refl_1 <- data.frame(x, (morlet_refl_1))
fit<-autopls(x ~.,  data=morlet_refl_1)
r2_morlet_refl_1 <- R2(fit)
rmse_morlet_refl_1 <- RMSEP(fit, estimate="CV")
# scale 2
morlet_refl_2 <- data.frame(x, (morlet_refl_2))
fit<-autopls(x ~.,  data=morlet_refl_2)
r2_morlet_refl_2 <- R2(fit)
rmse_morlet_refl_2 <- RMSEP(fit, estimate="CV")
# scale 4
morlet_refl_4 <- data.frame(x, (morlet_refl_4))
fit<-autopls(x ~.,  data=morlet_refl_4)
r2_morlet_refl_4 <- R2(fit)
rmse_morlet_refl_4 <- RMSEP(fit, estimate="CV")
# scale 8
morlet_refl_8 <- data.frame(x, (morlet_refl_8))
fit<-autopls(x ~.,  data=morlet_refl_8)
r2_morlet_refl_8 <- R2(fit)
rmse_morlet_refl_8 <- RMSEP(fit, estimate="CV")
# scale 16
morlet_refl_16 <- data.frame(x, (morlet_refl_16))
fit<-autopls(x ~.,  data=morlet_refl_16)
r2_morlet_refl_16 <- R2(fit)
rmse_morlet_refl_16 <- RMSEP(fit, estimate="CV")
# scale 32
morlet_refl_32 <- data.frame(x, (morlet_refl_32))
fit<-autopls(x ~.,  data=morlet_refl_32)
r2_morlet_refl_32 <- R2(fit)
rmse_morlet_refl_32 <- RMSEP(fit, estimate="CV")
#
###### paul
#
# scale 1
paul_refl_1 <- data.frame(x, (paul_refl_1))
fit<-autopls(x ~.,  data=paul_refl_1)
r2_paul_refl_1 <- R2(fit)
rmse_paul_refl_1 <- RMSEP(fit, estimate="CV")
# scale 2
paul_refl_2 <- data.frame(x, (paul_refl_2))
fit<-autopls(x ~.,  data=paul_refl_2)
r2_paul_refl_2 <- R2(fit)
rmse_paul_refl_2 <- RMSEP(fit, estimate="CV")
# scale 4
paul_refl_4 <- data.frame(x, (paul_refl_4))
fit<-autopls(x ~.,  data=paul_refl_4)
r2_paul_refl_4 <- R2(fit)
rmse_paul_refl_4 <- RMSEP(fit, estimate="CV")
# scale 8
paul_refl_8 <- data.frame(x, (paul_refl_8))
fit<-autopls(x ~.,  data=paul_refl_8)
r2_paul_refl_8 <- R2(fit)
rmse_paul_refl_8 <- RMSEP(fit, estimate="CV")
# scale 16
paul_refl_16 <- data.frame(x, (paul_refl_16))
fit<-autopls(x ~.,  data=paul_refl_16)
r2_paul_refl_16 <- R2(fit)
rmse_paul_refl_16 <- RMSEP(fit, estimate="CV")
# scale 32
paul_refl_32 <- data.frame(x, (paul_refl_32))
fit<-autopls(x ~.,  data=paul_refl_32)
r2_paul_refl_32 <- R2(fit)
rmse_paul_refl_32 <- RMSEP(fit, estimate="CV")
#
###### dog
#
# scale 1
dog_refl_1 <- data.frame(x, (dog_refl_1))
fit<-autopls(x ~.,  data=dog_refl_1)
r2_dog_refl_1 <- R2(fit)
rmse_dog_refl_1 <- RMSEP(fit, estimate="CV")
# scale 2
dog_refl_2 <- data.frame(x, (dog_refl_2))
fit<-autopls(x ~.,  data=dog_refl_2)
r2_dog_refl_2 <- R2(fit)
rmse_dog_refl_2 <- RMSEP(fit, estimate="CV")
# scale 4
dog_refl_4 <- data.frame(x, (dog_refl_4))
fit<-autopls(x ~.,  data=dog_refl_4)
r2_dog_refl_4 <- R2(fit)
rmse_dog_refl_4 <- RMSEP(fit, estimate="CV")
# scale 8
dog_refl_8 <- data.frame(x, (dog_refl_8))
fit<-autopls(x ~.,  data=dog_refl_8)
r2_dog_refl_8 <- R2(fit)
rmse_dog_refl_8 <- RMSEP(fit, estimate="CV")
# scale 16
dog_refl_16 <- data.frame(x, (dog_refl_16))
fit<-autopls(x ~.,  data=dog_refl_16)
r2_dog_refl_16 <- R2(fit)
rmse_dog_refl_16 <- RMSEP(fit, estimate="CV")
# scale 32
dog_refl_32 <- data.frame(x, (dog_refl_32))
fit<-autopls(x ~.,  data=dog_refl_32)
r2_dog_refl_32 <- R2(fit)
rmse_dog_refl_32 <- RMSEP(fit, estimate="CV")
#
###### cauchy
#
# scale 1
cauchy_refl_1 <- data.frame(x, (cauchy_refl_1))
fit<-autopls(x ~.,  data=cauchy_refl_1)
r2_cauchy_refl_1 <- R2(fit)
rmse_cauchy_refl_1 <- RMSEP(fit, estimate="CV")
# scale 2
cauchy_refl_2 <- data.frame(x, (cauchy_refl_2))
fit<-autopls(x ~.,  data=cauchy_refl_2)
r2_cauchy_refl_2 <- R2(fit)
rmse_cauchy_refl_2 <- RMSEP(fit, estimate="CV")
# scale 4
cauchy_refl_4 <- data.frame(x, (cauchy_refl_4))
fit<-autopls(x ~.,  data=cauchy_refl_4)
r2_cauchy_refl_4 <- R2(fit)
rmse_cauchy_refl_4 <- RMSEP(fit, estimate="CV")
# scale 8
cauchy_refl_8 <- data.frame(x, (cauchy_refl_8))
fit<-autopls(x ~.,  data=cauchy_refl_8)
r2_cauchy_refl_8 <- R2(fit)
rmse_cauchy_refl_8 <- RMSEP(fit, estimate="CV")
# scale 16
cauchy_refl_16 <- data.frame(x, (cauchy_refl_16))
fit<-autopls(x ~.,  data=cauchy_refl_16)
r2_cauchy_refl_16 <- R2(fit)
rmse_cauchy_refl_16 <- RMSEP(fit, estimate="CV")
# scale 32
cauchy_refl_32 <- data.frame(x, (cauchy_refl_32))
fit<-autopls(x ~.,  data=cauchy_refl_32)
r2_cauchy_refl_32 <- R2(fit)
rmse_cauchy_refl_32 <- RMSEP(fit, estimate="CV")
#
###### squiz
#
# scale 1
squiz_refl_1 <- data.frame(x, (squiz_refl_1))
fit<-autopls(x ~.,  data=squiz_refl_1)
r2_squiz_refl_1 <- R2(fit)
rmse_squiz_refl_1 <- RMSEP(fit, estimate="CV")
# scale 2
squiz_refl_2 <- data.frame(x, (squiz_refl_2))
fit<-autopls(x ~.,  data=squiz_refl_2)
r2_squiz_refl_2 <- R2(fit)
rmse_squiz_refl_2 <- RMSEP(fit, estimate="CV")
# scale 4
squiz_refl_4 <- data.frame(x, (squiz_refl_4))
fit<-autopls(x ~.,  data=squiz_refl_4)
r2_squiz_refl_4 <- R2(fit)
rmse_squiz_refl_4 <- RMSEP(fit, estimate="CV")
# scale 8
squiz_refl_8 <- data.frame(x, (squiz_refl_8))
fit<-autopls(x ~.,  data=squiz_refl_8)
r2_squiz_refl_8 <- R2(fit)
rmse_squiz_refl_8 <- RMSEP(fit, estimate="CV")
# scale 16
squiz_refl_16 <- data.frame(x, (squiz_refl_16))
fit<-autopls(x ~.,  data=squiz_refl_16)
r2_squiz_refl_16 <- R2(fit)
rmse_squiz_refl_16 <- RMSEP(fit, estimate="CV")
# scale 32
squiz_refl_32 <- data.frame(x, (squiz_refl_32))
fit<-autopls(x ~.,  data=squiz_refl_32)
r2_squiz_refl_32 <- R2(fit)
rmse_squiz_refl_32 <- RMSEP(fit, estimate="CV")
#
###### gabor
#
# scale 1
gabor_refl_1 <- data.frame(x, (gabor_refl_1))
fit<-autopls(x ~.,  data=gabor_refl_1)
r2_gabor_refl_1 <- R2(fit)
rmse_gabor_refl_1 <- RMSEP(fit, estimate="CV")
# scale 2
gabor_refl_2 <- data.frame(x, (gabor_refl_2))
fit<-autopls(x ~.,  data=gabor_refl_2)
r2_gabor_refl_2 <- R2(fit)
rmse_gabor_refl_2 <- RMSEP(fit, estimate="CV")
# scale 4
gabor_refl_4 <- data.frame(x, (gabor_refl_4))
fit<-autopls(x ~.,  data=gabor_refl_4)
r2_gabor_refl_4 <- R2(fit)
rmse_gabor_refl_4 <- RMSEP(fit, estimate="CV")
# scale 8
gabor_refl_8 <- data.frame(x, (gabor_refl_8))
fit<-autopls(x ~.,  data=gabor_refl_8)
r2_gabor_refl_8 <- R2(fit)
rmse_gabor_refl_8 <- RMSEP(fit, estimate="CV")
# scale 16
gabor_refl_16 <- data.frame(x, (gabor_refl_16))
fit<-autopls(x ~.,  data=gabor_refl_16)
r2_gabor_refl_16 <- R2(fit)
rmse_gabor_refl_16 <- RMSEP(fit, estimate="CV")
# scale 32
gabor_refl_32 <- data.frame(x, (gabor_refl_32))
fit<-autopls(x ~.,  data=gabor_refl_32)
r2_gabor_refl_32 <- R2(fit)
rmse_gabor_refl_32 <- RMSEP(fit, estimate="CV")

### 1rst derivative input ###

#
###### Morlet
#
# scale 1
morlet_deriv1_1 <- data.frame(x[1:length(x) - 1], (morlet_deriv1_1))
colnames(morlet_deriv1_1)[1] <-c("x")
fit<-autopls(x ~.,  data=morlet_deriv1_1)
r2_morlet_deriv1_1 <- R2(fit)
rmse_morlet_deriv1_1 <- RMSEP(fit, estimate="CV")
# scale 2
morlet_deriv1_2 <- data.frame(x[1:length(x) - 1], (morlet_deriv1_2))
colnames(morlet_deriv1_2)[1] <-c("x")
fit<-autopls(x ~.,  data=morlet_deriv1_2)
r2_morlet_deriv1_2 <- R2(fit)
rmse_morlet_deriv1_2 <- RMSEP(fit, estimate="CV")
# scale 4
morlet_deriv1_4 <- data.frame(x[1:length(x) - 1], (morlet_deriv1_4))
colnames(morlet_deriv1_4)[1] <-c("x")
fit<-autopls(x ~.,  data=morlet_deriv1_4)
r2_morlet_deriv1_4 <- R2(fit)
rmse_morlet_deriv1_4 <- RMSEP(fit, estimate="CV")
# scale 8
morlet_deriv1_8 <- data.frame(x[1:length(x) - 1], (morlet_deriv1_8))
colnames(morlet_deriv1_8)[1] <-c("x")
fit<-autopls(x ~.,  data=morlet_deriv1_8)
r2_morlet_deriv1_8 <- R2(fit)
rmse_morlet_deriv1_8 <- RMSEP(fit, estimate="CV")
# scale 16
morlet_deriv1_16 <- data.frame(x[1:length(x) - 1], (morlet_deriv1_16))
colnames(morlet_deriv1_16)[1] <-c("x")
fit<-autopls(x ~.,  data=morlet_deriv1_16)
r2_morlet_deriv1_16 <- R2(fit)
rmse_morlet_deriv1_16 <- RMSEP(fit, estimate="CV")
# scale 32
morlet_deriv1_32 <- data.frame(x[1:length(x) - 1], (morlet_deriv1_32))
colnames(morlet_deriv1_32)[1] <-c("x")
fit<-autopls(x ~.,  data=morlet_deriv1_32)
r2_morlet_deriv1_32 <- R2(fit)
rmse_morlet_deriv1_32 <- RMSEP(fit, estimate="CV")
#
###### paul
#
# scale 1
paul_deriv1_1 <- data.frame(x[1:length(x) - 1], (paul_deriv1_1))
colnames(paul_deriv1_1)[1] <-c("x")
fit<-autopls(x ~.,  data=paul_deriv1_1)
r2_paul_deriv1_1 <- R2(fit)
rmse_paul_deriv1_1 <- RMSEP(fit, estimate="CV")
# scale 2
paul_deriv1_2 <- data.frame(x[1:length(x) - 1], (paul_deriv1_2))
colnames(paul_deriv1_2)[1] <-c("x")
fit<-autopls(x ~.,  data=paul_deriv1_2)
r2_paul_deriv1_2 <- R2(fit)
rmse_paul_deriv1_2 <- RMSEP(fit, estimate="CV")
# scale 4
paul_deriv1_4 <- data.frame(x[1:length(x) - 1], (paul_deriv1_4))
colnames(paul_deriv1_4)[1] <-c("x")
fit<-autopls(x ~.,  data=paul_deriv1_4)
r2_paul_deriv1_4 <- R2(fit)
rmse_paul_deriv1_4 <- RMSEP(fit, estimate="CV")
# scale 8
paul_deriv1_8 <- data.frame(x[1:length(x) - 1], (paul_deriv1_8))
colnames(paul_deriv1_8)[1] <-c("x")
fit<-autopls(x ~.,  data=paul_deriv1_8)
r2_paul_deriv1_8 <- R2(fit)
rmse_paul_deriv1_8 <- RMSEP(fit, estimate="CV")
# scale 16
paul_deriv1_16 <- data.frame(x[1:length(x) - 1], (paul_deriv1_16))
colnames(paul_deriv1_16)[1] <-c("x")
fit<-autopls(x ~.,  data=paul_deriv1_16)
r2_paul_deriv1_16 <- R2(fit)
rmse_paul_deriv1_16 <- RMSEP(fit, estimate="CV")
# scale 32
paul_deriv1_32 <- data.frame(x[1:length(x) - 1], (paul_deriv1_32))
colnames(paul_deriv1_32)[1] <-c("x")
fit<-autopls(x ~.,  data=paul_deriv1_32)
r2_paul_deriv1_32 <- R2(fit)
rmse_paul_deriv1_32 <- RMSEP(fit, estimate="CV")
#
###### dog
#
# scale 1
dog_deriv1_1 <- data.frame(x[1:length(x) - 1], (dog_deriv1_1))
colnames(dog_deriv1_1)[1] <-c("x")
fit<-autopls(x ~.,  data=dog_deriv1_1)
r2_dog_deriv1_1 <- R2(fit)
rmse_dog_deriv1_1 <- RMSEP(fit, estimate="CV")
# scale 2
dog_deriv1_2 <- data.frame(x[1:length(x) - 1], (dog_deriv1_2))
colnames(dog_deriv1_2)[1] <-c("x")
fit<-autopls(x ~.,  data=dog_deriv1_2)
r2_dog_deriv1_2 <- R2(fit)
rmse_dog_deriv1_2 <- RMSEP(fit, estimate="CV")
# scale 4
dog_deriv1_4 <- data.frame(x[1:length(x) - 1], (dog_deriv1_4))
colnames(dog_deriv1_4)[1] <-c("x")
fit<-autopls(x ~.,  data=dog_deriv1_4)
r2_dog_deriv1_4 <- R2(fit)
rmse_dog_deriv1_4 <- RMSEP(fit, estimate="CV")
# scale 8
dog_deriv1_8 <- data.frame(x[1:length(x) - 1], (dog_deriv1_8))
colnames(dog_deriv1_8)[1] <-c("x")
fit<-autopls(x ~.,  data=dog_deriv1_8)
r2_dog_deriv1_8 <- R2(fit)
rmse_dog_deriv1_8 <- RMSEP(fit, estimate="CV")
# scale 16
dog_deriv1_16 <- data.frame(x[1:length(x) - 1], (dog_deriv1_16))
colnames(dog_deriv1_16)[1] <-c("x")
fit<-autopls(x ~.,  data=dog_deriv1_16)
r2_dog_deriv1_16 <- R2(fit)
rmse_dog_deriv1_16 <- RMSEP(fit, estimate="CV")
# scale 32
dog_deriv1_32 <- data.frame(x[1:length(x) - 1], (dog_deriv1_32))
colnames(dog_deriv1_32)[1] <-c("x")
fit<-autopls(x ~.,  data=dog_deriv1_32)
r2_dog_deriv1_32 <- R2(fit)
rmse_dog_deriv1_32 <- RMSEP(fit, estimate="CV")
#
###### cauchy
#
# scale 1
cauchy_deriv1_1 <- data.frame(x[1:length(x) - 1], (cauchy_deriv1_1))
colnames(cauchy_deriv1_1)[1] <-c("x")
fit<-autopls(x ~.,  data=cauchy_deriv1_1)
r2_cauchy_deriv1_1 <- R2(fit)
rmse_cauchy_deriv1_1 <- RMSEP(fit, estimate="CV")
# scale 2
cauchy_deriv1_2 <- data.frame(x[1:length(x) - 1], (cauchy_deriv1_2))
colnames(cauchy_deriv1_2)[1] <-c("x")
fit<-autopls(x ~.,  data=cauchy_deriv1_2)
r2_cauchy_deriv1_2 <- R2(fit)
rmse_cauchy_deriv1_2 <- RMSEP(fit, estimate="CV")
# scale 4
cauchy_deriv1_4 <- data.frame(x[1:length(x) - 1], (cauchy_deriv1_4))
colnames(cauchy_deriv1_4)[1] <-c("x")
fit<-autopls(x ~.,  data=cauchy_deriv1_4)
r2_cauchy_deriv1_4 <- R2(fit)
rmse_cauchy_deriv1_4 <- RMSEP(fit, estimate="CV")
# scale 8
cauchy_deriv1_8 <- data.frame(x[1:length(x) - 1], (cauchy_deriv1_8))
colnames(cauchy_deriv1_8)[1] <-c("x")
fit<-autopls(x ~.,  data=cauchy_deriv1_8)
r2_cauchy_deriv1_8 <- R2(fit)
rmse_cauchy_deriv1_8 <- RMSEP(fit, estimate="CV")
# scale 16
cauchy_deriv1_16 <- data.frame(x[1:length(x) - 1], (cauchy_deriv1_16))
colnames(cauchy_deriv1_16)[1] <-c("x")
fit<-autopls(x ~.,  data=cauchy_deriv1_16)
r2_cauchy_deriv1_16 <- R2(fit)
rmse_cauchy_deriv1_16 <- RMSEP(fit, estimate="CV")
# scale 32
cauchy_deriv1_32 <- data.frame(x[1:length(x) - 1], (cauchy_deriv1_32))
colnames(cauchy_deriv1_32)[1] <-c("x")
fit<-autopls(x ~.,  data=cauchy_deriv1_32)
r2_cauchy_deriv1_32 <- R2(fit)
rmse_cauchy_deriv1_32 <- RMSEP(fit, estimate="CV")
#
###### squiz
#
# scale 1
squiz_deriv1_1 <- data.frame(x[1:length(x) - 1], (squiz_deriv1_1))
colnames(squiz_deriv1_1)[1] <-c("x")
fit<-autopls(x ~.,  data=squiz_deriv1_1)
r2_squiz_deriv1_1 <- R2(fit)
rmse_squiz_deriv1_1 <- RMSEP(fit, estimate="CV")
# scale 2
squiz_deriv1_2 <- data.frame(x[1:length(x) - 1], (squiz_deriv1_2))
colnames(squiz_deriv1_2)[1] <-c("x")
fit<-autopls(x ~.,  data=squiz_deriv1_2)
r2_squiz_deriv1_2 <- R2(fit)
rmse_squiz_deriv1_2 <- RMSEP(fit, estimate="CV")
# scale 4
squiz_deriv1_4 <- data.frame(x[1:length(x) - 1], (squiz_deriv1_4))
colnames(squiz_deriv1_4)[1] <-c("x")
fit<-autopls(x ~.,  data=squiz_deriv1_4)
r2_squiz_deriv1_4 <- R2(fit)
rmse_squiz_deriv1_4 <- RMSEP(fit, estimate="CV")
# scale 8
squiz_deriv1_8 <- data.frame(x[1:length(x) - 1], (squiz_deriv1_8))
colnames(squiz_deriv1_8)[1] <-c("x")
fit<-autopls(x ~.,  data=squiz_deriv1_8)
r2_squiz_deriv1_8 <- R2(fit)
rmse_squiz_deriv1_8 <- RMSEP(fit, estimate="CV")
# scale 16
squiz_deriv1_16 <- data.frame(x[1:length(x) - 1], (squiz_deriv1_16))
colnames(squiz_deriv1_16)[1] <-c("x")
fit<-autopls(x ~.,  data=squiz_deriv1_16)
r2_squiz_deriv1_16 <- R2(fit)
rmse_squiz_deriv1_16 <- RMSEP(fit, estimate="CV")
# scale 32
squiz_deriv1_32 <- data.frame(x[1:length(x) - 1], (squiz_deriv1_32))
colnames(squiz_deriv1_32)[1] <-c("x")
fit<-autopls(x ~.,  data=squiz_deriv1_32)
r2_squiz_deriv1_32 <- R2(fit)
rmse_squiz_deriv1_32 <- RMSEP(fit, estimate="CV")
#
###### gabor
#
# scale 1
gabor_deriv1_1 <- data.frame(x[1:length(x) - 1], (gabor_deriv1_1))
colnames(gabor_deriv1_1)[1] <-c("x")
fit<-autopls(x ~.,  data=gabor_deriv1_1)
r2_gabor_deriv1_1 <- R2(fit)
rmse_gabor_deriv1_1 <- RMSEP(fit, estimate="CV")
# scale 2
gabor_deriv1_2 <- data.frame(x[1:length(x) - 1], (gabor_deriv1_2))
colnames(gabor_deriv1_2)[1] <-c("x")
fit<-autopls(x ~.,  data=gabor_deriv1_2)
r2_gabor_deriv1_2 <- R2(fit)
rmse_gabor_deriv1_2 <- RMSEP(fit, estimate="CV")
# scale 4
gabor_deriv1_4 <- data.frame(x[1:length(x) - 1], (gabor_deriv1_4))
colnames(gabor_deriv1_4)[1] <-c("x")
fit<-autopls(x ~.,  data=gabor_deriv1_4)
r2_gabor_deriv1_4 <- R2(fit)
rmse_gabor_deriv1_4 <- RMSEP(fit, estimate="CV")
# scale 8
gabor_deriv1_8 <- data.frame(x[1:length(x) - 1], (gabor_deriv1_8))
colnames(gabor_deriv1_8)[1] <-c("x")
fit<-autopls(x ~.,  data=gabor_deriv1_8)
r2_gabor_deriv1_8 <- R2(fit)
rmse_gabor_deriv1_8 <- RMSEP(fit, estimate="CV")
# scale 16
gabor_deriv1_16 <- data.frame(x[1:length(x) - 1], (gabor_deriv1_16))
colnames(gabor_deriv1_16)[1] <-c("x")
fit<-autopls(x ~.,  data=gabor_deriv1_16)
r2_gabor_deriv1_16 <- R2(fit)
rmse_gabor_deriv1_16 <- RMSEP(fit, estimate="CV")
# scale 32
gabor_deriv1_32 <- data.frame(x[1:length(x) - 1], (gabor_deriv1_32))
colnames(gabor_deriv1_32)[1] <-c("x")
fit<-autopls(x ~.,  data=gabor_deriv1_32)
r2_gabor_deriv1_32 <- R2(fit)
rmse_gabor_deriv1_32 <- RMSEP(fit, estimate="CV")


## Store results
# reflance imput
# r2
r2_all_ref <- data.frame(r2_morlet_refl_1$val[1], r2_morlet_refl_2$val[1], r2_morlet_refl_4$val[1], r2_morlet_refl_8$val[1], r2_morlet_refl_16$val[1], r2_morlet_refl_32$val[1],
                         r2_paul_refl_1$val[1], r2_paul_refl_2$val[1], r2_paul_refl_4$val[1], r2_paul_refl_8$val[1], r2_paul_refl_16$val[1], r2_paul_refl_32$val[1],
                         r2_dog_refl_1$val[1], r2_dog_refl_2$val[1], r2_dog_refl_4$val[1], r2_dog_refl_8$val[1], r2_dog_refl_16$val[1], r2_dog_refl_32$val[1],
                         r2_cauchy_refl_1$val[1], r2_cauchy_refl_2$val[1], r2_cauchy_refl_4$val[1], r2_cauchy_refl_8$val[1], r2_cauchy_refl_16$val[1], r2_cauchy_refl_32$val[1],
                         r2_squiz_refl_1$val[1], r2_squiz_refl_2$val[1], r2_squiz_refl_4$val[1], r2_squiz_refl_8$val[1], r2_squiz_refl_16$val[1], r2_squiz_refl_32$val[1],
                         r2_gabor_refl_1$val[1], r2_gabor_refl_2$val[1], r2_gabor_refl_4$val[1], r2_gabor_refl_8$val[1], r2_gabor_refl_16$val[1], r2_gabor_refl_32$val[1])
names(r2_all_ref) <- c("morlet_1", "morlet_2", "morlet_4", "morlet_8", "morlet_16", "morlet_32",
                       "paul_1", "paul_2", "paul_4", "paul_8", "paul_16", "paul_32",
                       "dog_1", "dog_2", "dog_4", "dog_8", "dog_16", "dog_32",
                       "cauchy_1", "cauchy_2", "cauchy_4", "cauchy_8", "cauchy_16", "cauchy_32",
                       "squiz_1", "squiz_2", "squiz_4", "squiz_8", "squiz_16", "squiz_32",
                       "gabor_1", "gabor_2", "gabor_4", "gabor_8", "gabor_16", "gabor_32")
abs(r2_all_ref)

# RMSE
rmse_all_ref <- data.frame(rmse_morlet_refl_1$val[1], rmse_morlet_refl_2$val[1], rmse_morlet_refl_4$val[1], rmse_morlet_refl_8$val[1], rmse_morlet_refl_16$val[1], rmse_morlet_refl_32$val[1],
                         rmse_paul_refl_1$val[1], rmse_paul_refl_2$val[1], rmse_paul_refl_4$val[1], rmse_paul_refl_8$val[1], rmse_paul_refl_16$val[1], rmse_paul_refl_32$val[1],
                         rmse_dog_refl_1$val[1], rmse_dog_refl_2$val[1], rmse_dog_refl_4$val[1], rmse_dog_refl_8$val[1], rmse_dog_refl_16$val[1], rmse_dog_refl_32$val[1],
                         rmse_cauchy_refl_1$val[1], rmse_cauchy_refl_2$val[1], rmse_cauchy_refl_4$val[1], rmse_cauchy_refl_8$val[1], rmse_cauchy_refl_16$val[1], rmse_cauchy_refl_32$val[1],
                         rmse_squiz_refl_1$val[1], rmse_squiz_refl_2$val[1], rmse_squiz_refl_4$val[1], rmse_squiz_refl_8$val[1], rmse_squiz_refl_16$val[1], rmse_squiz_refl_32$val[1],
                         rmse_gabor_refl_1$val[1], rmse_gabor_refl_2$val[1], rmse_gabor_refl_4$val[1], rmse_gabor_refl_8$val[1], rmse_gabor_refl_16$val[1], rmse_gabor_refl_32$val[1])
names(rmse_all_ref) <- names(r2_all_ref) ; abs(rmse_all_ref)

# derivative imput
# r2
r2_all_deriv <- data.frame(r2_morlet_deriv1_1$val[1], r2_morlet_deriv1_2$val[1], r2_morlet_deriv1_4$val[1], r2_morlet_deriv1_8$val[1], r2_morlet_deriv1_16$val[1], r2_morlet_deriv1_32$val[1],
                         r2_paul_deriv1_1$val[1], r2_paul_deriv1_2$val[1], r2_paul_deriv1_4$val[1], r2_paul_deriv1_8$val[1], r2_paul_deriv1_16$val[1], r2_paul_deriv1_32$val[1],
                         r2_dog_deriv1_1$val[1], r2_dog_deriv1_2$val[1], r2_dog_deriv1_4$val[1], r2_dog_deriv1_8$val[1], r2_dog_deriv1_16$val[1], r2_dog_deriv1_32$val[1],
                         r2_cauchy_deriv1_1$val[1], r2_cauchy_deriv1_2$val[1], r2_cauchy_deriv1_4$val[1], r2_cauchy_deriv1_8$val[1], r2_cauchy_deriv1_16$val[1], r2_cauchy_deriv1_32$val[1],
                         r2_squiz_deriv1_1$val[1], r2_squiz_deriv1_2$val[1], r2_squiz_deriv1_4$val[1], r2_squiz_deriv1_8$val[1], r2_squiz_deriv1_16$val[1], r2_squiz_deriv1_32$val[1],
                         r2_gabor_deriv1_1$val[1], r2_gabor_deriv1_2$val[1], r2_gabor_deriv1_4$val[1], r2_gabor_deriv1_8$val[1], r2_gabor_deriv1_16$val[1], r2_gabor_deriv1_32$val[1])
names(r2_all_deriv) <- names(r2_all_ref) ; abs(r2_all_deriv)

# RMSE
rmse_all_deriv <- data.frame(rmse_morlet_deriv1_1$val[1], rmse_morlet_deriv1_2$val[1], rmse_morlet_deriv1_4$val[1], rmse_morlet_deriv1_8$val[1], rmse_morlet_deriv1_16$val[1], rmse_morlet_deriv1_32$val[1],
                           rmse_paul_deriv1_1$val[1], rmse_paul_deriv1_2$val[1], rmse_paul_deriv1_4$val[1], rmse_paul_deriv1_8$val[1], rmse_paul_deriv1_16$val[1], rmse_paul_deriv1_32$val[1],
                           rmse_dog_deriv1_1$val[1], rmse_dog_deriv1_2$val[1], rmse_dog_deriv1_4$val[1], rmse_dog_deriv1_8$val[1], rmse_dog_deriv1_16$val[1], rmse_dog_deriv1_32$val[1],
                           rmse_cauchy_deriv1_1$val[1], rmse_cauchy_deriv1_2$val[1], rmse_cauchy_deriv1_4$val[1], rmse_cauchy_deriv1_8$val[1], rmse_cauchy_deriv1_16$val[1], rmse_cauchy_deriv1_32$val[1],
                           rmse_squiz_deriv1_1$val[1], rmse_squiz_deriv1_2$val[1], rmse_squiz_deriv1_4$val[1], rmse_squiz_deriv1_8$val[1], rmse_squiz_deriv1_16$val[1], rmse_squiz_deriv1_32$val[1],
                           rmse_gabor_deriv1_1$val[1], rmse_gabor_deriv1_2$val[1], rmse_gabor_deriv1_4$val[1], rmse_gabor_deriv1_8$val[1], rmse_gabor_deriv1_16$val[1], rmse_gabor_deriv1_32$val[1])
names(rmse_all_deriv) <- names(r2_all_ref) ; abs(rmse_all_deriv)

save.image("WAVELETS.RData")
}

###############################
### Choose the best dataset ###
##############################
# x = r2_all_ref
# y = rmse_all_ref
# z = r2_all_deriv
# r = rmse_all_deriv

bestWA <- function(x, y, z, r) {
output <- c(which.max(x), max(x), which.min(y), min(y),  which.max(y), max(y),  which.min(r), min(r))
return(output)
}

plotWA <- function(x, y, z, r,...){
## Plot result per family
# x = r2_all_ref
# y = rmse_all_ref
# z = r2_all_deriv
# r = rmse_all_deriv

Scale <- c(1,2,4,8,16,32)
par(mfrow=c(2,2),lend = 1, mai = c(1, 1, 0.4, 0.4))
#r2
# reflive input
plot (Scale, x[,1:6], type="o", ylab=expression(r^2), ylim=c(0,1), pch=1, lwd=1.5, main="", las=1, cex=1.1, cex.lab=1.1, cex.axis=1.1)
lines(Scale, x[,7:12], type="o", pch=2, lwd=1.5, lty=2)
lines(Scale, x[,13:18], type="o", pch=3, lwd=1.5, lty=3)
lines(Scale, x[,19:24], type="o", pch=4, lwd=1.5, lty=4)
lines(Scale, x[,25:30], type="o", pch=5, lwd=1.5, lty=5)
lines(Scale, x[,31:36], type="o", pch=6, lwd=1.5, lty=6)
mtext("A", side=3, line=0.5, adj=0, cex=1.3)
legend("bottomright", legend=c("Morlet", "Paul", "DOG", "Cauchy", "Squeezed", "Gabor"), lty=1:6, pch=1:6)
## derivative input
plot (Scale, z[,1:6], type="o", ylab=expression(r^2), ylim=c(0,1), pch=1, lwd=1.5, main="", las=1, cex=1.1, cex.lab=1.1, cex.axis=1.1)
lines(Scale, z[,7:12], type="o", pch=2, lwd=1.5, lty=2)
lines(Scale, z[,13:18], type="o", pch=3, lwd=1.5, lty=3)
lines(Scale, z[,19:24], type="o", pch=4, lwd=1.5, lty=4)
lines(Scale, z[,25:30], type="o", pch=5, lwd=1.5, lty=5)
lines(Scale, z[,31:36], type="o", pch=6, lwd=1.5, lty=6)
mtext("B", side=3, line=0.5, adj=0, cex=1.3)
# RMSE
# reflive input
plot (Scale, y[,1:6], type="o", ylab="RMSE", ylim=c(min(y), max(y)), pch=1, lwd=1.5, main="", las=1, cex=1.1, cex.lab=1.1, cex.axis=1.1)
lines(Scale, y[,7:12], type="o", pch=2, lwd=1.5, lty=2)
lines(Scale, y[,13:18], type="o", pch=3, lwd=1.5, lty=3)
lines(Scale, y[,19:24], type="o", pch=4, lwd=1.5, lty=4)
lines(Scale, y[,25:30], type="o", pch=5, lwd=1.5, lty=5)
lines(Scale, y[,31:36], type="o", pch=6, lwd=1.5, lty=6)
mtext("C", side=3, line=0.5, adj=0, cex=1.3)
## derivative input
plot (Scale, r[,1:6], type="o", ylab="RMSE", ylim=c(min(r), max(r)), pch=1, lwd=1.5, main="", las=1, cex=1.1, cex.lab=1.1, cex.axis=1.1)
lines(Scale, r[,7:12], type="o", pch=2, lwd=1.5, lty=2)
lines(Scale, r[,13:18], type="o", pch=3, lwd=1.5, lty=3)
lines(Scale, r[,19:24], type="o", pch=4, lwd=1.5, lty=4)
lines(Scale, r[,25:30], type="o", pch=5, lwd=1.5, lty=5)
lines(Scale, r[,31:36], type="o", pch=6, lwd=1.5, lty=6)
mtext("D", side=3, line=0.5, adj=0, cex=1.3)
}
