###################################################
## Obtain the first derivative of the spectra 
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
