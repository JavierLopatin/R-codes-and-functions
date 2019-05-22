
################################################################################
#
# Plot the relative sensitivity of leaf and canopy plant traits in the
# optical electromagnetic spectrum (400 - 2500 nm) according to PROSAIL
# radiative transfer model.
#
#
################################################################################

require(hsdar)
setwd("D:/googledrive/fieldcampaign_2016_KA_botgarden/1_analysis/")

# load traits values to be used as parameter ranges
traits <- read.csv("time_series/trait_summary_global.csv", sep="")

no_sim = 10000

# lai
parameter <- data.frame(N = c(rep(1.5, no_sim)),
                        LAI = rnorm(no_sim, mean = traits$lai[1], traits$lai[3]),
                        Cab = c(rep(traits$cab[1], no_sim)),
                        Car = c(rep(traits$car[1], no_sim)),
                        Cm = c(rep(traits$cm[1], no_sim)),
                        Cw = c(rep(traits$cw[1], no_sim)))
spec_lai = spectra(PROSAIL(parameterList = parameter))
sens_lai = apply(as.matrix(spec_lai), 2, FUN = mad)

plot(400:2500, sens_lai, ylim=c(0,1), type="l")
for(i in 1:no_sim){
  lines(400:2500,spec_lai[i,])
}
lines(400:2500, sens_lai, col="red")

# Cab
parameter <- data.frame(N = c(rep(1.5, no_sim)),
                        LAI = c(rep(traits$lai[1], no_sim)),
                        Cab = rnorm(no_sim, mean = traits$cab[1], traits$cab[3]),
                        Car = c(rep(traits$car[1], no_sim)),
                        Cm = c(rep(traits$cm[1], no_sim)),
                        Cw = c(rep(traits$cw[1], no_sim)))
spec_cab = spectra(PROSAIL(parameterList = parameter))
sens_cab = apply(as.matrix(spec_cab), 2, FUN = mad)

plot(400:2500, sens_cab, ylim=c(0,1), type="l")
for(i in 1:no_sim){
  lines(400:2500,spec_cab[i,])
}
lines(400:2500, sens_cab, col="red")

# Car
parameter <- data.frame(N = c(rep(1.5, no_sim)),
                        LAI = c(rep(traits$lai[1], no_sim)),
                        Cab = c(rep(traits$cab[1], no_sim)),
                        Car = rnorm(no_sim, mean = traits$car[1], traits$car[3]),
                        Cm = c(rep(traits$cm[1], no_sim)),
                        Cw = c(rep(traits$cw[1], no_sim)))
spec_car = spectra(PROSAIL(parameterList = parameter))
sens_car = apply(as.matrix(spec_car), 2, FUN = mad)

plot(400:2500, sens_car, ylim=c(0,1), type="l")
for(i in 1:no_sim){
  lines(400:2500,spec_car[i,])
}
lines(400:2500, sens_car, col="red")

# Cm
parameter <- data.frame(N = c(rep(1.5, no_sim)),
                        LAI = c(rep(traits$lai[1], no_sim)),
                        Cab = c(rep(traits$cab[1], no_sim)),
                        Car = c(rep(traits$car[1], no_sim)),
                        Cm = rnorm(no_sim, mean = traits$cm[1], traits$cm[3]),
                        Cw = c(rep(traits$cw[1], no_sim)))
spec_cm = spectra(PROSAIL(parameterList = parameter))
sens_cm = apply(as.matrix(spec_cm), 2, FUN = mad)

plot(400:2500, sens_cm, ylim=c(0,1), type="l")
for(i in 1:no_sim){
  lines(400:2500,spec_cm[i,])
}
lines(400:2500, sens_cm, col="red")

# Cw
parameter <- data.frame(N = c(rep(1.5, no_sim)),
                        LAI = c(rep(traits$lai[1], no_sim)),
                        Cab = c(rep(traits$cab[1], no_sim)),
                        Car = c(rep(traits$car[1], no_sim)),
                        Cm = c(rep(traits$cm[1], no_sim)),
                        Cw = rnorm(no_sim, mean = traits$cw[1], traits$cw[3]))
spec_cw = spectra(PROSAIL(parameterList = parameter))
sens_cw = apply(as.matrix(spec_cw), 2, FUN = mad)

plot(400:2500, sens_cw, ylim=c(0,1), type="l")
for(i in 1:no_sim){
  lines(400:2500,spec_cw[i,])
}
lines(400:2500, sens_cw, col="red")


# Leaf angles
parameter <- data.frame(N = c(rep(1.5, no_sim)),
                        LAI = c(rep(traits$lai[1], no_sim)),
                        Cab = c(rep(traits$cab[1], no_sim)),
                        Car = c(rep(traits$car[1], no_sim)),
                        Cm = c(rep(traits$cm[1], no_sim)),
                        Cw = c(rep(traits$cw[1], no_sim)),
                        lidfa = rnorm(no_sim, mean = traits$lidfa[1], traits$lidfa[3]),
                        lidfb = rnorm(no_sim, mean = traits$lidfb[1], traits$lidfb[3]))
spec_ang = spectra(PROSAIL(parameterList = parameter))
sens_ang = apply(as.matrix(spec_ang), 2, FUN = mad)

plot(400:2500, sens_ang, ylim=c(0,1), type="l")
for(i in 1:no_sim){
  lines(400:2500,spec_ang[i,])
}
lines(400:2500, sens_ang, col="red")



plot(400:2500, ylim=c(0,0.15))
lines(sens_lai, col="black")
lines(sens_cab, col="orange")
lines(sens_car, col="red")
lines(sens_cm, col="green")
lines(sens_cw, col="blue")
lines(sens_ang, col="purple")
legend("topright", c("lai", "cab", "car", "cm", "cw", "angle"), col=c("black", "orange", "red", "green", "blue", "purple"), lwd=1)


sens_all = rbind(sens_lai,
                 sens_cab,
                 sens_car,
                 sens_cm,
                 sens_cw,
                 sens_ang)


write.table(sens_all, file="simulate_PFT/prosail_sensivity_traits.csv")
