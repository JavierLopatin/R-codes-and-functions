## R-Script - Selecting hyperspectral bands for classification using Genetic Algorithm
## author: Javier Lopatin & Fabian Fassnacht
## mails: javierlopatin@gmail.com & fabian.fassnacht@kit.edu
## last changes: 12-07-2015
## info: - this script use the "galgo" package. This package require the R version 2.12.1
##       - in the input data file should be ordered rows (not in columns). The first row should be "id", then "class" (with all the classes that you have) and finally the hyperspectral bands "B1",..., "Bn"

# load galgo package
library(galgo)

# set working directory 
setwd("your/diretory")

# select your data file
data<-"galgo_20m_all.txt"

#preparar los especificaciones del algoritmo
Genetic.Algorith <-configBB.VarSel(file=data, classification.method="nearcent", chromosomeSize=4, 
                   maxSolutions=100, goalFitness= 0.99, saveVariable = "bb.5chr_1000sol_095goal", saveFrequency = 30)

# run galgo
blast(Genetic.Algorith)

# see confusion matrix
plot(Genetic.Algorith, type="confusion", set=c(1,0), splits=1)

# see stability of the predictors
plot(Genetic.Algorith, type="generankstability", xlim=c(0,41), xlab="Bands rank", ylab="Frequency chosen", main="")

# see models in a "stepwise" fashion
fsm <- forwardSelectionModels(Genetic.Algorith)
fsm$models


