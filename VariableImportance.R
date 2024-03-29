########################################################################
# Example getting variabe importance using MRPP y Delta of Kappa values

# Author: Javier Lopatin | javier.lopatin@uai.cl

# Justification of the approach:

# To assess the relative contribution of each indicator to classify certain classes, we used the Multi Response Permutation Procedure (MRPP; Mielke, 1991; McCune and Grace, 2002). 
# The MRPP is a multivariate non-parametric test of whether there is a significant difference between groups within each individual variable. 
# The MRPP provides change-corrected group agreement (A) and significance (P) values. Similar to the coefficient of determination, A ranges from 0 to 1, showing the level of discrimination between groups. 
# Accordingly, a hypothetical A value of 1 implies that an indicator thoroughly explains the variance between degradation levels. 
# In contrast, an A value of 0 implies that the indicator does not explain the degradation levels. 
# MRPP has shown robust results applied to biological and radiometric analyses, however, the potential synergies among multiple indicators are not considered in the indicator-wise MRPP-based analysis. 
# Therefore, we considered these synergies among variables by applying a partial least square (PLS) discriminant analysis. 
# We determined the contribution of each indicator to discern the level of degradation by applying a bootstrapping iteration procedure, wherein each iteration (100) we:

# 1. Fitted a general model using all indicators and observations available and stored the overall Kappa (Kall) value;
# 2. Fitted a indicator-wise partial models by randomizing the values of one indicator at the time in a stepwise procedure, storing one K value for each indicator replacement (Ki);
# 3. Estimated the relative contribution of each indicator by subtracting the indicator-wise partial model from the overall model (Kall ₋ Ki), generating a delta Kappa per indicator (ΔK).

# We used on average 63% samples for model training and 37% for validation during each iteration. We used PLS models on the training samples using a 5-fold cross-validation approach. We stored the results of the 100 iterations to present the distribution of ΔK and prevent stochastic biases. 
# We used the 'caret' and 'vegan' packages of the R statistical software (R Core Team, 2014).

########################################################################


setwd('~dir')

data <- read.csv('data.csv', sep = ',', dec=',', stringsAsFactors = TRUE)
str(data)

# ----------------------------------------------
# MRPP - Multi Response Permutation Procedure
# ----------------------------------------------

library(vegan)

var <- data[, 4:39]
str(var)

mat_mrpp = matrix(NA, ncol=ncol(var), nrow=2)
colnames(mat_mrpp) <- colnames(var)
rownames(mat_mrpp) <- c('A', 'P')
for(i in 1:ncol(var)){
  obj_mrpp = mrpp(dat = var[,i], grouping = data$degra, parallel = 8, distance = "mahalanobis")
  mat_mrpp[1,i] = obj_mrpp$A
  mat_mrpp[2,i] = obj_mrpp$Pvalue
}

mat_mrpp

# significance
which(mat_mrpp[2,] < 0.05)

svg(filename = 'MRPP.svg', width = 10, height = 5)
barplot(mat_mrpp[1,], las=2, ylab = 'A')
box()
dev.off()


# ----------------------------------------------
# PLS-DA - Partial Least Square Discriminant Analysis
# ----------------------------------------------

library(caret)

# Asegurar que la columna "clases" sea de tipo factor
data$degra <- as.factor(data$degra)
str(data)

var$class <- data$degra 
str(var)

# variables
X=data[,3:ncol(dta)]
Y=data$vegclass

# get best ncomp values with caret
tr_control <- trainControl(method = 'cv', number = 10)

# train model using vegetation clases
PLS <- train(x=X, y=Y , method = "pls",  tuneLength=20, 
             trControl = tr_control, preProcess = c("center", "scale"), 
             metric = "Kappa", maximize = T, na.action=na.omit)
PLS

# function
relPLSimp = function(X, Y, ncomp, iter=100){
    #############################################################
    #
    # Estimate relative importance using PLS-DA by using Kappa 
    # and predicting values in random variables
    #
    # Y = observed classes (vector of Factor class)
    # X = predictors (data frame of numerical variables)
    # iter = number of bootstrap iterations (intener)
    #
    ##############################################################
     
    kappa <- matrix(nrow = iter, ncol = ncol(X))
    colnames(kappa) = colnames(X)
    
    boot <- createResample(Y, times = iter, list = TRUE)
    
    pb <- txtProgressBar(min = 0, max = iter, style = 3, width = 50, char = "=") 

   for (i in 1:iter){ # bootstraps
          
      train_X <- X[boot[[i]],]
      train_Y <- Y[boot[[i]]]  
      val_X   <- X[-boot[[i]],]
      val_Y   <- Y[-boot[[i]]]

      PLS  <- caret::plsda(x = train_X, y = train_Y, ncomp = ncomp, probMethod = 'softmax')
      pred <- predict(PLS, val_X)

      mat_pls <- confusionMatrix(pred, val_Y)

      for (k in 1:(ncol(X))){ # loop through variables
        validar_rand  <- val_X
        kappp = c()
         for (j in 1:10){ # 10 random replaces per variable
             validar_rand[,k]  <- runif(nrow(validar_rand), min=min(validar_rand[,k]), max=max(validar_rand[,k]))
             pred_rand <- predict(PLS, validar_rand)
             mat_pls_rand <- confusionMatrix(pred_rand, val_Y)
             kappp[j] <- mat_pls_rand$overall[2]
          }
        # delta kappa
        kappa[i,k] <- mat_pls$overall[2] - median(kappp)
      }
    setTxtProgressBar(pb, i)
    }
    close(pb)
    return(kappa)
}


imp = relPLSimp(X=X, Y=Y, ncomp=10, iter = 100)
imp

save.image('varImp.RData')

svg(filename = 'kappa.svg', width = 10, height = 5)
boxplot((kappa), las=2, outline=F)
abline(0,0, lty=2)
dev.off()

# Significance alpha = 0.05
for (i in 1:36) print(quantile(kappa[,i]))
