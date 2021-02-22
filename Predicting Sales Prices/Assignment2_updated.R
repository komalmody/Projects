graphics.off() # This closes all of R's graphics windows.
rm(list = ls())  # Careful! This clears all of R's memory!

# installing packages and calling libraries

#install.packages("dplyr")

library("dplyr") # For Data pre-processing manipulation

library(ggplot2)

library(ggpubr)

library(ks) 

library(rjags) # For initiating the Jags model

library(runjags) # For executing the Jags model

#--------------------------------------------------

# setting the working directory

setwd(
  "C:/Users/komal/OneDrive/Documents/Update RMIT Study Material/3rd Semester/Applied Baysein Statistics/Assignment2"
)

#--------------------------------------------------

# Importing the utilities function

source("DBDA2E-utilities.R")

#===============PRELIMINARY FUNCTIONS FOR POSTERIOR INFERENCES====================

smryMCMC_Sale_Price = function(codaSamples ,
                               compVal = NULL,
                               saveName = NULL) {
  summaryInfo = NULL
  mcmcMat = as.matrix(codaSamples, chains = TRUE)
  paramName = colnames(mcmcMat)
  for (pName in paramName) {
    if (pName %in% colnames(compVal)) {
      if (!is.na(compVal[pName])) {
        summaryInfo = rbind(summaryInfo ,
                            summarizePost(
                              paramSampleVec = mcmcMat[, pName] ,
                              compVal = as.numeric(compVal[pName])
                            ))
      }
      else {
        summaryInfo = rbind(summaryInfo , summarizePost(paramSampleVec = mcmcMat[, pName]))
      }
    } else {
      summaryInfo = rbind(summaryInfo , summarizePost(paramSampleVec = mcmcMat[, pName]))
    }
  }
  rownames(summaryInfo) = paramName
  
  
  if (!is.null(saveName)) {
    write.csv(summaryInfo , file = paste(saveName, "SummaryInfo.csv", sep =
                                           ""))
  }
  return(summaryInfo)
}

  #=========================Defining plot method for construction of posterior distribution======================================================
  plotMCMC_Sale_Price = function( codaSamples , data, xName="x" , yName="y" ,
                        showCurve=FALSE ,  pairsPlot=FALSE , compVal = NULL,
                        saveName=NULL , saveType="jpg" ) {
  #-----------------------------------------------------------------------------
  
  #====================Separating Independent and dependent variables=====================
  
  y = data[, yName]
  x = as.matrix(data[, xName])
  
  # ======================Initializing the priors and predicion variables=================
  
  mcmcMat = as.matrix(codaSamples, chains = TRUE)
  chainLength = NROW(mcmcMat)
  zbeta0 = mcmcMat[, "zbeta0"]
  zbeta  = mcmcMat[, grep("^zbeta$|^zbeta\\[", colnames(mcmcMat))]
  if (ncol(x) == 1) {
    zbeta = matrix(zbeta , ncol = 1)
  }
  zVar = mcmcMat[, "zVar"]
  beta0 = mcmcMat[, "beta0"]
  beta  = mcmcMat[, grep("^beta$|^beta\\[", colnames(mcmcMat))]
  if (ncol(x) == 1) {
    beta = matrix(beta , ncol = 1)
  }
  tau = mcmcMat[, "tau"]
  pred1 = mcmcMat[, "pred[1]"] # Prediction 1 for 1st Settings
  pred2 = mcmcMat[, "pred[2]"] # Prediction 2 for 2nd settings
  pred3 = mcmcMat[, "pred[3]"] # Prediction 3 for 3rd settings
  pred4 = mcmcMat[, "pred[4]"] # Prediction 4 for 4th settings
  pred5 = mcmcMat[, "pred[5]"] # Prediction 5 for 5th settings
  
  #-----------------------------------------------------------------------------
  
  # Compute R^2 for credible parameters:
  
  YcorX = cor(y , x) # correlation of y with each x predictor
  
  Rsq = zbeta %*% matrix(YcorX , ncol = 1)
  
  
  #-----------------------------------------------------------------------------
  
  # Marginal histograms:
  
  decideopenGraph_Sale_Price = function( panelCount , saveName , finished=FALSE , 
                              nRow=2 , nCol=3 ) {
    # If finishing a set:
    if ( finished==TRUE ) {
      if ( !is.null(saveName) ) {
        saveGraph( file=paste0(saveName,ceiling((panelCount-1)/(nRow*nCol))), 
                   type=saveType)
      }
      panelCount = 1 # re-set panelCount
      return(panelCount)
    } else {
      # If this is first panel of a graph:
      if ( ( panelCount %% (nRow*nCol) ) == 1 ) {
        # If previous graph was open, save previous one:
        if ( panelCount>1 & !is.null(saveName) ) {
          saveGraph( file=paste0(saveName,(panelCount%/%(nRow*nCol))), 
                     type=saveType)
        }
        # Open new graph
        openGraph(width=nCol*7.0/3,height=nRow*2.0)
        layout( matrix( 1:(nRow*nCol) , nrow=nRow, byrow=TRUE ) )
        par( mar=c(4,4,2.5,0.5) , mgp=c(2.5,0.7,0) )
      }
      # Increment and return panel count:
      panelCount = panelCount+1
      return(panelCount)
    }
  }
  
  # Original scale:
  panelCount = 1
  if (!is.na(compVal["beta0"])){
    panelCount = decideopenGraph_Sale_Price( panelCount , saveName=paste0(saveName,"PostMarg") )
    histInfo = plotPost( beta0 , cex.lab = 1.75 , showCurve=showCurve ,
                         xlab=bquote(beta[0]) , main="Intercept", compVal = as.numeric(compVal["beta0"] ))
  } else {  
    histInfo = plotPost( beta0 , cex.lab = 1.75 , showCurve=showCurve ,
                         xlab=bquote(beta[0]) , main="Intercept")
  }
  for ( bIdx in 1:ncol(beta) ) {
    panelCount = decideopenGraph_Sale_Price( panelCount , saveName=paste0(saveName,"PostMarg") )
    if (!is.na(compVal[paste0("beta[",bIdx,"]")])) {
      histInfo = plotPost( beta[,bIdx] , cex.lab = 1.75 , showCurve=showCurve ,
                           xlab=bquote(beta[.(bIdx)]) , main=xName[bIdx],
                           compVal = as.numeric(compVal[paste0("beta[",bIdx,"]")]))
    } else{
      histInfo = plotPost( beta[,bIdx] , cex.lab = 1.75 , showCurve=showCurve ,
                           xlab=bquote(beta[.(bIdx)]) , main=xName[bIdx])
    }
  }
  panelCount = decideopenGraph_Sale_Price( panelCount , saveName=paste0(saveName,"PostMarg") )
  histInfo = plotPost( tau , cex.lab = 1.75 , showCurve=showCurve ,
                       xlab=bquote(tau) , main=paste("Scale") )
  panelCount = decideopenGraph_Sale_Price( panelCount , saveName=paste0(saveName,"PostMarg") )
  histInfo = plotPost( Rsq , cex.lab = 1.75 , showCurve=showCurve ,
                       xlab=bquote(R^2) , main=paste("Prop Var Accntd") )
  panelCount = decideopenGraph_Sale_Price( panelCount ,  saveName=paste0(saveName,"PostMarg") )
  histInfo = plotPost( pred1 , cex.lab = 1.75 , showCurve=showCurve ,
                       xlab="pred1" , main="Prediction 1" ) # Added by Demirhan
  panelCount = decideopenGraph_Sale_Price( panelCount , saveName=paste0(saveName,"PostMarg") )
  histInfo = plotPost( pred2 , cex.lab = 1.75 , showCurve=showCurve ,
                       xlab="pred2" , main="Prediction 2" ) # Added by Demirhan
  panelCount = decideopenGraph_Sale_Price( panelCount ,  saveName=paste0(saveName,"PostMarg") )
  histInfo = plotPost( pred3 , cex.lab = 1.75 , showCurve=showCurve ,
                       xlab="pred3" , main="Prediction 3" ) # Added by Demirhan
  panelCount = decideopenGraph_Sale_Price( panelCount ,  saveName=paste0(saveName,"PostMarg") )
  histInfo = plotPost( pred4 , cex.lab = 1.75 , showCurve=showCurve ,
                       xlab="pred4" , main="Prediction 4" ) # Added by Demirhan
  panelCount = decideopenGraph_Sale_Price( panelCount , finished=TRUE , saveName=paste0(saveName,"PostMarg") )
  histInfo = plotPost( pred5 , cex.lab = 1.75 , showCurve=showCurve ,
                       xlab="pred5" , main="Prediction 5" ) # Added by Demirhan
  
  # Standardized scale:
  panelCount = 1
  panelCount = decideopenGraph_Sale_Price( panelCount , saveName=paste0(saveName,"PostMargZ") )
  histInfo = plotPost( zbeta0 , cex.lab = 1.75 , showCurve=showCurve ,
                       xlab=bquote(z*beta[0]) , main="Intercept" )
  for ( bIdx in 1:ncol(beta) ) {
    panelCount = decideopenGraph_Sale_Price( panelCount , saveName=paste0(saveName,"PostMargZ") )
    histInfo = plotPost( zbeta[,bIdx] , cex.lab = 1.75 , showCurve=showCurve ,
                         xlab=bquote(z*beta[.(bIdx)]) , main=xName[bIdx] )
  }
  panelCount = decideopenGraph_Sale_Price( panelCount , saveName=paste0(saveName,"PostMargZ") )
  histInfo = plotPost( zVar , cex.lab = 1.75 , showCurve=showCurve ,
                       xlab=bquote(z*tau) , main=paste("Scale") )
  panelCount = decideopenGraph_Sale_Price( panelCount , saveName=paste0(saveName,"PostMargZ") )
  histInfo = plotPost( Rsq , cex.lab = 1.75 , showCurve=showCurve ,
                       xlab=bquote(R^2) , main=paste("Prop Var Accntd") )
  panelCount = decideopenGraph_Sale_Price( panelCount , finished=TRUE , saveName=paste0(saveName,"PostMargZ") )
  
  #-----------------------------------------------------------------------------
}


#===============PRELIMINARY FUNCTIONS FOR POSTERIOR INFERENCES====================

#-------------- Importing the Assignment2PropertySalePrices in Data_Sale_Price ---

Data_Sale_Price <-
  read.csv(
    "C:/Users/komal/OneDrive/Documents/Update RMIT Study Material/3rd Semester/Applied Baysein Statistics/Assignment2/Assignment2PropertyPrices.csv"
  )
head(Data_Sale_Price)


# Setting the seed sample
set.seed(130500)

Sale_Price_Sample = sample_n(Data_Sale_Price, 1500)
Sale_Price_Sample

# Scatter plots to examine correlation between
# dependent and Independent variables

p1 <- ggplot(Sale_Price_Sample, aes(x = Area, y = SalePrice)) +
  geom_point() +
  xlab("Area") +
  ylab("SalePrice")

p2 <-
  ggplot(
    Sale_Price_Sample,
    aes(x = Bedrooms, y = SalePrice)
  ) +
  geom_point() +
  xlab("Bedrooms") +
  ylab("SalePrice")


p3 <- ggplot(Sale_Price_Sample, aes(x = Bathrooms, y = SalePrice)) +
  geom_point() +
  xlab("Bathrooms") +
  ylab("SalePrice")


p4 <- ggplot(Sale_Price_Sample, aes(x = CarParks, y = SalePrice)) +
  geom_point() +
  xlab("CarParks") +
  ylab("SalePrice")

p5 <-
  ggplot(Sale_Price_Sample, aes(x = PropertyType, y = SalePrice)) +
  geom_point() +
  xlab("PropertyType") +
  ylab("SalePrice")

figure <- ggarrange(p1, p2, p3, p4, p5, nrow = 2, ncol = 3)
figure

# Histogram of the dependent variable
hist(Data_Sale_Price$SalePrice, main = " Histogram of the Sale Price(Dependent Variable)", xlab = "SalePrice")

# Kernel density estimation of Dependent Variable SalePrice
plot(kde(Data_Sale_Price$SalePrice), xlab = "SalePrice") # with default settings


# Separating the Dependent and Independent variable from the sample data
y = Data_Sale_Price[, "SalePrice"]
x = as.matrix(Data_Sale_Price[, c("Area", "Bedrooms", "Bathrooms", "CarParks", "PropertyType")])

# Some more descriptives
cat("\nCORRELATION MATRIX OF PREDICTORS:\n ")
show(round(cor(x), 3))
cat("\n")

xPred = array(NA, dim = c(5, 5))
xPred[1,] = c(600, 2, 2, 1, 1)
xPred[2,] = c(800, 3, 1, 2, 0)
xPred[3,] = c(1500, 2, 1, 1, 0)
xPred[4,] = c(2500, 5, 4, 4, 0)
xPred[5,] = c(250, 3, 2, 1, 1)

# Recollect Independent and Dependent variables in list

dataList <- list(
  x = x ,
  y = y ,
  xPred = xPred ,
  Nx = dim(x)[2] ,
  Ntotal = dim(x)[1]
)

# First run without initials!
initsList <- list(
  zbeta0 = 200,
  zbeta = c(100, 1, 1, 1, 1),
  Var = 1200
)

# Execution of Jags Model

# THE preliminary for model.

modelString = "
  
  # Standardize the data:
  
  data {
    ysd <- sd(y)
    for ( i in 1:Ntotal ) {
      zy[i] <- y[i] / ysd
    }
    for ( j in 1:Nx ) {
      xsd[j] <-   sd(x[,j])
      for ( i in 1:Ntotal ) {
        zx[i,j] <- x[i,j] / xsd[j]
      }
    }
  }

  # Model Specification for scaled data:

model {
for (i in 1:Ntotal) {
zy[i] ~ dgamma((mu[i] ^ 2) / zVar , mu[i] / zVar)
    mu[i] <- zbeta0 + sum(zbeta[1:Nx] * zx[i, 1:Nx])
  }
  # Priors on standardized scale:
  zbeta0 ~ dnorm(0 , 1 / 2 ^ 2)
  zbeta[1] ~ dnorm(0.0009 / xsd[1] , 1 / (0.2 / xsd[1] ^ 2))
  zbeta[2] ~ dnorm(1 / xsd[2] , 1 / (30 / xsd[2] ^ 2))
  zbeta[3] ~ dnorm(0 , 1 / 2 ^ 2)
  zbeta[4] ~ dnorm(1.2 / xsd[4] , 1 / (2 / xsd[4] ^ 2))
  zbeta[5] ~ dnorm(-1.5 / xsd[5] , 1 / (0.2 / xsd[5] ^ 2))

  zVar ~ dgamma(0.01 , 0.01)

  # Transform to original scale:

  beta[1:Nx] <- (zbeta[1:Nx] / xsd[1:Nx]) * ysd
  beta0 <- zbeta0 * ysd
  tau <- zVar * (ysd) ^ 2

  # Compute predictions at every step of the MCMC

  for (i in 1:5) {
pred[i] <-
 beta0 + beta[1] * xPred[i, 1] + beta[2] * xPred[i, 2] + beta[3] * xPred[i, 3] + beta[4] * xPred[i, 4] + beta[5] * xPred[i, 5]
  }

}
" # close quote for modelString

# Write out modelString to a text file

writeLines(modelString , con = "SalePriceModel.txt")

# Setting for MCMC

adaptSteps = 1800  # Number of steps to "tune" the samplers
burnInSteps = 6000 # Number of Burn-in steps
nChains = 3
thinSteps = 48 # First run for 3
numSavedSteps = 1200

# Parallel run of the jags model

STARTTIME <- Sys.time()
runJagsOut <- run.jags(
  method = "parallel" ,
  model = "SalePriceModel.txt" ,
  monitor = c("zbeta0" ,  "zbeta" , "beta0" ,  "beta" ,  "tau", "zVar", "pred")  ,
  data = dataList  ,
  inits = initsList,
  n.chains = nChains ,
  adapt = adaptSteps ,
  burnin = burnInSteps ,
  sample = numSavedSteps ,
  thin = thinSteps ,
  summarise = FALSE ,
  plots = FALSE
)
codaSamples = as.mcmc.list(runJagsOut)

ENDTIME <- Sys.time()
print(ENDTIME - STARTTIME)
save.image(file = "rEnvironment.RData")
#load(file = "rEnvironment.RData") # Load the results with 27600 iterations

# plotting the MCMC diagnostics plots

diagMCMC(codaSamples , parName = "beta0")
diagMCMC(codaSamples , parName = "beta[1]")
diagMCMC(codaSamples , parName = "beta[2]")
diagMCMC(codaSamples , parName = "beta[3]")
diagMCMC(codaSamples , parName = "beta[4]")
diagMCMC(codaSamples , parName = "beta[5]")
diagMCMC(codaSamples , parName = "tau")
diagMCMC(codaSamples , parName = "pred[1]")
diagMCMC(codaSamples , parName = "pred[2]")
diagMCMC(codaSamples , parName = "pred[3]")
diagMCMC(codaSamples , parName = "pred[4]")
diagMCMC(codaSamples , parName = "pred[5]")
#diagMCMC(codaSamples , parName = "zbeta0")
#diagMCMC(codaSamples , parName = "zbeta[1]")

compVal <-
  data.frame(
    "beta0" = NA,
    "beta[1]" = NA,
    "beta[2]" = NA,
    "beta[3]" = NA,
    "beta[4]" =  NA,
    "beta[5]" =  NA,
    "tau" = NA ,
    check.names = FALSE
  )

# Printing the summary of MCMC diagnostics

summaryInfo_Sale_Price <-
  smryMCMC_Sale_Price(codaSamples = codaSamples , compVal = compVal)

print(summaryInfo_Sale_Price)



# plotting posterior distribution
plotMCMC_Sale_Price(
  codaSamples = codaSamples ,
  data = Data_Sale_Price,
  xName = c("Area", "Bedrooms", "Bathrooms", "CarParks", "PropertyType") ,
  yName = "SalePrice",
  compVal = compVal
)

# ============ Predictive check ============

coefficients <-  summaryInfo_Sale_Price[8:13, 3] # Get the model coefficients out
Variance <- summaryInfo_Sale_Price[14, 3] # Get the variance out
# Since we imposed the regression model on the mean of the gamma likelihood,
# we use the model (X*beta) to generate the mean of gamma population for each
# observed x vector.

meanGamma <- as.matrix(cbind(rep(1, nrow(x)),  x)) %*% as.vector(coefficients)

# Generate random data from the posterior distribution. Here I take the
# reparameterisation back to alpha and beta.

randomData <-
  rgamma(n = 231,
         shape = meanGamma ^ 2 / Variance,
         rate = meanGamma / Variance)


# Display the density plot of observed data and posterior distribution:

predicted <- data.frame(Sale_Price = randomData)
observed <- data.frame(Sale_Price = y)
predicted$type <- "Predicted"
observed$type <- "Observed"
dataPred <- rbind(predicted, observed)

ggplot(dataPred, aes(elapsed, fill = type)) + 
  geom_density(alpha = 0.2)+ 
  ggtitle("Density Distribution of observed and predicted posterior distribution")




failed.jags('model')
failed.jags()
cleanup.jags()


