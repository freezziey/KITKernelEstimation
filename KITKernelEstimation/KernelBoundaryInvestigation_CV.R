# _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _
#
# Boundary Kernel Density Estimation
# Investigation of KDE behavior near the boundary
# 15.06.16
#
# _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _


#Include the neccessary packages and sourcefiles
install.packages("Rcpp")
library(Rcpp)
sourceCpp(
  "/Users/Christian/Documents/Christian/XCode workspace/KITKernelEstimation/KITKernelEstimation/CPPKernelEstimation.cpp"
)

install.packages("xtable")

library(xtable)

#Implementation of neccessary functions
#______________________________________________________

#MSE
MSE <- function(real, predicted) {
  result <- ((predicted - real) ^ 2)
  return(result)
}

#MISE
MISE <- function(real, predicted) {
  result <- ((predicted - real) ^ 2)
  result <- mean(result)
  return(result)
}


arg_min <- function(x, y) {
  argmin <- y[1]
  min <- x[1]
  for (i in 1:length(x)) {
    if (y[i] < argmin) {
      argmin <- y[i]
      min <- x[i]
    }
  }
  return(min)
}


#Generate samples doubled Norm
generateSamples_doubledNorm <- function(numberOfSamples, sd) {
  #seed<-6795999
  result <- list()
  for (i in 1:numberOfSamples) {
    #set.seed(seed)
    result[[i]] <- abs(rnorm(sizeOfSimulatedVector, sd = sd))
    #seed<-seed+5
  }
  return(result)
}

#Generate samples gamma
generateSamples_gamma <- function(numberOfSamples, p, b) {
  #seed<-7685865
  result <- list()
  for (i in 1:numberOfSamples) {
    # set.seed(seed)
    result[[i]] <-
      rgamma(sizeOfSimulatedVector, shape = p, rate = b)
    #seed<-seed+5
  }
  return(result)
}

#Generate samples sinus
generateSamples_sinus <- function(numberOfSamples) {
  #seed<-7685865
  result <- list()
  for (i in 1:numberOfSamples) {
    # set.seed(seed)
    result[[i]] <- generateSinusDistribution()
    #seed<-seed+5
  }
  return(result)
}


#Bias gamma_MOD
Bias_gamma_MOD <- function(real, samples) {
  result <- numeric()
  temp <- list()
  for (s in 1:length(samples)) {
    temp[[s]] <-
      CPP_KDE_gamma_MOD(as.numeric(real$x), as.numeric(unlist(samples[s])), h_gamma_MOD)
  }
  mean <- numeric(length(real$x))
  for (i in 1:length(temp)) {
    mean <- mean + as.numeric(unlist(temp[i]))
  }
  mean <- mean / length(samples)
  result <- mean - real$y
  return(list(x = real$x, y = result))
}

#Bias gamma_EXP
Bias_gamma_EXP <- function(real, samples) {
  result <- numeric()
  temp <- list()
  for (s in 1:length(samples)) {
    temp[[s]] <-
      CPP_KDE_gamma_EXP(as.numeric(real$x), as.numeric(unlist(samples[s])), h_gamma_EXP)
  }
  mean <- numeric(length(real$x))
  for (i in 1:length(temp)) {
    mean <- mean + as.numeric(unlist(temp[i]))
  }
  mean <- mean / length(samples)
  result <- mean - real$y
  return(list(x = real$x, y = result))
}

#Bias gammaChen_MOD
Bias_gammaChen_MOD <- function(real, samples) {
  result <- numeric()
  temp <- list()
  for (s in 1:length(samples)) {
    temp[[s]] <-
      CPP_KDE_gammaChen_MOD(as.numeric(real$x), as.numeric(unlist(samples[s])), h_gammaChen_MOD)
  }
  mean <- numeric(length(real$x))
  for (i in 1:length(temp)) {
    mean <- mean + as.numeric(unlist(temp[i]))
  }
  mean <- mean / length(samples)
  result <- mean - real$y
  return(list(x = real$x, y = result))
}

#Bias gammaChen_EXP
Bias_gammaChen_EXP <- function(real, samples) {
  result <- numeric()
  temp <- list()
  for (s in 1:length(samples)) {
    temp[[s]] <-
      CPP_KDE_gammaChen_EXP(as.numeric(real$x), as.numeric(unlist(samples[s])), h_gammaChen_EXP)
  }
  mean <- numeric(length(real$x))
  for (i in 1:length(temp)) {
    mean <- mean + as.numeric(unlist(temp[i]))
  }
  mean <- mean / length(samples)
  result <- mean - real$y
  return(list(x = real$x, y = result))
}

#Bias rectangular
Bias_rectangular <- function(real, samples) {
  result <- numeric()
  temp <- list()
  for (s in 1:length(samples)) {
    temp[[s]] <-
      CPP_KDE_rectanuglar(as.numeric(real$x),
                          as.numeric(unlist(samples[s])),
                          h_rectangular,
                          lowerBoundary = 0)
  }
  mean <- numeric(length(real$x))
  for (i in 1:length(temp)) {
    mean <- mean + as.numeric(unlist(temp[i]))
  }
  mean <- mean / length(samples)
  result <- mean - real$y
  return(list(x = real$x, y = result))
}

#Bias gamma_MOD_EXP
Bias_gamma_MOD_EXP <- function(real, samples) {
  result <- numeric()
  temp <- list()
  for (s in 1:length(samples)) {
    temp[[s]] <-
      0.5 * (
        CPP_KDE_gamma_MOD(as.numeric(real$x), as.numeric(unlist(samples[s])), h_gamma_MOD) +
          CPP_KDE_gamma_EXP(as.numeric(real$x), as.numeric(unlist(samples[s])), h_gamma_EXP)
      )
  }
  mean <- numeric(length(real$x))
  for (i in 1:length(temp)) {
    mean <- mean + as.numeric(unlist(temp[i]))
  }
  mean <- mean / length(samples)
  result <- mean - real$y
  return(list(x = real$x, y = result))
}


includeLegend <- function() {
  legend(
    "topright",
    legend = c(
      "true density",
      "gamma_MOD",
      "gamma_EXP" ,
      "gamma_Chen_MOD",
      "gamma_Chen_EXP",
      "rectangular",
      "mean MOD_EXP"
    ),
    col = c("black", "red", "blue", "green", "cyan", "grey", "purple"),
    lty = c(1, 1, 1, 1, 1, 1, 1),
    text.width = 7.5
  )
}

includeMSELegend <- function(position = "topright") {
  legend(
    position,
    legend = c(
      "gamma_MOD",
      "gamma_EXP" ,
      "gamma_Chen_MOD",
      "gamma_Chen_EXP",
      "rectangular",
      "mean MOD_EXP"
    ),
    col = c("red", "blue", "green", "cyan", "grey", "purple"),
    lty = c(1, 1, 1, 1, 1, 1),
    text.width = 1,
    cex = 1
  )
}

includeSmallMSELegend <- function(position = "topright") {
  legend(
    position,
    legend = c(
      "gamma_MOD",
      "gamma_EXP" ,
      "gamma_Chen_MOD",
      "gamma_Chen_EXP",
      "rectangular",
      "mean MOD_EXP"
    ),
    col = c("red", "blue", "green", "cyan", "grey", "purple"),
    lty = c(1, 1, 1, 1, 1, 1),
    text.width = 0.5,
    cex = 1
  )
}



#Initilizing Parameters for the investigation
sizeOfSimulatedVector <- 10    #Number of simulated values
grainSize <-
  1000             #how fine the aproximation shall be displayed
numberOfBiasSamples <-
  100   #How many samples are used to calculate the bias


#Define the distrtibutions for the simulation
#_____________________________________________

#_____________________________________________________________________________
#
#Doubled normal distribution
#
#_____________________________________________________________________________


#Generating sample
#set.seed(7685865) #first sample
#set.seed(3457676) #second sample
set.seed(6795999) #third sample
sdNomr = 1
normSimulation <- abs(rnorm(sizeOfSimulatedVector, sd = sdNomr))


#Bandwidth calculation
realDensity <-
  curve(2 * dnorm(x, sd = sdNomr), n = grainSize, xlim = c(0, 4))
h <-
  curve(
    CPP_CrossValidation_rectangular(
      bandwidth = x,
      xValues = realDensity$x,
      X_i = normSimulation
    ),
    xlim = c(0.1, 10),
    n = 1000,
    col = 'grey'
  )
h_rectangular <- arg_min(h$x, h$y)
h <-
  curve(
    CPP_CrossValidation_gamma_MOD(
      bandwidth = x,
      xValues = realDensity$x,
      X_i = normSimulation
    ),
    xlim = c(0.01, 2),
    n = 1000,
    add = TRUE,
    col = 'red'
  )
h_gamma_MOD <- arg_min(h$x, h$y)
h <-
  curve(
    CPP_CrossValidation_gamma_EXP(
      bandwidth = x,
      xValues = realDensity$x,
      X_i = normSimulation
    ),
    xlim = c(0.01, 2),
    n = 1000,
    add = TRUE,
    col = 'blue'
  )
h_gamma_EXP <- arg_min(h$x, h$y)
h <-
  curve(
    CPP_CrossValidation_gammaChen_MOD(
      bandwidth = x,
      xValues = realDensity$x,
      X_i = normSimulation
    ),
    xlim = c(0.01, 2),
    n = 1000,
    add = TRUE,
    col = 'green'
  )
h_gammaChen_MOD <- arg_min(h$x, h$y)
h <-
  curve(
    CPP_CrossValidation_gammaChen_EXP(
      bandwidth = x,
      xValues = realDensity$x,
      X_i = normSimulation
    ),
    xlim = c(0.01, 2),
    n = 1000,
    add = TRUE,
    col = 'cyan'
  )
h_gammaChen_EXP <- arg_min(h$x, h$y)

#Calculation and visualisation
plot(
  x = normSimulation,
  y = rep(0, length(normSimulation)) ,
  xlim = c(0, 4),
  xlab = "",
  ylab = "",
  ylim = c(0, 0.8)
)
realDensity <-
  curve(2 * dnorm(x, sd = sdNomr), n = grainSize, add = TRUE)
result_CPP_KDE_gamma_MOD <-
  curve(
    CPP_KDE_gamma_MOD(x, normSimulation, h_gamma_MOD),
    n = grainSize,
    add = TRUE,
    col = 'red'
  )
result_CPP_KDE_gamma_EXP <-
  curve(
    CPP_KDE_gamma_EXP(x, normSimulation, h_gamma_EXP),
    n = grainSize,
    add = TRUE,
    col = 'blue'
  )
result_CPP_KDE_gammaChen_MOD <-
  curve(
    CPP_KDE_gammaChen_MOD(x, normSimulation, h_gammaChen_MOD),
    n = grainSize,
    add = TRUE,
    col = 'green'
  )
result_CPP_KDE_gammaChen_EXP <-
  curve(
    CPP_KDE_gammaChen_EXP(x, normSimulation, h_gammaChen_EXP),
    n = grainSize,
    add = TRUE,
    col = 'cyan'
  )
result_CPP_KDE_rectangular <-
  curve(
    CPP_KDE_rectanuglar(x, normSimulation, h_rectangular, lowerBoundary = 0),
    n = grainSize,
    add = TRUE,
    col = 'grey',
    type = "s"
  )
result_CPP_KDE_gamma_MOD_EXP <-
  curve((0.5 * (
    CPP_KDE_gamma_MOD(x, normSimulation, h_gamma_MOD) + CPP_KDE_gamma_EXP(x, normSimulation, h_gamma_EXP)
  )),
  n = grainSize,
  add = TRUE,
  col = 'purple')
includeLegend()

#Analysing the results
#_____________________________________________________________________________

#MSE
MSE_result_gamma_MOD <-
  curve((
    2 * dnorm(x, sd = sdNomr) - CPP_KDE_gamma_MOD(x, normSimulation, h_gamma_MOD)
  ) ^
    2,
  n = grainSize,
  col = 'red',
  xlim = c(0, 4),
  ylim = c(0, 0.0005),
  ylab = "MSE(x)"
  )
MSE_result_gamma_EXP <-
  curve((
    2 * dnorm(x, sd = sdNomr) - CPP_KDE_gamma_EXP(x, normSimulation, h_gamma_EXP)
  ) ^
    2,
  n = grainSize,
  add = TRUE,
  col = 'blue'
  )
MSE_result_gammaChen_MOD <-
  curve((
    2 * dnorm(x, sd = sdNomr) - CPP_KDE_gammaChen_MOD(x, normSimulation, h_gammaChen_MOD)
  ) ^ 2,
  n = grainSize,
  add = TRUE,
  col = 'green'
  )
MSE_result_gammaChen_EXP <-
  curve((
    2 * dnorm(x, sd = sdNomr) - CPP_KDE_gammaChen_EXP(x, normSimulation, h_gammaChen_EXP)
  ) ^ 2,
  n = grainSize,
  add = TRUE,
  col = 'cyan'
  )
MSE_result_rectangular <-
  curve((
    2 * dnorm(x, sd = sdNomr) - CPP_KDE_rectanuglar(x, normSimulation, h_rectangular, lowerBoundary = 0)
  ) ^ 2,
  n = grainSize,
  add = TRUE,
  col = 'grey'
  )
MSE_result_MOD_EXP <-
  curve((
    2 * dnorm(x, sd = sdNomr) - 0.5 * (
      CPP_KDE_gamma_MOD(x, normSimulation, h_gamma_MOD) + CPP_KDE_gamma_EXP(x, normSimulation, h_gamma_EXP)
    )
  ) ^ 2,
  n = grainSize,
  add = TRUE,
  col = 'purple')
includeMSELegend()

#MISE
MISE_result <-
  MISE(real = realDensity$y, predicted = result_CPP_KDE_gamma_MOD$y)
MISE_result <-
  c(MISE_result,
    MISE(real = realDensity$y, predicted = result_CPP_KDE_gamma_EXP$y))
MISE_result <-
  c(MISE_result,
    MISE(real = realDensity$y, predicted = result_CPP_KDE_gammaChen_MOD$y))
MISE_result <-
  c(MISE_result,
    MISE(real = realDensity$y, predicted = result_CPP_KDE_gammaChen_EXP$y))
MISE_result <-
  c(MISE_result,
    MISE(real = realDensity$y, predicted = result_CPP_KDE_rectangular$y))
MISE_result <-
  c(MISE_result,
    MISE(real = realDensity$y, predicted = result_CPP_KDE_gamma_MOD_EXP$y))


#Create dataframe
result <- data.frame(MISE_result)

#Add to dataframe
result <- cbind(result, MISE_result)
rownames(result) <-
  c(
    "gamma_MOD",
    "gamma_EXP",
    "gammaChen_MOD",
    "gammaChen_EXP",
    "rectangular",
    "gamma_MOD&EXP"
  )
colnames(result) <- c("n=10", "n=100", "n=1000", "n=10000")

xtable(result, display = rep("e", 5), digits = 3)

#Calculate general mean bandwidth
h_rectangular <- 0                             #Reset the varibales
h_gamma_MOD <- 0
h_gamma_EXP <- 0
h_gammaChen_MOD <- 0
h_gammaChen_EXP <- 0
realDensity <-
  curve(2 * dnorm(x, sd = sdNomr), n = grainSize, xlim = c(0, 20)) #Define the base for the calculation
sizeOfSimulatedVector <- 100                    #Set size of samples
samples <-
  generateSamples_doubledNorm(10, sd = sdNomr)   #Generate 10 samples for the bandwidth calulation

for (i in 1:length(samples)) {
  print(i)
  h <-
    curve(
      CPP_CrossValidation_rectangular(
        bandwidth = x,
        xValues = realDensity$x,
        X_i = as.numeric(unlist(samples[i]))
      ),
      xlim = c(0.1, 10),
      n = 1000,
      col = 'grey'
    )
  h_rectangular <- h_rectangular + arg_min(h$x, h$y)
  h <-
    curve(
      CPP_CrossValidation_gamma_MOD(
        bandwidth = x,
        xValues = realDensity$x,
        X_i = as.numeric(unlist(samples[i]))
      ),
      xlim = c(0.01, 2),
      n = 1000,
      add = TRUE,
      col = 'red'
    )
  h_gamma_MOD <- h_gamma_MOD + arg_min(h$x, h$y)
  h <-
    curve(
      CPP_CrossValidation_gamma_EXP(
        bandwidth = x,
        xValues = realDensity$x,
        X_i = as.numeric(unlist(samples[i]))
      ),
      xlim = c(0.01, 2),
      n = 1000,
      add = TRUE,
      col = 'blue'
    )
  h_gamma_EXP <- h_gamma_EXP + arg_min(h$x, h$y)
  h <-
    curve(
      CPP_CrossValidation_gammaChen_MOD(
        bandwidth = x,
        xValues = realDensity$x,
        X_i = as.numeric(unlist(samples[i]))
      ),
      xlim = c(0.01, 2),
      n = 1000,
      add = TRUE,
      col = 'green'
    )
  h_gammaChen_MOD <- h_gammaChen_MOD + arg_min(h$x, h$y)
  h <-
    curve(
      CPP_CrossValidation_gammaChen_EXP(
        bandwidth = x,
        xValues = realDensity$x,
        X_i = as.numeric(unlist(samples[i]))
      ),
      xlim = c(0.01, 2),
      n = 1000,
      add = TRUE,
      col = 'cyan'
    )
  h_gammaChen_EXP <- h_gammaChen_EXP + arg_min(h$x, h$y)
}

h_rectangular <- h_rectangular / length(samples)
h_gamma_MOD <- h_gamma_MOD / length(samples)
h_gamma_EXP <- h_gamma_EXP / length(samples)
h_gammaChen_MOD <- h_gammaChen_MOD / length(samples)
h_gammaChen_EXP <- h_gammaChen_EXP / length(samples)


#Bias analysation
samples <-
  generateSamples_doubledNorm(numberOfBiasSamples, sd = sdNomr)
result_Bias_gamma_MOD <-
  Bias_gamma_MOD(real = realDensity, samples = samples)
result_Bias_gamma_EXP <-
  Bias_gamma_EXP(real = realDensity, samples = samples)
result_Bias_gammaChen_MOD <-
  Bias_gammaChen_MOD(real = realDensity, samples = samples)
result_Bias_gammaChen_EXP <-
  Bias_gammaChen_EXP(real = realDensity, samples = samples)
result_Bias_rectangular <-
  Bias_rectangular(real = realDensity, samples = samples)
result_Bias_gamma_MOD_EXP <-
  Bias_gamma_MOD_EXP(real = realDensity, samples = samples)
#Plot bias
plot(
  x = result_Bias_gamma_MOD$x ,
  y = result_Bias_gamma_MOD$y,
  type = 'l',
  col = 'red',
  ylim = c(-0.02, 0.01),
  xlab = "x",
  ylab = "Bias(x)"
)
lines(
  x = result_Bias_gamma_EXP$x ,
  y = result_Bias_gamma_EXP$y,
  type = 'l',
  col = 'blue'
)
lines(
  x = result_Bias_gammaChen_MOD$x ,
  y = result_Bias_gammaChen_MOD$y,
  type = 'l',
  col = 'green'
)
lines(
  x = result_Bias_gammaChen_EXP$x ,
  y = result_Bias_gammaChen_EXP$y,
  type = 'l',
  col = 'cyan'
)
lines(
  x = result_Bias_rectangular$x ,
  y = result_Bias_rectangular$y,
  type = 'l',
  col = 'grey'
)
lines(
  x = result_Bias_gamma_MOD_EXP$x ,
  y = result_Bias_gamma_MOD_EXP$y,
  type = 'l',
  col = 'purple'
)
abline(h = 0, lty = 3)
includeMSELegend("bottomright")



#_____________________________________________________________________________
#
#Gamma distribution
#
#_____________________________________________________________________________

#Calculation and visualisation
sizeOfSimulatedVector <- 100    #Number of simulated values
#set.seed(76865) #first sample
#set.seed(3457676) #second sample
set.seed(677659) #third sample
p = 1
rate <- 0.7
gammaSimulation <-
  rgamma(sizeOfSimulatedVector, shape = p, rate = rate)

#Bandwidth calculation
realDensity <-
  curve(dgamma(x, shape = p, rate = rate),
        n = grainSize,
        xlim = c(0, 20),ylim=c(0,0.17),xlab="",ylab="")
h <-
  curve(
    CPP_CrossValidation_rectangular(
      bandwidth = x,
      xValues = realDensity$x,
      X_i = gammaSimulation
    ),
    xlim = c(0.1, 8),
    n = 1000,
    col = 'grey'
  )
h_rectangular <- arg_min(h$x, h$y)
h <-
  curve(
    CPP_CrossValidation_gamma_MOD(
      bandwidth = x,
      xValues = realDensity$x,
      X_i = gammaSimulation
    ),
    xlim = c(0.01, 8),
    n = 1000,
    add = TRUE,
    col = 'red'
  )
h_gamma_MOD <- arg_min(h$x, h$y)
h <-
  curve(
    CPP_CrossValidation_gamma_EXP(
      bandwidth = x,
      xValues = realDensity$x,
      X_i = gammaSimulation
    ),
    xlim = c(0.01, 8),
    n = 1000,
    add = TRUE,
    col = 'blue'
  )
h_gamma_EXP <- arg_min(h$x, h$y)
h <-
  curve(
    CPP_CrossValidation_gammaChen_MOD(
      bandwidth = x,
      xValues = realDensity$x,
      X_i = gammaSimulation
    ),
    xlim = c(0.01, 8),
    n = 1000,
    add = TRUE,
    col = 'green'
  )
h_gammaChen_MOD <- arg_min(h$x, h$y)
h <-
  curve(
    CPP_CrossValidation_gammaChen_EXP(
      bandwidth = x,
      xValues = realDensity$x,
      X_i = gammaSimulation
    ),
    xlim = c(0.01, 8),
    n = 1000,
    add = TRUE,
    col = 'cyan'
  )
h_gammaChen_EXP <- arg_min(h$x, h$y)

plot(
  x = gammaSimulation,
  y = rep(0, length(gammaSimulation)) ,
  xlim = c(0, 3),
  xlab = "",
  ylab = "",
  ylim = c(0, 0.12)
)
realDensity <-
  curve(dgamma(x, shape = p, rate = rate), n = grainSize, add = TRUE)
result_CPP_KDE_gamma_MOD <-
  curve(
    CPP_KDE_gamma_MOD(x, gammaSimulation, h_gamma_MOD),
    n = grainSize,
    add = TRUE,
    col = 'red'
  )
result_CPP_KDE_gamma_EXP <-
  curve(
    CPP_KDE_gamma_EXP(x, gammaSimulation, h_gamma_EXP),
    n = grainSize,
    add = TRUE,
    col = 'blue'
  )
result_CPP_KDE_gammaChen_MOD <-
  curve(
    CPP_KDE_gammaChen_MOD(x, gammaSimulation, h_gammaChen_MOD),
    n = grainSize,
    add = TRUE,
    col = 'green'
  )
result_CPP_KDE_gammaChen_EXP <-
  curve(
    CPP_KDE_gammaChen_EXP(x, gammaSimulation, h_gammaChen_EXP),
    n = grainSize,
    add = TRUE,
    col = 'cyan'
  )
result_CPP_KDE_rectangular <-
  curve(
    CPP_KDE_rectanuglar(x, gammaSimulation, h_rectangular, lowerBoundary = 0),
    n = grainSize,
    add = TRUE,
    col = 'grey',
    type = "s"
  )
result_CPP_KDE_gamma_MOD_EXP <-
  curve((0.5 * (
    CPP_KDE_gamma_MOD(x, gammaSimulation, h_gamma_MOD) + CPP_KDE_gamma_EXP(x, gammaSimulation, h_gamma_EXP)
  )),
  n = grainSize,
  add = TRUE,
  col = 'purple')
includeLegend()


#Analysing the results
#_____________________________________________________________________________

#MSE
MSE_result_gamma_MOD <-
  curve((
    dgamma(x, shape = p, rate = rate) - CPP_KDE_gamma_MOD(x, gammaSimulation, h_gamma_MOD)
  ) ^ 2,
  n = grainSize,
  col = 'red',
  xlim = c(0, 4),
  ylab = "MSE(x)",
  ylim = c(0, 0.00085)
  )
MSE_result_gamma_EXP <-
  curve((
    dgamma(x, shape = p, rate = rate) - CPP_KDE_gamma_EXP(x, gammaSimulation, h_gamma_EXP)
  ) ^ 2,
  n = grainSize,
  add = TRUE,
  col = 'blue'
  )
MSE_result_gammaChen_MOD <-
  curve((
    dgamma(x, shape = p, rate = rate) - CPP_KDE_gammaChen_MOD(x, gammaSimulation, h_gammaChen_MOD)
  ) ^ 2,
  n = grainSize,
  add = TRUE,
  col = 'green'
  )
MSE_result_gammaChen_EXP <-
  curve((
    dgamma(x, shape = p, rate = rate) - CPP_KDE_gammaChen_EXP(x, gammaSimulation, h_gammaChen_EXP)
  ) ^ 2,
  n = grainSize,
  add = TRUE,
  col = 'cyan'
  )
MSE_result_rectangular <-
  curve((
    dgamma(x, shape = p, rate = rate) - CPP_KDE_rectanuglar(x, gammaSimulation, h_rectangular, lowerBoundary = 0)
  ) ^ 2,
  n = grainSize,
  add = TRUE,
  col = 'grey'
  )
MSE_result_MOD_EXP <-
  curve((
    dgamma(x, shape = p, rate = rate) - 0.5 * (
      CPP_KDE_gamma_MOD(x, gammaSimulation, h_gamma_MOD) + CPP_KDE_gamma_EXP(x, gammaSimulation, h_gamma_EXP)
    )
  ) ^ 2,
  n = grainSize,
  add = TRUE,
  col = 'purple')
includeMSELegend()


#MISE
MISE_result <-
  MISE(real = realDensity$y, predicted = result_CPP_KDE_gamma_MOD$y)
MISE_result <-
  c(MISE_result,
    MISE(real = realDensity$y, predicted = result_CPP_KDE_gamma_EXP$y))
MISE_result <-
  c(MISE_result,
    MISE(real = realDensity$y, predicted = result_CPP_KDE_gammaChen_MOD$y))
MISE_result <-
  c(MISE_result,
    MISE(real = realDensity$y, predicted = result_CPP_KDE_gammaChen_EXP$y))
MISE_result <-
  c(MISE_result,
    MISE(real = realDensity$y, predicted = result_CPP_KDE_rectangular$y))
MISE_result <-
  c(MISE_result,
    MISE(real = realDensity$y, predicted = result_CPP_KDE_gamma_MOD_EXP$y))

#Create dataframe
result <- data.frame(MISE_result)

#Add to dataframe
result <- cbind(result, MISE_result)
rownames(result) <-
  c(
    "gamma_MOD",
    "gamma_EXP",
    "gammaChen_MOD",
    "gammaChen_EXP",
    "rectangular",
    "gamma_MOD&EXP"
  )
colnames(result) <- c("n=10", "n=100", "n=1000", "n=10000")

xtable(result, display = rep("e", 5), digits = 3)

#Calculate general mean bandwidth
h_rectangular <- 0                             #Reset the varibales
h_gamma_MOD <- 0
h_gamma_EXP <- 0
h_gammaChen_MOD <- 0
h_gammaChen_EXP <- 0
realDensity <-
  curve(dgamma(x, shape = p, rate = rate),
        n = grainSize,
        xlim = c(0, 20)) #Define the base for the calculation
sizeOfSimulatedVector <- 100                    #Set size of samples
samples <-
  generateSamples_gamma(10, p = p, b = rate)   #Generate 10 samples for the bandwidth calulation

for (i in 1:length(samples)) {
  print(i)
  h <-
    curve(
      CPP_CrossValidation_rectangular(
        bandwidth = x,
        xValues = realDensity$x,
        X_i = as.numeric(unlist(samples[i]))
      ),
      xlim = c(0.1, 8),
      n = 1000,
      col = 'grey'
    )
  h_rectangular <- h_rectangular + arg_min(h$x, h$y)
  h <-
    curve(
      CPP_CrossValidation_gamma_MOD(
        bandwidth = x,
        xValues = realDensity$x,
        X_i = as.numeric(unlist(samples[i]))
      ),
      xlim = c(0.01, 8),
      n = 1000,
      add = TRUE,
      col = 'red'
    )
  h_gamma_MOD <- h_gamma_MOD + arg_min(h$x, h$y)
  h <-
    curve(
      CPP_CrossValidation_gamma_EXP(
        bandwidth = x,
        xValues = realDensity$x,
        X_i = as.numeric(unlist(samples[i]))
      ),
      xlim = c(0.01, 8),
      n = 1000,
      add = TRUE,
      col = 'blue'
    )
  h_gamma_EXP <- h_gamma_EXP + arg_min(h$x, h$y)
  h <-
    curve(
      CPP_CrossValidation_gammaChen_MOD(
        bandwidth = x,
        xValues = realDensity$x,
        X_i = as.numeric(unlist(samples[i]))
      ),
      xlim = c(0.01, 8),
      n = 1000,
      add = TRUE,
      col = 'green'
    )
  h_gammaChen_MOD <- h_gammaChen_MOD + arg_min(h$x, h$y)
  h <-
    curve(
      CPP_CrossValidation_gammaChen_EXP(
        bandwidth = x,
        xValues = realDensity$x,
        X_i = as.numeric(unlist(samples[i]))
      ),
      xlim = c(0.01, 8),
      n = 1000,
      add = TRUE,
      col = 'cyan'
    )
  h_gammaChen_EXP <- h_gammaChen_EXP + arg_min(h$x, h$y)
}

h_rectangular <- h_rectangular / length(samples)
h_gamma_MOD <- h_gamma_MOD / length(samples)
h_gamma_EXP <- h_gamma_EXP / length(samples)
h_gammaChen_MOD <- h_gammaChen_MOD / length(samples)
h_gammaChen_EXP <- h_gammaChen_EXP / length(samples)


#Bias analysation
samples <- generateSamples_gamma(numberOfBiasSamples, p, b = rate)
result_Bias_gamma_MOD <-
  Bias_gamma_MOD(real = realDensity, samples = samples)
result_Bias_gamma_EXP <-
  Bias_gamma_EXP(real = realDensity, samples = samples)
result_Bias_gammaChen_MOD <-
  Bias_gammaChen_MOD(real = realDensity, samples = samples)
result_Bias_gammaChen_EXP <-
  Bias_gammaChen_EXP(real = realDensity, samples = samples)
result_Bias_rectangular <-
  Bias_rectangular(real = realDensity, samples = samples)
result_Bias_gamma_MOD_EXP <-
  Bias_gamma_MOD_EXP(real = realDensity, samples = samples)
#Plot bias
plot(
  x = result_Bias_gamma_MOD$x ,
  y = result_Bias_gamma_MOD$y,
  type = 'l',
  col = 'red',
  ylim = c(-0.025, 0.03),
  xlab = "x",
  ylab = "Bias(x)"
)
lines(
  x = result_Bias_gamma_EXP$x ,
  y = result_Bias_gamma_EXP$y,
  type = 'l',
  col = 'blue'
)
lines(
  x = result_Bias_gammaChen_MOD$x ,
  y = result_Bias_gammaChen_MOD$y,
  type = 'l',
  col = 'green'
)
lines(
  x = result_Bias_gammaChen_EXP$x ,
  y = result_Bias_gammaChen_EXP$y,
  type = 'l',
  col = 'cyan'
)
lines(
  x = result_Bias_rectangular$x ,
  y = result_Bias_rectangular$y,
  type = 'l',
  col = 'grey'
)
lines(
  x = result_Bias_gamma_MOD_EXP$x ,
  y = result_Bias_gamma_MOD_EXP$y,
  type = 'l',
  col = 'purple'
)
abline(h = 0, lty = 3)
includeMSELegend()

#_____________________________________________________________________________
#
#Sinus distribution
#
#_____________________________________________________________________________

#Functions needed
sinusDistribution <- function(x) {
  result <- (abs(sin(x / 2)) / 8)
  for (v in 1:length(result)) {
    if (x[v] < 0) {
      result[v] <- 0
    }
    if (x[v] > (pi * 4)) {
      result[v] <- 0
    }
  }
  return(result)
}

generateSinusDistribution <- function() {
  result <- seq(1:sizeOfSimulatedVector)
  for (i in 1:sizeOfSimulatedVector) {
    t <- 1
    x <- runif(1, 0, 13)
    while (t > sinusDistribution(x)) {
      x <- runif(1, 0, 13)
      t <- runif(1, 0, 0.125)
    }
    result[i] <- x
  }
  return(result)
}

#Calculation and visualisation
sizeOfSimulatedVector <- 10    #Number of simulated values
#set.seed(768654) #first sample
#set.seed(3457676) #second sample
set.seed(677659) #third sample
sinusSimulation <- generateSinusDistribution()

#Bandwidth calculation
realDensity <-
  curve(sinusDistribution(x), xlim = c(0, 20), n = grainSize)
h <-
  curve(
    CPP_CrossValidation_rectangular(
      bandwidth = x,
      xValues = realDensity$x,
      X_i = sinusSimulation
    ),
    xlim = c(0.1, 20),
    n = 1000,
    col = 'grey'
  )
h_rectangular <- arg_min(h$x, h$y)
h <-
  curve(
    CPP_CrossValidation_gamma_MOD(
      bandwidth = x,
      xValues = realDensity$x,
      X_i = sinusSimulation
    ),
    xlim = c(0.01, 20),
    n = 1000,
    add = TRUE,
    col = 'red'
  )
h_gamma_MOD <- arg_min(h$x, h$y)
h <-
  curve(
    CPP_CrossValidation_gamma_EXP(
      bandwidth = x,
      xValues = realDensity$x,
      X_i = sinusSimulation
    ),
    xlim = c(0.01, 20),
    n = 1000,
    add = TRUE,
    col = 'blue'
  )
h_gamma_EXP <- arg_min(h$x, h$y)
h <-
  curve(
    CPP_CrossValidation_gammaChen_MOD(
      bandwidth = x,
      xValues = realDensity$x,
      X_i = sinusSimulation
    ),
    xlim = c(0.01, 20),
    n = 1000,
    add = TRUE,
    col = 'green'
  )
h_gammaChen_MOD <- arg_min(h$x, h$y)
h <-
  curve(
    CPP_CrossValidation_gammaChen_EXP(
      bandwidth = x,
      xValues = realDensity$x,
      X_i = sinusSimulation
    ),
    xlim = c(0.01, 20),
    n = 1000,
    add = TRUE,
    col = 'cyan'
  )
h_gammaChen_EXP <- arg_min(h$x, h$y)


plot(
  x = sinusSimulation,
  y = rep(0, length(sinusSimulation)) ,
  xlim = c(0, 1),
  xlab = "",
  ylab = "",
  ylim = c(0, 0.2)
)
realDensity <- curve(sinusDistribution(x),
                     add = TRUE,
                     n = grainSize)
result_CPP_KDE_gamma_MOD <-
  curve(
    CPP_KDE_gamma_MOD(x, sinusSimulation, h_gamma_MOD),
    n = grainSize,
    add = TRUE,
    col = 'red'
  )
result_CPP_KDE_gamma_EXP <-
  curve(
    CPP_KDE_gamma_EXP(x, sinusSimulation, h_gamma_EXP),
    n = grainSize,
    add = TRUE,
    col = 'blue'
  )
result_CPP_KDE_gammaChen_MOD <-
  curve(
    CPP_KDE_gammaChen_MOD(x, sinusSimulation, h_gammaChen_MOD),
    n = grainSize,
    add = TRUE,
    col = 'green'
  )
result_CPP_KDE_gammaChen_EXP <-
  curve(
    CPP_KDE_gammaChen_EXP(x, sinusSimulation, h_gammaChen_EXP),
    n = grainSize,
    add = TRUE,
    col = 'cyan'
  )
result_CPP_KDE_rectangular <-
  curve(
    CPP_KDE_rectanuglar(x, sinusSimulation, h_rectangular, lowerBoundary = 0),
    n = grainSize,
    add = TRUE,
    col = 'grey',
    type = "s"
  )
result_CPP_KDE_gamma_MOD_EXP <-
  curve((0.5 * (
    CPP_KDE_gamma_MOD(x, sinusSimulation, h_gamma_MOD) + CPP_KDE_gamma_EXP(x, sinusSimulation, h_gamma_EXP)
  )),
  n = grainSize,
  add = TRUE,
  col = 'purple')
includeLegend()

#Analysing the results
#_____________________________________________________________________________

MSE_result_gamma_MOD <-
  curve((
    sinusDistribution(x) - CPP_KDE_gamma_MOD(x, sinusSimulation, h_gamma_MOD)
  ) ^ 2,
  n = grainSize,
  col = 'red',
  xlim = c(0, 1),
  ylim = c(0, 0.001),
  ylab = "MSE(x)"
  )
MSE_result_gamma_EXP <-
  curve((
    sinusDistribution(x) - CPP_KDE_gamma_EXP(x, sinusSimulation, h_gamma_EXP)
  ) ^ 2,
  n = grainSize,
  add = TRUE,
  col = 'blue'
  )
MSE_result_gammaChen_MOD <-
  curve((
    sinusDistribution(x) - CPP_KDE_gammaChen_MOD(x, sinusSimulation, h_gammaChen_MOD)
  ) ^ 2,
  n = grainSize,
  add = TRUE,
  col = 'green'
  )
MSE_result_gammaChen_EXP <-
  curve((
    sinusDistribution(x) - CPP_KDE_gammaChen_EXP(x, sinusSimulation, h_gammaChen_EXP)
  ) ^ 2,
  n = grainSize,
  add = TRUE,
  col = 'cyan'
  )
MSE_result_rectangular <-
  curve((
    sinusDistribution(x) - CPP_KDE_rectanuglar(x, sinusSimulation, h_rectangular, lowerBoundary = 0)
  ) ^ 2,
  n = grainSize,
  add = TRUE,
  col = 'grey'
  )
MSE_result_MOD_EXP <-
  curve((
    sinusDistribution(x) - 0.5 * (
      CPP_KDE_gamma_MOD(x, sinusSimulation, h_gamma_MOD) + CPP_KDE_gamma_EXP(x, sinusSimulation, h_gamma_EXP)
    )
  ) ^ 2,
  n = grainSize,
  add = TRUE,
  col = 'purple'
  )
includeSmallMSELegend()


#MISE
MISE_result <-
  MISE(real = realDensity$y, predicted = result_CPP_KDE_gamma_MOD$y)
MISE_result <-
  c(MISE_result,
    MISE(real = realDensity$y, predicted = result_CPP_KDE_gamma_EXP$y))
MISE_result <-
  c(MISE_result,
    MISE(real = realDensity$y, predicted = result_CPP_KDE_gammaChen_MOD$y))
MISE_result <-
  c(MISE_result,
    MISE(real = realDensity$y, predicted = result_CPP_KDE_gammaChen_EXP$y))
MISE_result <-
  c(MISE_result,
    MISE(real = realDensity$y, predicted = result_CPP_KDE_rectangular$y))
MISE_result <-
  c(MISE_result,
    MISE(real = realDensity$y, predicted = result_CPP_KDE_gamma_MOD_EXP$y))



#Create dataframe
result <- data.frame(MISE_result)

#Add to dataframe
result <- cbind(result, MISE_result)
rownames(result) <-
  c(
    "gamma_MOD",
    "gamma_EXP",
    "gammaChen_MOD",
    "gammaChen_EXP",
    "rectangular",
    "gamma_MOD&EXP"
  )
colnames(result) <- c("n=10", "n=100", "n=1000", "n=10000")

xtable(result, display = rep("e", 5), digits = 3)

#Calculate general mean bandwidth
h_rectangular <- 0                             #Reset the varibales
h_gamma_MOD <- 0
h_gamma_EXP <- 0
h_gammaChen_MOD <- 0
h_gammaChen_EXP <- 0
realDensity <-
  curve(sinusDistribution(x), xlim = c(0, 20), n = grainSize) #Define the base for the calculation
sizeOfSimulatedVector <- 100                    #Set size of samples
samples <-
  generateSamples_sinus(10)   #Generate 10 samples for the bandwidth calulation

for (i in 1:length(samples)) {
  print(i)
  h <-
    curve(
      CPP_CrossValidation_rectangular(
        bandwidth = x,
        xValues = realDensity$x,
        X_i = as.numeric(unlist(samples[i]))
      ),
      xlim = c(0.1, 8),
      n = 1000,
      col = 'grey'
    )
  h_rectangular <- h_rectangular + arg_min(h$x, h$y)
  h <-
    curve(
      CPP_CrossValidation_gamma_MOD(
        bandwidth = x,
        xValues = realDensity$x,
        X_i = as.numeric(unlist(samples[i]))
      ),
      xlim = c(0.01, 8),
      n = 1000,
      add = TRUE,
      col = 'red'
    )
  h_gamma_MOD <- h_gamma_MOD + arg_min(h$x, h$y)
  h <-
    curve(
      CPP_CrossValidation_gamma_EXP(
        bandwidth = x,
        xValues = realDensity$x,
        X_i = as.numeric(unlist(samples[i]))
      ),
      xlim = c(0.01, 8),
      n = 1000,
      add = TRUE,
      col = 'blue'
    )
  h_gamma_EXP <- h_gamma_EXP + arg_min(h$x, h$y)
  h <-
    curve(
      CPP_CrossValidation_gammaChen_MOD(
        bandwidth = x,
        xValues = realDensity$x,
        X_i = as.numeric(unlist(samples[i]))
      ),
      xlim = c(0.01, 8),
      n = 1000,
      add = TRUE,
      col = 'green'
    )
  h_gammaChen_MOD <- h_gammaChen_MOD + arg_min(h$x, h$y)
  h <-
    curve(
      CPP_CrossValidation_gammaChen_EXP(
        bandwidth = x,
        xValues = realDensity$x,
        X_i = as.numeric(unlist(samples[i]))
      ),
      xlim = c(0.01, 8),
      n = 1000,
      add = TRUE,
      col = 'cyan'
    )
  h_gammaChen_EXP <- h_gammaChen_EXP + arg_min(h$x, h$y)
}

h_rectangular <- h_rectangular / length(samples)
h_gamma_MOD <- h_gamma_MOD / length(samples)
h_gamma_EXP <- h_gamma_EXP / length(samples)
h_gammaChen_MOD <- h_gammaChen_MOD / length(samples)
h_gammaChen_EXP <- h_gammaChen_EXP / length(samples)



#Bias analysation
samples <- generateSamples_sinus(numberOfBiasSamples)
result_Bias_gamma_MOD <-
  Bias_gamma_MOD(real = realDensity, samples = samples)
result_Bias_gamma_EXP <-
  Bias_gamma_EXP(real = realDensity, samples = samples)
result_Bias_gammaChen_MOD <-
  Bias_gammaChen_MOD(real = realDensity, samples = samples)
result_Bias_gammaChen_EXP <-
  Bias_gammaChen_EXP(real = realDensity, samples = samples)
result_Bias_rectangular <-
  Bias_rectangular(real = realDensity, samples = samples)
result_Bias_gamma_MOD_EXP <-
  Bias_gamma_MOD_EXP(real = realDensity, samples = samples)
#Plot bias
plot(
  x = result_Bias_gamma_MOD$x ,
  y = result_Bias_gamma_MOD$y,
  type = 'l',
  col = 'red',
  ylim = c(-0.01, 0.07),
  xlab = "x",
  ylab = "Bias(x)"
)
lines(
  x = result_Bias_gamma_EXP$x ,
  y = result_Bias_gamma_EXP$y,
  type = 'l',
  col = 'blue'
)
lines(
  x = result_Bias_gammaChen_MOD$x ,
  y = result_Bias_gammaChen_MOD$y,
  type = 'l',
  col = 'green'
)
lines(
  x = result_Bias_gammaChen_EXP$x ,
  y = result_Bias_gammaChen_EXP$y,
  type = 'l',
  col = 'cyan'
)
lines(
  x = result_Bias_rectangular$x ,
  y = result_Bias_rectangular$y,
  type = 'l',
  col = 'grey'
)
lines(
  x = result_Bias_gamma_MOD_EXP$x ,
  y = result_Bias_gamma_MOD_EXP$y,
  type = 'l',
  col = 'purple'
)
abline(h = 0, lty = 3)
includeSmallMSELegend()

#Speed test
system.time(density())
