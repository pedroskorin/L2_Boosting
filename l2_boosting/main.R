#L2 BOOSTING ALGORITHM
#
#Implementation of the L2 Boosting algorithm
#
#theoretical source: "Jing Zang, 2014, Forecasting Aggregates with Disaggregate Variables: 
#Does boosting help to select the most informative predictors?"
#
#Authors: guislinden, pedroskorin
#UFRGS
#
#2019-6-19

# Getting data from github without the first column
data = read.csv("https://raw.githubusercontent.com/pedroskorin/L2_Boosting/master/l2_boosting/dados.csv")[,-1]

# Selecting data
Y = data$Y # Response variable

X = data[,-ncol(data)] # Predictors variables

# Creating function L2_boost:
# INPUT: Y - response ; X - predictors ; M - Iterations ; v - shrinkage parameter (standart 0.1)
# OUTPUT : MSE of algorithm ; times each predictor was choosed
v_out = 0.1
M_out = 9
AICc_out = F
Y_out = Y
X_out = X

L2_boost = function(Y_out = Y, X_out = X, M_out = 9, v_out=0.1, AICc_out = F) {

Fitting = function(Y = Y_out, X = X_out, M = M_out, v = v_out, AICc = AICc_out){
  
  # Defining variables
  ft = mean(Y) # First appearence of function ft
  m = 0 # m index
  
  u = 0 # error vector (Y-ft)
  
  b = 0 # optimum coefficient for each predictor
  teta = c() # vector of optimum coefficient for all possible predictors
  
  sum_squared_resid = 0 # SSR for regression of ut for each predictor
  SSR = c() # list of all SSR
  
  g = c() # vector of selected predictor * optimum coefficient
  coef_control = list() # list of all g's
  index = 0 # index of selected predictor in teta's vector
  
  X_optimum = c() # Vector of selected predictor
  Matrix = c() # Matrix of all selected predictors*coefficient*v
  f_optimum = c() # final prediction vector
  
  choosed_predictors = rep(0,ncol(X)) # controlling choosed predictors
  names(choosed_predictors) = names(X)
  
  coefficient_control = rep(0,ncol(X))
  
  while (m < M) { # Iterate M times
    u = Y - ft # Calculating error in m
    
    # calculating optimum coefficients
    teta = c()
    
    for (i in X){
      b = solve(t(i)%*%i)%*%t(i)%*%u
      
      teta = c(teta, b)
    }
    
    # calculating SSR
    
    SSR = c()
    
    for (i in 1:length(teta)){
      sum_squared_resid = sum((u - teta[i]*X[,i])^2)
      
      SSR = c(SSR, sum_squared_resid)
    }
    
    # Minimizing SSR
    
    index = which.min(SSR)
    
    X_optimum = X[,index]
    
    # Calculating g (vector of selected predictor * optimum coefficient)
    g = teta[index]*X_optimum
    coefficient_control[index] = coefficient_control[index] + teta[index]
    
    choosed_predictors[index] = choosed_predictors[index] + 1
    
    # Updating ft
    ft = ft + v*g
    
    # Iterating
    m = m + 1
    
    coef_control[[m]] <- coefficient_control
    
  }
  
  # Matrix of all coefficients vectors
  Matrix_coef = do.call(rbind, coef_control)
  
  # Final function
  f_optimum = ft
  
  # MSE of L2_boost
  MSE = (sum((f_optimum - Y)^2)/nrow(X)) 
  
  # AIC
  llg = sum(dnorm(Y, mean = ft, sd = MSE, log = TRUE)) # calculate the log likelihood
  k = sum(choosed_predictors) + 2 # quantity of parameters plus SD and Mean
  
  AIC = (-2*llg + 2*k)
  
  # AICc
  
  if (AICc) {
   
    AICC = AIC + ((2*k)^(2) + 2*k)/(length(Y) - k - 1)
    
    specific_IC = AICC
    
  } else {
    
    specific_IC = AIC
    
  }
  
  final_vector = list(MSE, specific_IC, Matrix_coef, ft, Y-ft)
  return(final_vector)
}

# Function for minimizing AIC 

L2_min = function(Y_in = Y_out, X_in = X_out, AICc = AICc_out){
  
  x = Fitting(Y_in, X_in, 2, AICc = AICc)[[2]]
  IC = c(Fitting(Y_in, X_in, 1, AICc = AICc)[[2]], x)
  i = 3
  
  while(IC[length(IC)-1] > x) {
  
    x = Fitting(Y_in , X_in ,i,AICc = AICc)[[2]]
    
    IC = c(IC, x)
    
    i = i + 1
    
  }
  
  return(i-1)
}

optimal = Fitting(Y, X, M = L2_min())

coefficients_in = round(optimal[[3]][nrow(optimal[[3]]),], 4)
names(coefficients_in) = names(X)

message = function(coef_in = coefficients_in, kylo = optimal) {
  cat("\n","L2_Boosting", "\n", "\n", "Coefficients ","\n")
  print(coef_in)
  cat("\n", "AIC:", kylo[[2]])
  cat("\n", "MSE:",kylo[[1]])
}


p = matplot(optimal[[3]],type = c("S") ,col = 1:nrow(optimal[[3]]), xlab = "Iterations", ylab = "Coefficients")
p = ggplot(X, aes(X1, X2)) + geom_point()
op <- list(coefficients = coefficients_in, fitted = optimal[[4]],
           residuals = optimal[[5]], IC = optimal[[2]], MSE = optimal[[1]],
           plot = p, a = message)

#my_func <- function(x, y, z, t){
#  library(ggplot2)
#  df <- data.frame(x, y, z)
#  p <- ggplot(df, aes(x, y)) + geom_point()
#  ds <- sapply(df, summary)
#  
#  op <- list(data = df, plot = p, summary = ds, dale = t)
class(op) <- 'my_list'
#  return(op)
#}

op
#return(my_func(1:5, 5:1, rnorm(5), p))
}

print.my_list <- function(h){
  print(h$a())
}

b = L2_boost()

L2_boost()



# Function for MSE printing
print_L2_boost = function(start, end, by = 1, Y_in = Y, X_in = X) {
  
  m_vector = c()
  
  for (i in seq(start, end, by)) {
    m_vector_in = L2_boost(Y_in,X_in,i)[1]
    
    m_vector = c(m_vector, m_vector_in)
  }
  
  plot(seq(start, end, by), m_vector, type = "l", xlab = "m", ylab = "MSE")
}

# Example
par(mfrow = c(1,2))

print_L2_boost(0,250,10)

# Function for AIC printing

print_AIC = function(start, end, by = 1, Y_in = Y, X_in = X, AICc = F) {
  
  m_vector = c()
  
  for (i in seq(start, end, by)) {
    m_vector_in = L2_boost(Y_in,X_in,i,AICc = AICc)[2]
    
    m_vector = c(m_vector, m_vector_in)
  }
  
  plot(seq(start, end, by), m_vector, type = "l", xlab = "m", ylab = "IC")
}

# Example

print_AIC(0,250,10,AICc = F)



# Example
L2_min(AICc = F)

# Final Example with selected m

L2_boost(Y,X,L2_min(),AICc = F)

