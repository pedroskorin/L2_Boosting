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

# Getting functions
# Module for L2 functions

Fitting = function(Y, X, M, v, AICc){
  
  # Defining variables
  ft = mean(Y) # First appearence of function ft
  m = 1 # m index
  
  #X = cbind(X,rep(1,nrow(X)))
  
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
  
  s = c() # vector of choosed predictors
  
  while (m <= M) { # Iterate M times
    u = Y - ft # Calculating error in m
    
    # calculating optimum coefficients
    teta = c()
    
    for (i in X){
      
      b = t(i) %*% u %*% 1/sum(i^2)
      
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
    
    s = append(s, index)
    
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
  
  # Calculating hat matrix
  
  for (j in M:1){
    
    H = as.matrix(X[,s[j]]) %*% (1/sum((X[,s[j]])^2)) %*% t(as.matrix(X[,s[j]]))
    
    if (j == M){
      
      C = diag(nrow(X)) - v * H
      
    } else {
      
      C = C %*% (diag(nrow(X)) - v * H)
      
    }
    
  }
  
  # Calculating Big Beta
  
  B = diag(nrow(X)) - C
  
  # Calculating sigma squared
  
  sigma = nrow(X)^(-1) * sum((Y - ( B %*% Y ))^2)
  
  # Calculating trace
  
  df_m = sum(diag(B))
  
  # Calculating AIC
  
  AIC = log(sigma) + (1 + df_m/nrow(X))/(1 - (df_m + 2)/nrow(X))
  
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

L2_min = function(Y, X, v, AICc){
  
  x = Fitting(Y, X, 2, v , AICc)[[2]]
  IC = c(Fitting(Y, X, 1, v, AICc)[[2]], x)
  i = 3
  
  while(IC[length(IC)-1] > x) {
    
    x = Fitting(Y, X, i, v, AICc)[[2]]
    
    IC = c(IC, x)
    
    i = i + 1
    
  }
  
  return(i-1)
}

L2_boost = function(Y, X, v, AICc) {
  
  # Function for minimizing AIC 
  
  # Optimal M
  m_opt = L2_min(Y, X, v, AICc)
  
  # Fitting
  optimal = Fitting(Y, X, m_opt, v, AICc)
  
  coefficients_in = optimal[[3]][nrow(optimal[[3]]),]
  names(coefficients_in) = names(X)
  
  message = function(coef_in = coefficients_in, kylo = optimal) {
    cat("\n","L2_Boosting", "\n", "\n", "Selected M ","\n")
    print(m_opt)
    cat("\n", "AIC:", kylo[[2]])
    cat("\n", "MSE:",kylo[[1]])
    cat("\n")
  }
  
  op <- list(coefficients = coefficients_in, fitted = optimal[[4]],
             residuals = optimal[[5]], IC = optimal[[2]], MSE = optimal[[1]],
             coef_history = optimal[[3]], a = message)
  
  class(op) <- 'L2_boost'
  
  op
  
}

# Creating function L2_boost:
# INPUT: Y - response ; X - predictors ; M - Iterations ; v - shrinkage parameter (standart 0.1)
# OUTPUT : MSE of algorithm ; times each predictor was choosed
v = 0.1
AICc = F

b = L2_boost(Y, X, v, AICc)

b$a()

# Function for Historical Coefficients Plotting

coef_hist = function(arg) {
  
  matplot(arg$coef_history, type = c("S") ,col = 1:nrow(arg$coef_history), xlab = "Iterations", ylab = "Coefficients")
}

## Example

coef_hist(b)

# Function for MSE printing
print_L2_boost = function(start, end, by = 1, Y, X, v, AICc) {
  
  m_vector = c()
  
  for (i in seq(start, end, by)) {
    m_vector_in = Fitting(Y, X, i, v, AICc)[1]
    
    m_vector = c(m_vector, m_vector_in)
  }
  
  plot(seq(start, end, by), m_vector, type = "l", xlab = "m", ylab = "MSE")
}

## Example

print_L2_boost(0, 250, 10, Y, X, v, AICc)

# Function for AIC printing

print_AIC = function(start, end, by = 1, Y, X, v, AICc) {
  
  m_vector = c()
  
  for (i in seq(start, end, by)) {
    m_vector_in = Fitting(Y, X, i, v, AICc)[2]
    
    m_vector = c(m_vector, m_vector_in)
    print(i/end)
  }
  
  plot(seq(start, end, by), m_vector, type = "l", xlab = "m", ylab = "IC")
}

## Example

print_AIC(1, 500, 50, Y, X, 0.1, F)

