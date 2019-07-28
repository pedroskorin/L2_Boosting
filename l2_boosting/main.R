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

L2_boost = function(Y, X, M, v=0.1, AICc = F){
  
  # Defining variables
  ft = mean(Y) # First appearence of function ft
  m = 0 # m index
  
  u = 0 # error vector (Y-ft)
  
  b = 0 # optimum coefficient for each predictor
  teta = c() # vector of optimum coefficient for all possible predictors
  
  sum_squared_resid = 0 # SSR for regression of ut for each predictor
  SSR = c() # list of all SSR
  
  g = c() # vector of selected predictor * optimum coefficient
  g_optimum = list() # list of all g's
  index = 0 # index of selected predictor in teta's vector
  
  X_optimum = c() # Vector of selected predictor
  Matrix = c() # Matrix of all selected predictors*coefficient*v
  f_optimum = c() # final prediction vector
  
  choosed_predictors = rep(0,ncol(X)) # controlling choosed predictors
  names(choosed_predictors) = names(X)
  
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
    
    choosed_predictors[index] = choosed_predictors[index] + 1
    
    # Updating ft
    ft = ft + v*g
    
    # Iterating
    m = m + 1
    
    g_optimum[[m]] <- g
    
  }
  
  # Matrix of all g vectors
  Matrix = do.call(cbind, g_optimum)
  
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
    
    final_vector = c(MSE, AICC)
    names(final_vector) = c("MSE", "IC") 
    
  } else {
    
    final_vector = c(MSE, AIC)
    names(final_vector) = c("MSE", "IC") 
    
  }
  
  #return(ft)
  #print(choosed_predictors)
  
  return(final_vector)
}

#Example
L2_boost(Y, X, 9, AICc = F)

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

# Function for minimizing AIC 
L2_min = function(Y_in = Y, X_in=X, ending = (length(Y)/2), AICc = F){
  
  IC = c()
  for(i in 1:ending){
    x = L2_boost(Y,X,i,AICc = AICc)[2]
    IC = c(IC, x)
  }
    min = unname(which.min(IC))
    
    #return(L2_boost(Y,X,M=min))
    return(min) # return the m* with the smallest AIC
}

# Example
L2_min(AICc = F)

# Final Example with selected m

L2_boost(Y,X,L2_min(),AICc = F)

