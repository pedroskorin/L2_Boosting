#L2 BOOSTING ALGORITHM
#
#Forecasting using the L2 Boosting algorithm 
#
#theoretical source: "Jing Zang, 2014, Forecasting Aggregates with Disaggregate Variables: 
#Does boosting help to select the most informative predictors?"
#
#Authors: guislinden, pedroskorin
#UFRGS
#
#2019-6-19

# Getting funcitons

Fitting = function(Y, X, M, v){
  
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
  
  #for (j in M:1){
  
  #  H = as.matrix(X[,s[j]]) %*% (1/sum((X[,s[j]])^2)) %*% t(as.matrix(X[,s[j]]))
  
  #  if (j == M){
  
  #    C = diag(nrow(X)) - v * H
  
  #  } else {
  
  #    C = C %*% (diag(nrow(X)) - v * H)
  
  #  }
  
  #}
  
  ## Calculating Big Beta
  
  #B = diag(nrow(X)) - C
  
  ## Calculating sigma squared_squared
  
  #sigma_squared = nrow(X)^(-1) * sum((Y - ( B %*% Y ))^2)
  
  ## Calculating trace
  
  #df_m = sum(diag(B))
  
  ## Calculating AIC
  
  #AIC = log(sigma_squared) + (1 + df_m/nrow(X))/(1 - (df_m + 2)/nrow(X))
  
  ## Computing final vector
  
  final_vector = list(MSE,
                      1,
                      Matrix_coef, ft, Y-ft)
  return(final_vector)
}

L2_min = function(Y, X, v, m_start = 3){
  
  i = m_start
  x = Fitting(Y, X, i-1, v)[[2]]
  IC = c(Fitting(Y, X, i-2, v)[[2]], x)
  
  while(IC[length(IC)-1] > x) {
    
    x = Fitting(Y, X, i, v)[[2]]
    
    IC = c(IC, x)
    
    i = i + 1
    
  }
  
  return(i-2)
}

L2_boost = function(Y, X, v) {
  
  # Function for minimizing AIC 
  
  # Optimal M
  m_opt = L2_min(Y, X, v)
  
  # Fitting
  optimal = Fitting(Y, X, m_opt, v)
  
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

library(forecast)
prediciton_boost = function(Y, X, v, h, ratio_start = 0.75, Mstop = 100) {
  
  edge = seq(1, nrow(X), 3)[which.min(abs(seq(1, nrow(X), 3) - nrow(X)*ratio_start))]
  
  X_train = X[1:(edge-(h-1)),]
  Y_train = Y[1:(edge-(h-1))]
  
  Y_predicted = c()
  Y_arima = c()
  
  X_test = X[(edge+1):nrow(X),]
  Y_test = Y[(edge+1):nrow(X)]
  
  management = list()
  
  k = 0
  m_start = 3
  
  for(i in (edge+1):(nrow(X))){
    
    #m_otimo = L2_min(Y_train, X_train, v, m_start)
    
    #m_start = ceiling((1/2)*m_otimo)
    
    L2_predicted = Fitting(Y_train, X_train, Mstop, v)  
    
    y_predicted = as.matrix(X[i,]) %*% as.matrix(L2_predicted[[3]][nrow(L2_predicted[[3]]),])
    
    Y_predicted = append(Y_predicted, y_predicted)
    
    bench = arima(Y_train, c(1,0,0)
                  #, seasonal = list(order = c(1,0,0), period = 12)
    )
    
    #bench = auto.arima(Y_train)
    print(bench$coef)
    print(bench$arma)
    forecast_bench = forecast(bench, h)
    y_predicted_bench = forecast_bench$mean[h]
    Y_arima = append(Y_arima, y_predicted_bench)
    
    management = append(management, list(L2_predicted))
    
    k = k + 1  
    
    X_train = X[1:(edge-(h-1)+k),]
    Y_train = Y[1:(edge-(h-1)+k)]  
    print(k/length(Y_test))
    
  }
  results <- list(forecast = Y_predicted, benchmark = Y_arima, mng = management, 
                  FE_boost = Y_test - Y_predicted, FE_bench = Y_test - Y_arima, test = Y_test)
  return(results)
}

# Getting data from github without the first column
X = read.csv("https://raw.githubusercontent.com/pedroskorin/L2_Boosting/master/l2_boosting/dados/regressors.csv")[-c(1,2,3),-c(1,2)]
Y = read.csv("https://raw.githubusercontent.com/pedroskorin/L2_Boosting/master/l2_boosting/dados/objective_seasonal_adj.csv")[,2]

# Selecting data ####

# Selecting Y
Y = X[,]

# Selecting X

X_en = X[,-c(65,66,67,68,69,70,71,72,73,74,75,76,77,78,
             660,
             719,720,721,722,723,724,725,726,727,728,729,
             730,731,732,733,734,735,736,737,738,739)]

#X_en = X  

X_en = X_en[,!(colnames(X_en) %in% colnames(metereologicos))]



# Avaliacao específica ####

# Selecao de parametros
v = 0.5
h = 1
Mstop = 25
ratio_start = 0.75

# Forecasting
a = prediciton_boost(Y_en_rs,X_en,v, h, ratio_start, Mstop)

# Calculando Erros
RMSFE_boost = sqrt(sum((a$FE_boost)^2) *
                     (1/length(a$test)))
RMSFE_arima = sqrt(sum((a$FE_bench)^2) *
                     (1/length(a$test)))
rRMSFE = RMSFE_boost/RMSFE_arima
print(rRMSFE)

# Plotando grafico
plot(a$test, type = "l", col = "red", ylim = c(-0.1,0.1)
     , main = paste("v =", v, "| h =", h, "| M =", Mstop, "| rRMSFE =", round(rRMSFE,2)),
     ylab = "Y")

lines(a$benchmark, col = "blue")
lines(a$forecast, col = "green")

# Apresentacao de variaveis escolhidas

choosed = rep(0, ncol(X_en))
tamanho = Mstop

for(i in 1:length(a$test)) {

  w = which(a$mng[[i]][[3]][tamanho,] != 0)
      
  choosed[w] = choosed[w] + 1
        
} 

names(choosed) = colnames(X_en)

var_choosed = sort(choosed[which(choosed != 0)]/length(a$test),
              decreasing = T)*100

View(var_choosed)

# Avaliacao entre diferentes valores de m e de v ####
# Selecao dos parametros
v = c(0.1,0.5,0.7,1)
m = c(1,5,15,25)
h_f = 1

# Criacao das matrizes

results <- as.data.frame(matrix(ncol=length(m),nrow=0))
names(results) <- m

var <- as.data.frame(matrix(ncol=length(m),nrow=0))
names(var) <- m

# Loop

for (k in 1:length(m)){

cat("Estamos em", round((k/length(m))*100 ,2), "% | dalers")
  
#  par(mfrow = c(5,1), mar = c(1,2,1,2))
  
for (i in 1:length(v)){
  
  print(i/length(v)) 
  
  a = prediciton_boost(Y_en,X_en,v[i], h=h_f, Mstop = m[k])
  
#  plot(a$test, type = "l", col = "red",ylim = c(-.3,.3), main = paste("h =", h[i]))
  
#  lines(a$benchmark, col = "blue")
#  lines(a$forecast, col = "green")
  
  RMSFE_boost = sqrt(sum((a$FE_boost)^2) *
                       (1/length(a$test)))
  RMSFE_arima = sqrt(sum((a$FE_bench)^2) *
                       (1/length(a$test)))
  rRMSFE = RMSFE_boost/RMSFE_arima
  results[i,k] = rRMSFE
  
  var[i,k] = var(a$FE_boost)
#if(i ==1){
#  mtext(paste("mstop =", m[k]), outer = F, cex = 1, side = 4, adj = 1)
#
#}  
   
}
  
}

rownames(results) <- v
rownames(var) <- v

#write.csv(results, "C:\\Users\\Pedro.Pedro-PC\\Documents\\GitHub\\results.csv")
#write.csv(var, "C:\\Users\\Pedro.Pedro-PC\\Documents\\GitHub\\var3.csv")

