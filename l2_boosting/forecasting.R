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
    
    coefficient_control[index] = coefficient_control[index] + v*teta[index]
    
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
  
  n_tot <- length(Y)
  n_out <- ceiling(n_tot - ratio_start*n_tot)
  
  ind_out <- seq(to = n_tot, by = 1, length = n_out)
  
  Y_predicted = c()
  Y_arima = c()

  Y_test = Y[ind_out]
  
  management = list()
  
  #m_start = 3
  
  for(i in 1:n_out){
    
    #m_otimo = L2_min(Y_train, X_train, v, m_start)
    
    #m_start = ceiling((1/2)*m_otimo)
    
    ind_in <- seq(from = 1, to = ind_out[i] - h, by = 1) # expanding window
    
    y_dep <- Y[tail(ind_in,-h)] # variavel dependente t = h+1, ..., T.in
    y_ind <- Y[head(ind_in,-h)] # y independente t = 1, ..., T.in-h
    x_ind <- X[head(ind_in,-h),] # x independente t = 1, ..., T.in-h
    
    x_reg <- cbind(y_ind,x_ind) # modelo auto-regressivo em y e defasado em x
    y_reg <- as.matrix(y_dep)
    x0_reg <- matrix(c(Y[tail(ind_in,1)],X[tail(ind_in,1),]), nrow = 1)
    
    L2_predicted = Fitting(y_reg, x_reg, Mstop, v)  
    
    y_predicted = (as.vector(unlist(x0_reg)) %*% as.vector(L2_predicted[[3]][nrow(L2_predicted[[3]]),])) + mean(y_reg)
    
    Y_predicted = append(Y_predicted, y_predicted)
    
    bench = arima(Y[1:(ind_out[i]-h)], c(1,0,0)
                  , seasonal = list(order = c(1,0,0), period = 12)
    )
    
    #bench = auto.arima(Y[1:(ind_out[i]-h)], seasonal = T)
    print(bench$coef)
    print(bench$arma)
    forecast_bench = forecast(bench, h)
    y_predicted_bench = forecast_bench$mean[h]
    Y_arima = append(Y_arima, y_predicted_bench)
    
    management = append(management, list(L2_predicted))
    
    print(i/n_out)
    
  }
  results <- list(forecast = Y_predicted, benchmark = Y_arima, mng = management, 
                  FE_boost = Y_test - Y_predicted, FE_bench = Y_test - Y_arima, test = Y_test)
  return(results)
}

# Getting data from github without the first column
X = read.csv("https://raw.githubusercontent.com/pedroskorin/L2_Boosting/master/l2_boosting/dados/regressors_saz.csv", encoding = "UTF-8")[,-c(1,2)]
Y = read.csv("https://raw.githubusercontent.com/pedroskorin/L2_Boosting/master/l2_boosting/dados/objective_seasonal_adj.csv")[,2]

# Selecting data ####

as.vector(unlist(matrix(c(Y[tail(ind_in,3)][1],X[tail(ind_in,3),][1,]), nrow = 1))) %*% as.vector(L2_predicted[[3]][nrow(L2_predicted[[3]]),])

# Selecting Y
Y = diff(log(consumo_energia_RS))[-1]

# Selecting X

# Avaliacao específica ####

# Selecao de parametros
v = 0.5
h_in = 1
Mstop = 30
ratio_start = 0.8

# Forecasting
a = prediciton_boost(Y[],
                     X[,],
                     v, h_in, ratio_start, Mstop)

# Calculando Erros
RMSFE_boost = sqrt(sum((a$FE_boost)^2) *
                     (1/length(a$test)))
RMSFE_arima = sqrt(sum((a$FE_bench)^2) *
                     (1/length(a$test)))
rRMSFE = RMSFE_boost/RMSFE_arima
print(rRMSFE)

# Plotando grafico ####

plot((a$test[1:35]), type = "l", col = "black",
     ylim = c(1600000, max(a$forecast)*1.05),
     main = paste("v =", v, "| h =", h_in, "| M =", Mstop, "| rRMSFE =", round(rRMSFE,2)),
     ylab = "Y")

lines((a$benchmark)[1:35], col = "blue")
lines((a$forecast)[1:35], col = "red")

# Calculando MAPE
MAE_boost = mean(abs((exp(a$test)-exp(a$forecast))))
MAE_boost

MAPE_boost_exp = mean(abs((exp(a$test)-exp(a$forecast))/exp((a$test))))*100
MAPE_boost_exp

MAPE_boost = mean(abs(((a$test[1:35])-(a$forecast[1:35]))/(a$test[1:35])))*100
MAPE_boost

MAPE_bench_exp = mean(abs((exp(a$test)-exp(a$benchmark))/exp((a$test))))*100
MAPE_bench_exp

MAPE_bench = mean(abs(((a$test)-(a$benchmark))/(a$test)))*100
MAPE_bench 

MAPE_boost/MAPE_bench
MAPE_boost_exp/MAPE_bench_exp


MAE_bench = mean(abs((a$test-a$benchmark)))
MAPE_bench = mean(abs((a$test-a$benchmark)/(a$test)))*100

MAPE_nos = mean(abs((serie_normal)-serie_nivel_boost)/serie_normal)*100
print(MAPE_nos)

MAPE_eles = mean(abs((serie_normal)-serie_nivel_bench)/serie_normal)*100
print(MAPE_eles)

print(MAPE_nos/MAPE_eles)

print(MAE_bench)
print(MAPE_bench)

print(MAE_boost/MAE_bench)
print(MAPE_boost/MAPE_bench)

# Apresentacao de variaveis escolhidas

choosed = rep(0, ncol(cbind(Y,X)))
tamanho = Mstop

for(i in 1:(length(a$test)-1)) {

  w = which(a$mng[[i]][[3]][tamanho,] != 0)
      
  choosed[w] = choosed[w] + 1
        
} 

names(choosed) = colnames(cbind(Y,X))

var_choosed = sort(choosed[which(choosed != 0)]/length(a$test),
              decreasing = T)*100

View(var_choosed)

# Avaliacao entre diferentes valores de m e de v ####
# Selecao dos parametros
v = seq(0.1,1,0.1)
m = seq(1,150,10)
h_f = seq(1,12)

# Loop H

for (j in h_f) {

# Criacao das matrizes
  
results <- as.data.frame(matrix(ncol=length(m),nrow=0))
names(results) <- m

var <- as.data.frame(matrix(ncol=length(m),nrow=0))
names(var) <- m

# Loop

for (k in 1:length(m)){
  
#  par(mfrow = c(5,1), mar = c(1,2,1,2))
  
for (i in 1:length(v)){
  
  print(i/length(v)) 
  
  a = prediciton_boost(Y,
                       X,
                       v[i],
                       ratio_start = 0.8, h=j,
                       Mstop = m[k])
  
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

cat("Estamos em", round((k/length(m))*100 ,2), "% | dalers")
  
}

rownames(results) <- v
rownames(var) <- v

write.csv(results, paste("C:\\Users\\Pedro.Pedro-PC\\Documents\\GitHub\\results_h",as.character(j),".csv", sep = ""))

}
#write.csv(var, "C:\\Users\\Pedro.Pedro-PC\\Documents\\GitHub\\var3.csv")

# Retornando ao nível sem estacionaridade ####

serie_nivel_boost = c()

for (i in 1:length(a$FE_boost)) {
  
  j = length(a$FE_boost)-i+2
  nivel_orig = tail(consumo_energia_RS, j)[1]
  
  nivel_prev = nivel_orig * exp(a$FE_boost[i])
  
  serie_nivel_boost = append(serie_nivel_boost, nivel_prev)
  
}

serie_nivel_bench = c()

for (i in 1:length(a$FE_bench)) {
  
  j = length(a$FE_bench)-i+2
  nivel_orig = tail(consumo_energia_RS, j)[1]
  
  nivel_prev = nivel_orig * exp(a$FE_bench[i])
  
  serie_nivel_bench = append(serie_nivel_bench, nivel_prev)
  
}

serie_nivel_boost = serie_nivel_boost[-1]
serie_nivel_bench = serie_nivel_bench[-1]
serie_normal = tail(consumo_energia_RS,length(a$FE_boost))[-length(tail(consumo_energia_RS,length(a$FE_boost)))]

plot(serie_normal,type = "l", ylim = c(1500000, 2500000))

lines(serie_nivel_boost, col = "red")
lines(serie_nivel_bench, col = "blue")
lines(ser_n_bench_2, col = "green")
lines(ser_n_boost_2, col = "purple")

# Retornando ao nível sem est método 2

primeiro = tail(consumo_energia_RS, length(a$FE_boost)+1)[1]

edge = length(consumo_energia_RS) - length(a$FE_boost) + 1

ind_out = seq(edge, length(consumo_energia_RS))

for (i in 1:length(a$FE_boost)) {
  
  x_in = ser_n_boost_2[i]
  y_in = exp(a$FE_boost[i])*x_in
  
  ser_n_boost_2 = append(ser_n_boost_2, y_in)
  
}

ser_n_boost_2 = ser_n_boost_2[-1]

