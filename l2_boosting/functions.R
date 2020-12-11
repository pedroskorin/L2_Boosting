# L2-Boosting Function

prediciton_boost_2 = function(Y_or, Y, X, v, h, ratio_start = 0.8, Mstop = 3500) {
  
  n_tot <- length(Y)
  n_out <- ceiling(n_tot - ratio_start*n_tot)
  
  ind_out <- seq(to = n_tot, by = 1, length = n_out)
  
  Y_predicted_uni = c(Y_or[ind_out[1]])
  Y_predicted_mes = c(Y_or[ind_out[1]])
  Y_arima = c(Y_or[ind_out[1]])
  
  varimp_df = data.frame(rep(0, (ncol(X)+1)))
  selected_var = c()
  
  for(i in 1:n_out){
    
    ind_in <- seq(from = 1, to = ind_out[i] - h, by = 1)
    
    y_extra_uni = c()
    
    x_reg <- X[head(ind_in,-1),] # x independent t = 1, ..., T.in-h
    
    x0_reg <- matrix(X[tail(ind_in,1),], nrow = 1)
    
    for(j in 1:h) {
      
      # expanding window
      
      y_dep_uni <- append(Y[tail(ind_in,-j)], y_extra_uni)
      
      y_reg_uni <- as.matrix(y_dep_uni)
      
      # finding m*
      
      model_1 = glmboost(y_reg_uni ~ ., data = x_reg,
                         family = Gaussian(),
                         control = boost_control(mstop = Mstop, nu = v),
                         center = T)
      
      AIC = AIC(model_1, method = "corrected" , df = "actset")
      
      x0_reg_df = data.frame(t(data.frame(unlist(x0_reg))))
      colnames(x0_reg_df) = colnames(x_reg)
      
      y_predicted_uni = unname(predict(model_1[mstop(AIC)], newdata = x0_reg_df,
                                       type = "response")[1,1])
      
      cat("Selected M is: ", mstop(AIC), "\n")
      
      # visualizing selected predictors varimp
      
      varimp_df_partial = data.frame(varimp(model_1))
      sum_reduction = sum(varimp_df_partial[,1])
      varimp_partial = varimp_df_partial[,1]/sum_reduction
      varimp_df = cbind(varimp_df, varimp_partial)
      
      # visualizing selected predictors frequency
      
      selected_var = append(selected_var, list(model_1$xselect()))
      
      # output
      
      y_extra_uni = append(y_extra_uni, y_predicted_uni)
      
      
    }
    
    Y_predicted_uni = append(Y_predicted_uni, Y_or[(ind_out[1]+i-(h))] + sum(y_extra_uni))
    
    print(i/n_out)
    
  }
  results <- list(forecast_uni = Y_predicted_uni,
                  varimp = varimp_df[,-1],
                  selected = selected_var
                  
  )
  return(results)
}

# SARIMA

SARIMA_bench = function(Y_or, Y, h, ratio_start = 0.8) {
  
  n_tot <- length(Y)
  n_out <- ceiling(n_tot - ratio_start*n_tot)
  
  ind_out <- seq(to = n_tot, by = 1, length = n_out)

  Y_arima = c(Y_or[ind_out[1]])
  
  for(i in 1:n_out){
    
    ind_in <- seq(from = 1, to = ind_out[i] - h, by = 1)
    
    bench = arima(exp(Y_or[ 1:(ind_out[i] - h + 1) ]), c(1,1,0)
                  , seasonal = list(order = c(1,1,0), period = 12)
    )
    
    forecast_bench = forecast(bench, h)
    y_predicted_bench = forecast_bench$mean[h]
    
    y_predicted_arima = log(y_predicted_bench)
    
    Y_arima = append(Y_arima, (y_predicted_arima))
    
    print(i/n_out)
    
  }
  results <- list(benchmark = Y_arima)
  return(results)
}

# Análise

evaluation = function(Z, W, index, texto) {
  
  cat("Evaluation of", texto)
  
  MAPE = mean((abs(exp(W[index+1])[]-exp(Z[-1])[])/exp(W[index+1])[]))*100
  cat("\n MAPE:", MAPE)
  
  MPE = max((exp(Z[-1]) - exp(W[index+1]))/(exp(W[index+1])))*100
  MNE = min((exp(Z[-1]) - exp(W[index+1]))/(exp(W[index+1])))*100
  cat("\n MPE: ", MPE)
  paste("\n MNE: ", MNE)
  
  P90 = quantile(abs((exp(Z[-1]) - exp(W[index+1]))/(exp(W[index+1])))*100, 0.9)
  P95 = quantile(abs((exp(Z[-1]) - exp(W[index+1]))/(exp(W[index+1])))*100, 0.95)
  cat("\n P90: ", P90)
  cat("\n P95: ", P95)
  
  RMSFE = mean((exp(Z[-1]) - exp(W[index+1]))^2)
  cat("\n RMFE: ", RMSFE)

}
