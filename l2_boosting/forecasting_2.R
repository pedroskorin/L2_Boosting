
Y_or = log(consumo_energia_RS)[-1]
Y = diff(log(consumo_energia_RS))[-1]


library(forecast)

prediciton_boost_2 = function(Y_or, Y, X, v, h, ratio_start = 0.75, Mstop = 100) {
  
  n_tot <- length(Y)
  n_out <- ceiling(n_tot - ratio_start*n_tot)
  
  ind_out <- seq(to = n_tot, by = 1, length = n_out)
  
  Y_predicted = c(Y_or[ind_out[1]])
  Y_arima = c(Y_or[ind_out[1]])
  
  for(i in 1:n_out){
    
    ind_in <- seq(from = 1, to = ind_out[i] - h, by = 1)
    
    y_extra = c()
    
    y_ind <- Y[head(ind_in,-1)] # y independente t = 1, ..., T.in-h
    x_ind <- X[head(ind_in,-1),] # x independente t = 1, ..., T.in-h
    x_reg <- cbind(y_ind,x_ind) # modelo auto-regressivo em y e defasado em x

    x0_reg <- matrix(c(Y[tail(ind_in,1)],X[tail(ind_in,1),]), nrow = 1)
    
    for(j in 1:h) {
      
      # expanding window
      
      y_dep <- append(Y[tail(ind_in,-j)], y_extra) # variavel dependente t = h+1, ..., T.in

      y_reg <- as.matrix(y_dep)     
      
      L2_predicted = Fitting(y_reg, x_reg, Mstop, v)  
      
      y_predicted = as.vector(unlist(x0_reg)) %*% as.vector(L2_predicted[[3]][nrow(L2_predicted[[3]]),]) + mean(y_reg)
      
      y_extra = append(y_extra, y_predicted)
      

    }
    
      bench = arima(exp(Y_or[ 1:(ind_out[i] - h + 1) ]), c(1,1,0)
                    , seasonal = list(order = c(1,0,0), period = 12)
      )
      
      #bench = auto.arima(Y[1:(ind_out[i]-h)], seasonal = T)
      forecast_bench = forecast(bench, h)
      y_predicted_bench = forecast_bench$mean[h]
      
      y_predicted_arima = log(y_predicted_bench)
      
      Y_arima = append(Y_arima, (y_predicted_arima))
      
      Y_predicted = append(Y_predicted, Y_or[(ind_out[1]+i-(h))] + sum(y_extra))

      print(i/n_out)
      
  }
  results <- list(forecast = Y_predicted,
                  benchmark = Y_arima
  #                mng = management, 
  #                FE_boost = Y_test - Y_predicted,
  #                FE_bench = Y_test - Y_arima,
  #                test = Y_test)
  )
  return(results)
}

b = prediciton_boost_2(Y_or,
                       Y,
                       X,
                       v = 0.2,
                       h = 4,
                       ratio_start = 0.8,
                       Mstop = 70)


plot(exp(Y_or[153:191]), type = "l", ylim = c(exp(14.2), exp(14.8)))
# 
lines(exp(b$benchmark[]), col = "blue")

lines(exp(b$forecast), col = "red")

# MAPE
MAPE_boost_2 = mean((abs(exp(Y_or[ind_out+1])-exp(b$forecast[-1]))/exp(Y_or[ind_out+1])))*100
MAPE_boost_2

MAPE_bench_2 = mean((abs(exp(Y_or[ind_out+1])-exp(b$benchmark[-1]))/exp(Y_or[ind_out+1])))*100
MAPE_bench_2
