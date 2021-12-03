
##########################################
## Point Estimate Forecasting Functions ##
##########################################

# L2-Boosting AIC function
boosting_reg_aic <- function(Y_or, Y, X, v, h, ratio_start = 0.8, Mstop = 3500) {
  
  n_tot <- length(Y)
  n_out <- ceiling(n_tot - ratio_start*n_tot)
  ind_out <- seq(to = n_tot, by = 1, length = n_out)
  Y_predicted <- c(Y_or[ind_out[1]])
  varimp_df <- data.frame(rep(0, (ncol(X)+1)))
  selected_var <- c()
  
  for(i in 1:n_out){
    ind_in <- seq(from = 1, to = ind_out[i] - h, by = 1)
    y_extra <- c()
    x_reg <- X[head(ind_in,-1),] # x independent t = 1, ..., T.in-h
    x0_reg <- matrix(X[tail(ind_in,1),], nrow = 1)
    
    for(j in 1:h) {

      # expanding window
      y_dep <- append(Y[tail(ind_in,-j)], y_extra)
      y_reg <- as.matrix(y_dep)
      
      # finding m*
      model_1 <- glmboost(y_reg ~ ., data = x_reg,
                         family = Gaussian(),
                         control = boost_control(mstop = Mstop, nu = v),
                         center = T)
      
      AIC <- AIC(model_1, method = "corrected" , df = "actset")
      x0_reg_df <- data.frame(t(data.frame(unlist(x0_reg))))
      colnames(x0_reg_df) <- colnames(x_reg)
      y_predicted <- unname(predict(model_1[mstop(AIC)], newdata = x0_reg_df,
                                       type = "response")[1,1])
      cat("Selected M is: ", mstop(AIC), "\n")
      
      # visualizing selected predictors varimp
      varimp_df_partial <- data.frame(varimp(model_1))
      sum_reduction <- sum(varimp_df_partial[,1])
      varimp_partial <- varimp_df_partial[,1]/sum_reduction
      varimp_df <- cbind(varimp_df, varimp_partial)
      
      # visualizing selected predictors frequency
      selected_var <- append(selected_var, list(model_1$xselect()))
      
      # output
      y_extra <- append(y_extra, y_predicted)
      
    }
    
    Y_predicted <- append(Y_predicted, Y_or[(ind_out[1]+i-(h))] + sum(y_extra))
    print(i/n_out)
    
  }
  results <- list(forecast = Y_predicted,
                  varimp = varimp_df[,-1],
                  selected = selected_var
  )
  return(results)
}

# L2-Boosting K-fold function
boosting_reg_kfold <- function(Y_or, Y, X, v, h, ratio_start = 0.8, Mstop = 3500) {
  
  n_tot <- length(Y)
  n_out <- ceiling(n_tot - ratio_start*n_tot)
  ind_out <- seq(to = n_tot, by = 1, length = n_out)
  Y_predicted <- c(Y_or[ind_out[1]])
  varimp_df <- data.frame(rep(0, (ncol(X)+1)))
  selected_var <- c()
  
  for(i in 1:n_out){
    ind_in <- seq(from = 1, to = ind_out[i] - h, by = 1)
    y_extra <- c()
    x_reg <- X[head(ind_in,-1),] # x independent t = 1, ..., T.in-h
    x0_reg <- matrix(X[tail(ind_in,1),], nrow = 1)
    
    for(j in 1:h) {
      
      # expanding window
      y_dep <- append(Y[tail(ind_in,-j)], y_extra)
      y_reg <- as.matrix(y_dep)
      
      # finding m*
      model_1 <- glmboost(y_reg ~ ., data = x_reg,
                         family = Gaussian(),
                         control = boost_control(mstop = Mstop, nu = v),
                         center = T)
      cv10f <- cv(model.weights(model_1), type = "kfold")
      cvm <- cvrisk(model_1, folds = cv10f, papply = lapply)
      
      # AIC = AIC(model_1, method = "corrected" , df = "actset")
      x0_reg_df <- data.frame(t(data.frame(unlist(x0_reg))))
      colnames(x0_reg_df) <- colnames(x_reg)
      y_predicted <- unname(predict(model_1[mstop(cvm)], newdata = x0_reg_df,
                                       type = "response")[1,1])
      cat("Selected M is: ", mstop(cvm), "\n")
      
      # visualizing selected predictors varimp
      varimp_df_partial <- data.frame(varimp(model_1))
      sum_reduction <- sum(varimp_df_partial[,1])
      varimp_partial <- varimp_df_partial[,1]/sum_reduction
      varimp_df <- cbind(varimp_df, varimp_partial)
      
      # visualizing selected predictors frequency
      selected_var <- append(selected_var, list(model_1$xselect()))
      
      # output
      y_extra <- append(y_extra, y_predicted)
      
    }
    Y_predicted <- append(Y_predicted, Y_or[(ind_out[1]+i-(h))] + sum(y_extra))
    print(i/n_out)
    
  }
  results <- list(forecast = Y_predicted,
                  varimp = varimp_df[,-1],
                  selected = selected_var)
  return(results)
}

# SARIMA Function
SARIMA_bench <- function(Y_or, Y, h, ratio_start = 0.8) {
  
  n_tot <- length(Y)
  n_out <- ceiling(n_tot - ratio_start*n_tot)
  ind_out <- seq(to = n_tot, by = 1, length = n_out)
  Y_arima <- c(Y_or[ind_out[1]])
  
  for(i in 1:n_out){
    ind_in <- seq(from = 1, to = ind_out[i] - h, by = 1)
    bench <- arima(exp(Y_or[ 1:(ind_out[i] - h + 1) ]), c(1,1,0)
                  , seasonal = list(order = c(1,1,0), period = 12)
    )
    
    forecast_bench <- forecast(bench, h)
    y_predicted_bench <- forecast_bench$mean[h]
    y_predicted_arima <- log(y_predicted_bench)
    Y_arima <- append(Y_arima, (y_predicted_arima))
    print(i/n_out)
    
  }
  results <- list(benchmark = Y_arima)
  return(results)
}

#############################################
## Interval Estimate Forecasting Functions ##
#############################################

# L2-Boosting quantile function
boosting_reg_quantile <- function(Y_or, Y, X, v, h, ratio_start = 0.8, Mstop = 3500, tau_in = 0.5, offset_in = 0.5, m_mult = 4) {
  
  n_tot <- length(Y)
  n_out <- ceiling(n_tot - ratio_start*n_tot)
  ind_out <- seq(to = n_tot, by = 1, length = n_out)
  Y_predicted <- c(Y_or[ind_out[1]])
  varimp_df <- data.frame(rep(0, (ncol(X)+1)))
  selected_var <- c()
  
  for(i in 1:n_out){
    ind_in <- seq(from = 1, to = ind_out[i] - h, by = 1)
    y_extra <- c()
    x_reg <- X[head(ind_in,-1),] # x independent t = 1, ..., T.in-h
    x0_reg <- matrix(X[tail(ind_in,1),], nrow = 1)
    
    for(j in 1:h) {
      
      # expanding window
      y_dep <- append(Y[tail(ind_in,-j)], y_extra)
      y_reg <- as.matrix(y_dep)
    
      # finding m*
      model_1 <- glmboost(y_reg ~ ., data = x_reg,
                         family = QuantReg(tau=tau_in, qoffset = offset_in),
                         control = boost_control(mstop = m_mult*Mstop, nu = v),
                         center = T)
      
      model_2 <- glmboost(y_reg ~ ., data = x_reg,
                         family = Gaussian(),
                         control = boost_control(mstop = Mstop, nu = v),
                         center = T)
      
      AIC <- AIC(model_2, method = "corrected" , df = "actset")
      x0_reg_df <- data.frame(t(data.frame(unlist(x0_reg))))
      colnames(x0_reg_df) <- colnames(x_reg)
      y_predicted <- unname(predict(model_1[m_mult*mstop(AIC)], newdata = x0_reg_df
                                       #,type = "link"
      )[1, 1])
      cat("Selected M is: ", m_mult*mstop(AIC), "\n")

      # visualizing selected predictors varimp
      varimp_df_partial <- data.frame(varimp(model_1))
      sum_reduction <- sum(varimp_df_partial[,1])
      varimp_partial <- varimp_df_partial[,1]/sum_reduction
      varimp_df <- cbind(varimp_df, varimp_partial)
      
      # visualizing selected predictors frequency
      selected_var <- append(selected_var, list(model_1$xselect()))
      
      # output
      y_extra <- append(y_extra, y_predicted)
      
    }
    Y_predicted <- append(Y_predicted, Y_or[(ind_out[1]+i-(h))] + sum(y_extra))
    print(i/n_out)
    
  }
  results <- list(forecast = Y_predicted,
                  varimp = varimp_df[,-1],
                  selected = selected_var
  )
  return(results)
}

# Sarima quantile function
quantile_sarima <- function(Y_or, Y, h, ratio_start = 0.79, tau_in) {
  Y_lag <- tail(Y,-12)
  Y_1 <- head(tail(Y,-11),-1)
  Y_12 <- head(Y,-12)  
  X <- as.data.frame(cbind(Y_1,Y_12))
  
  n_tot <- length(Y_lag)
  n_out <- ceiling(n_tot - ratio_start*n_tot)
  Y_or_lag <- tail(Y_or,-12)
  ind_out <- seq(to = n_tot, by = 1, length = n_out)
  Y_predicted <- c(Y_or[ind_out[1]])
  
  for(i in 1:n_out){
    ind_in <- seq(from = 1, to = ind_out[i] - h, by = 1)
    y_extra <- c()
    x_reg <- X[head(ind_in,-1),] # x independent t = 1, ..., T.in-h
    x0_reg <- matrix(X[tail(ind_in,1),], nrow = 1)
    
    for(j in 1:h) {
      
      # expanding window
      y_dep <- append(Y[tail(ind_in,-j)], y_extra)
      y_reg <- as.matrix(y_dep)
      q <- rq(Y_lag ~ Y_1 + Y_12, tau = tau_in)
      x0_reg_df <- data.frame(t(data.frame(unlist(x0_reg))))
      colnames(x0_reg_df) <- colnames(x_reg)
      y_predicted <- unname(predict(q, newdata = x0_reg_df,
                                   type = "percentile"))
      # output
      y_extra <- append(y_extra, y_predicted)
    }
    Y_predicted <- append(Y_predicted, Y_or[(ind_out[1]+i-(h))] + sum(y_extra))
    print(i/n_out)
    
  }
  results <- list(forecast = Y_predicted)
  
  return(results)
}

#####################
## Other Functions ##
#####################

# Performance measures
evaluation <- function(Z, W, index, texto) {
  
  cat("Evaluation of", texto)
  
  MAPE = mean((abs(exp(W[index+1])[]-exp(Z[-1])[])/exp(W[index+1])[]))*100
  cat("\n MAPE:", MAPE)
  
  MPE = max((exp(Z[-1]) - exp(W[index+1]))/(exp(W[index+1])))*100
  MNE = min((exp(Z[-1]) - exp(W[index+1]))/(exp(W[index+1])))*100
  cat("\n MPE: ", MPE)
  cat("\n MNE: ", MNE)
  
  P90 = quantile(abs((exp(Z[-1]) - exp(W[index+1]))/(exp(W[index+1])))*100, 0.9)
  P95 = quantile(abs((exp(Z[-1]) - exp(W[index+1]))/(exp(W[index+1])))*100, 0.95)
  cat("\n P90: ", P90)
  cat("\n P95: ", P95)
  
  RMSFE = mean((exp(Z[-1]) - exp(W[index+1]))^2)
  cat("\n RMSFE: ", RMSFE)
  return(RMSFE)

}

# Adding Lags (There is definitely an easier way to do this)
add_lags <- function(X,Y) {
  X <- cbind(Y, X)
  
  X_sem_1 <- head(X,-1)
  X_sem_2 <- head(X_sem_1,-1)
  X_sem_3 <- head(X_sem_2,-1)
  X_sem_4 <- head(X_sem_3,-1)
  X_sem_5 <- head(X_sem_4,-1)
  X_sem_6 <- head(X_sem_5,-1)
  X_sem_7 <- head(X_sem_6,-1)
  X_sem_8 <- head(X_sem_7,-1)
  X_sem_9 <- head(X_sem_8,-1)
  X_sem_10 <- head(X_sem_9,-1)
  X_sem_11 <- head(X_sem_10,-1)
  
  X_s_0 <- tail(X,nrow(X_sem_11))
  X_s_1 <- tail(X_sem_1,nrow(X_sem_11))
  X_s_2 <- tail(X_sem_2,nrow(X_sem_11))
  X_s_3 <- tail(X_sem_3,nrow(X_sem_11))
  X_s_4 <- tail(X_sem_4,nrow(X_sem_11))
  X_s_5 <- tail(X_sem_5,nrow(X_sem_11))
  X_s_6 <- tail(X_sem_6,nrow(X_sem_11))
  X_s_7 <- tail(X_sem_7,nrow(X_sem_11))
  X_s_8 <- tail(X_sem_8,nrow(X_sem_11))
  X_s_9 <- tail(X_sem_9,nrow(X_sem_11))
  X_s_10 <- tail(X_sem_10,nrow(X_sem_11))
  X_s_11 <- tail(X_sem_11,nrow(X_sem_11))

  X_total <- cbind(X_s_0,
                  X_s_1,
                  X_s_2,
                  X_s_3,
                  X_s_4,
                  X_s_5,
                  X_s_6,
                  X_s_7,
                  X_s_8,
                  X_s_9,
                  X_s_10,
                  X_s_11
                  )
  
  troca_nome <- function(x, numero = n) {return(paste("L",n," - ", x, sep = ""))}
  nome_all <- c()
  for (i in 1:12) {
    n <- i
    nome_in <- unname(lapply(colnames(X), troca_nome))
    nome_all <- c(nome_all, nome_in)
  }
  colnames(X_total) <- unlist(nome_all)
  return(X_total)
}
