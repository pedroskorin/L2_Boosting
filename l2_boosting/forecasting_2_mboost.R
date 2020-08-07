
Y_or = log(consumo_energia_RS)
Y = diff(log(consumo_energia_RS))

Y_or = log(nacional_mensal[,70]*1000)[-1]
Y = diff(log(nacional_mensal[,70]*1000))[-1]

X = X[,-1]

Y_or = read.csv("https://raw.githubusercontent.com/pedroskorin/L2_Boosting/master/l2_boosting/dados/target.csv",
                               encoding = "UTF-8")[,4]
Y = read.csv("https://raw.githubusercontent.com/pedroskorin/L2_Boosting/master/l2_boosting/dados/target.csv",
                                    encoding = "UTF-8")[,3]
X = read.csv("https://raw.githubusercontent.com/pedroskorin/L2_Boosting/master/l2_boosting/dados/predictors.csv",
                               encoding = "UTF-8")[,-c(1,2)]

library(forecast)
library(mboost)

prediciton_boost_2 = function(Y_or, Y, X, v, h, ratio_start = 0.8, Mstop = 3500) {
  
  n_tot <- length(Y)
  n_out <- ceiling(n_tot - ratio_start*n_tot)
  
  ind_out <- seq(to = n_tot, by = 1, length = n_out)
  
  Y_predicted_uni = c(Y_or[ind_out[1]])
  Y_predicted_mes = c(Y_or[ind_out[1]])
  Y_arima = c(Y_or[ind_out[1]])
  
  varimp_df = data.frame(rep(0, (ncol(X)+2)))
  
  for(i in 1:n_out){
    
    ind_in <- seq(from = 1, to = ind_out[i] - h, by = 1)
    
    y_extra_uni = c()
    
    y_ind <- Y[head(ind_in,-1)] # y independente t = 1, ..., T.in-h
    x_ind <- X[head(ind_in,-1),] # x independente t = 1, ..., T.in-h
    x_reg <- cbind(y_ind,x_ind) # modelo auto-regressivo em y e defasado em x

    x0_reg <- matrix(c(Y[tail(ind_in,1)],X[tail(ind_in,1),]), nrow = 1)
    
    for(j in 1:h) {
      
      # expanding window
      
      y_dep_uni <- append(Y[tail(ind_in,-j)], y_extra_uni) # variavel dependente t = h+1, ..., T.in

      y_reg_uni <- as.matrix(y_dep_uni)
      
      # Encontrando M Ótimo
      
      model_1 = glmboost(y_reg_uni ~ ., data = x_reg,
                         family = Gaussian(),
                         control = boost_control(mstop = 2500, nu = v),
                         center = T)
      
      AIC = AIC(model_1, method = "corrected" , df = "actset")
      
      x0_reg_df = data.frame(t(data.frame(unlist(x0_reg))))
      colnames(x0_reg_df) = colnames(x_reg)
      
      y_predicted_uni = unname(predict(model_1[mstop(AIC)], newdata = x0_reg_df,
                                       type = "response")[1,1])
      
      cat("O M ÓTIMO É", mstop(AIC), "\n")
      
      # Visualização de variáveis escolhidas

      varimp_df_partial = data.frame(varimp(model_1))
      
      sum_reduction = sum(varimp_df_partial[,1])
      
      varimp_partial = varimp_df_partial[,1]/sum_reduction
      
      varimp_df = cbind(varimp_df, varimp_partial)
            
      # Feitoria
      
      y_extra_uni = append(y_extra_uni, y_predicted_uni)
      

    }
    
      bench = arima(exp(Y_or[ 1:(ind_out[i] - h + 1) ]), c(1,1,0)
                    , seasonal = list(order = c(1,1,0), period = 12)
      )

      #bench = auto.arima(ts(exp(Y_or[ 1:(ind_out[i] - h + 1) ]), frequency = 12))
      forecast_bench = forecast(bench, h)
      y_predicted_bench = forecast_bench$mean[h]
      
      y_predicted_arima = log(y_predicted_bench)
      
      Y_arima = append(Y_arima, (y_predicted_arima))
      
      Y_predicted_uni = append(Y_predicted_uni, Y_or[(ind_out[1]+i-(h))] + sum(y_extra_uni))
      
      print(i/n_out)
      
  }
  results <- list(forecast_uni = Y_predicted_uni,
                  benchmark = Y_arima,
                  varimp = varimp_df[,-1]
  #                mng = management, 
  #                FE_boost = Y_test - Y_predicted,
  #                FE_bench = Y_test - Y_arima,
  #                test = Y_test)
  )
  return(results)
}

v_in = 0.1
h_in = 2
Mstop_in = 2500
ratio_start_in = 0.805

b = prediciton_boost_2(Y_or,
                       Y,
                       X,
                       v = v_in,
                       h = h_in,
                       ratio_start = ratio_start_in,
                       Mstop = Mstop_in)

# Ver importancia total

var_imp = rowSums(b$varimp)/ncol(b$varimp)

names(var_imp) = c("Intercept", "Y_{t-1}", colnames(X))

df = data.frame(imp = var_imp[order(var_imp, decreasing = T)][1:5],
                name = names(var_imp[order(var_imp, decreasing = T)][1:5]))

library(ggplot2)
library(ggthemes)

ggplot(df, aes(x=name, y=imp, fill = name)) + 
  geom_bar(stat = "identity") +
  coord_flip() +
  theme_wsj()+ scale_colour_wsj("colors6")+
  ggtitle("Iris data")

# MAPE
MAPE_boost_uni = mean((abs(exp(Y_or[ind_out+1])[]-exp(b$forecast_uni[-1])[])/exp(Y_or[ind_out+1])[]))*100
MAPE_boost_uni

MPE_boost = max((exp(b$forecast_uni[-1]) - exp(Y_or[ind_out+1]))/(exp(Y_or[ind_out+1])))*100
MNE_boost = min((exp(b$forecast_uni[-1]) - exp(Y_or[ind_out+1]))/(exp(Y_or[ind_out+1])))*100
print(MPE_boost)
print(MNE_boost)

P90_boost = quantile(abs((exp(b$forecast_uni[-1]) - exp(Y_or[ind_out+1]))/(exp(Y_or[ind_out+1])))*100, 0.9)
P95_boost = quantile(abs((exp(b$forecast_uni[-1]) - exp(Y_or[ind_out+1]))/(exp(Y_or[ind_out+1])))*100, 0.95)
print(P90_boost)
print(P95_boost)

#MAPE_boost_mes = mean((abs(exp(Y_or[ind_out+1])-exp(b$forecast_mes[-1]))/exp(Y_or[ind_out+1])))*100
#MAPE_boost_mes

MAPE_bench = mean((abs(exp(Y_or[ind_out+1])[]-exp(b$benchmark[-1])[])/exp(Y_or[ind_out+1])[]))*100
MAPE_bench

MPE_bench = max((exp(b$benchmark[-1]) - exp(Y_or[ind_out+1]))/(exp(Y_or[ind_out+1])))*100
MNE_bench = min((exp(b$benchmark[-1]) - exp(Y_or[ind_out+1]))/(exp(Y_or[ind_out+1])))*100
print(MPE_bench)
print(MNE_bench)

P90_benchmark = quantile(abs((exp(b$benchmark[-1]) - exp(Y_or[ind_out+1]))/(exp(Y_or[ind_out+1])))*100, 0.9)
P95_benchmark = quantile(abs((exp(b$benchmark[-1]) - exp(Y_or[ind_out+1]))/(exp(Y_or[ind_out+1])))*100, 0.95)
print(P90_benchmark)
print(P95_benchmark)

#

plot(exp(Y_or[154:192]), type = "l", ylim = c(exp(14.2), exp(14.7)),
     main = paste("h=",h_in," | v=",v_in," | M_a =", round(MAPE_bench,2), "% | M_b_uni = ",round(MAPE_boost_uni,2),"%"))
# 
lines(exp(b$benchmark[]), col = "blue")

lines(exp(b$forecast_uni), col = "red")
