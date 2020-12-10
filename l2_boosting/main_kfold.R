
Y_or = read.csv("https://raw.githubusercontent.com/pedroskorin/L2_Boosting/master/l2_boosting/data/target.csv",
                encoding = "UTF-8")[,4]
Y = read.csv("https://raw.githubusercontent.com/pedroskorin/L2_Boosting/master/l2_boosting/data/target.csv",
             encoding = "UTF-8")[-1,3]
X = read.csv("https://raw.githubusercontent.com/pedroskorin/L2_Boosting/master/l2_boosting/data/predictors.csv",
             encoding = "UTF-8")[,-c(1,2)]

library(forecast)
library(mboost)

prediciton_boost = function(Y_or, Y, X, v, h, ratio_start = 0.8, Mstop = 3500) {
  
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
    
    y_ind <- Y[head(ind_in,-1)] # y independent t = 1, ..., T.in-h
    x_ind <- X[head(ind_in,-1),] # x independent t = 1, ..., T.in-h
    x_reg <- cbind(y_ind,x_ind)
    
    x0_reg <- matrix(c(Y[tail(ind_in,1)],X[tail(ind_in,1),]), nrow = 1)
    
    for(j in 1:h) {
      
      # expanding window
      
      y_dep_uni <- append(Y[tail(ind_in,-j)], y_extra_uni)
      
      y_reg_uni <- as.matrix(y_dep_uni)
      
      # finding m*
      
      model_1 = glmboost(y_reg_uni ~ ., data = x_reg,
                         family = Gaussian(),
                         control = boost_control(mstop = Mstop, nu = v),
                         center = T)
      cv10f <- cv(model.weights(model_1), type = "kfold")
      cvm <- cvrisk(model_1, folds = cv10f, papply = lapply)

     # AIC = AIC(model_1, method = "corrected" , df = "actset")
      
      x0_reg_df = data.frame(t(data.frame(unlist(x0_reg))))
      colnames(x0_reg_df) = colnames(x_reg)
      
      y_predicted_uni = unname(predict(model_1[mstop(cvm)], newdata = x0_reg_df,
                                       type = "response")[1,1])
      
      cat("O M ÓTIMO É", mstop(cvm), "\n")
      
      # visualizing selected predictors
      
      varimp_df_partial = data.frame(varimp(model_1))
      
      sum_reduction = sum(varimp_df_partial[,1])
      
      varimp_partial = varimp_df_partial[,1]/sum_reduction
      
      varimp_df = cbind(varimp_df, varimp_partial)
      
      # output
      
      y_extra_uni = append(y_extra_uni, y_predicted_uni)
      
      
    }
    
    bench = arima(exp(Y_or[ 1:(ind_out[i] - h + 1) ]), c(1,1,0)
                  , seasonal = list(order = c(1,1,0), period = 12)
    )
    
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
                  
  )
  return(results)
}

# Example with v=0.1, h=1, M_max = 2500 and ratio_start = 80%

v_in = 0.1
h_in = 1
Mstop_in = 100
ratio_start_in = 0.79

b = prediciton_boost(Y_or,
                     Y,
                     X,
                     v = v_in,
                     h = h_in,
                     ratio_start = ratio_start_in,
                     Mstop = Mstop_in)

n_tot <- length(Y)
n_out <- ceiling(n_tot - ratio_start_in*n_tot)

ind_out <- seq(to = n_tot, by = 1, length = n_out)

# RESULTS FOR BOOSTING

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

RMSFE_boost = mean((exp(b$forecast_uni[-1]) - exp(Y_or[ind_out+1]))^2)
RMSFE_boost

# RESULTS FOR BENCH

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

RMSFE_benchmark = mean((exp(b$benchmark[-1]) - exp(Y_or[ind_out+1]))^2)
RMSFE_benchmark

rRMSFE = RMSFE_boost/RMSFE_benchmark
rRMSFE

# Variable Importance

var_imp = rowSums(b$varimp)/ncol(b$varimp)

names(var_imp) = c("Intercept", "Y_{t-1}", colnames(X))

df = data.frame(imp = var_imp[order(var_imp, decreasing = T)],
                name = names(var_imp[order(var_imp, decreasing = T)]))

# Plotting graph

library(ggplot2)
library(reshape2)

# Plot

date = seq(as.Date("2002/1/1"), as.Date("2017/10/1"), "months")

intervalo = 135:190

Boosting = exp(Y_or[intervalo])
SARIMA = exp(Y_or[intervalo])
Original = exp(Y_or[intervalo])
date = date[intervalo]

n = length(tail(exp(b$forecast_uni),-1))

Boosting[c(rep(F,length(date)-n),rep(T,n))] = tail(exp(b$forecast_uni),-1)
SARIMA[c(rep(F,length(date)-n),rep(T,n))] = tail(exp(b$benchmark),-1)

df = data.frame(date = date, SARIMA = SARIMA, 
                Boosting  = Boosting, Original = Original)

meltdf = melt(df, id = "date")
colnames(meltdf)[2] = "Model"

meltdf_1$date = seq(as.Date("2013/3/1"), as.Date("2017/10/1"), "months")
meltdf_2$date = seq(as.Date("2013/3/1"), as.Date("2017/10/1"), "months")
meltdf_3$date = seq(as.Date("2013/3/1"), as.Date("2017/10/1"), "months")

# Multiple line plot
p1 = ggplot(meltdf_1, aes(x = date, y = value)) + 
  geom_line(aes(color = Model), size = 0.3) +
  scale_color_manual(values = c("#00AFBB", "#FC4E07", "black")) +
  theme_minimal() +
  xlab("") +
  ggtitle("h = 1") +
  ylab("")

p2 = ggplot(meltdf_2, aes(x = date, y = value)) + 
  geom_line(aes(color = Model), size = 0.3) +
  scale_color_manual(values = c("#00AFBB", "#FC4E07", "black")) +
  theme_minimal() +
  xlab("") +
  ylab("Consumption in MWh") +
  ggtitle("h = 2")

p3 = ggplot(meltdf_3, aes(x = date, y = value)) + 
  geom_line(aes(color = Model), size = 0.3) +
  scale_color_manual(values = c("#00AFBB", "#FC4E07", "black")) +
  theme_minimal() +
  xlab("Time") +
  ylab("") +
  ggtitle("h = 3")

library(cowplot)

prow <- plot_grid( p1 + theme(legend.position="none"),
                   p2 + theme(legend.position="none"),
                   p3 + theme(legend.position="none"),
                   align = 'vh',
                   labels = c("", "", ""),
                   hjust = -1,
                   nrow = 3
)

# extract the legend from one of the plots
# (clearly the whole thing only makes sense if all plots
# have the same legend, so we can arbitrarily pick one.)
legend <- get_legend(p1)

# add the legend to the row we made earlier. Give it one-third of the width
# of one plot (via rel_widths).
p <- plot_grid( prow, legend, rel_widths = c(3, .3))
p