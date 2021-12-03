
library(forecast)
library(mboost)
library(quantreg)
library(scoringRules)

source("functions.R")

Y_or = read.csv("https://raw.githubusercontent.com/pedroskorin/L2_Boosting/master/l2_boosting/data/target.csv",
                               encoding = "UTF-8")[,4]
Y = read.csv("https://raw.githubusercontent.com/pedroskorin/L2_Boosting/master/l2_boosting/data/target.csv",
                                    encoding = "UTF-8")[-1,3]
X = read.csv("https://raw.githubusercontent.com/pedroskorin/L2_Boosting/master/l2_boosting/data/predictors.csv",
                               encoding = "UTF-8")[,-c(1,2)]

X_lag = add_lags(X,Y)

Y_lag = tail(Y,-11)
Y_or_lag = tail(Y_or,-11)
names = colnames(X_lag)
colnames(X_lag) = 1:(823*12)

# Example with v=0.1, h=1, M_max = 2500 and ratio_start = 80%

v_in = 0.1
h_in = 1
Mstop_in = 50
ratio_start_lag = 0.79

b_aic_1 = boosting_reg_aic(Y_or = Y_or_lag,Y = Y_lag, X_lag, v = v_in, h = 1, ratio_start = ratio_start_lag,
                       Mstop = Mstop_in)
b_aic_2 = boosting_reg_aic(Y_or = Y_or_lag,Y = Y_lag, X_lag, v = v_in, h = 2, ratio_start = ratio_start_lag,
                         Mstop = Mstop_in)
b_aic_3 = boosting_reg_aic(Y_or = Y_or_lag,Y = Y_lag, X_lag, v = v_in, h = 3, ratio_start = ratio_start_lag,
                         Mstop = Mstop_in)

b_kfold_50_1 = boosting_reg_kfold(Y_or_lag, Y_lag, X_lag, v = v_in, h = 1, ratio_start = ratio_start_lag,
                           Mstop = Mstop_in)
b_kfold_50_2 = boosting_reg_kfold(Y_or_lag, Y_lag, X_lag, v = v_in, h = 2, ratio_start = ratio_start_lag,
                             Mstop = Mstop_in)
b_kfold_50_3 = boosting_reg_kfold(Y_or_lag, Y_lag, X_lag, v = v_in, h = 3, ratio_start = ratio_start_lag,
                             Mstop = Mstop_in)

b_quantile_5_1 = boosting_reg_quantile(Y_or_lag, Y_lag, X_lag, v = v_in, h = 1, ratio_start = ratio_start_lag,
                           Mstop = 50, tau_in = 0.5)
b_quantile_5_2 = boosting_reg_quantile(Y_or_lag, Y_lag, X_lag, v = v_in, h = 2, ratio_start = ratio_start_lag,
                           Mstop = 50, tau_in = 0.5)
b_quantile_5_3 = boosting_reg_quantile(Y_or_lag, Y_lag, X_lag, v = v_in, h = 3, ratio_start = ratio_start_lag,
                           Mstop = 50, tau_in = 0.5)
b_quantile_95_3 = boosting_reg_quantile(Y_or_lag, Y_lag, X_lag, v = v_in, h = 3, ratio_start = ratio_start_lag,
                           Mstop = 50, tau_in = 0.95)
b_quantile_95_3 = boosting_reg_quantile(Y_or_lag, Y_lag, X_lag, v = v_in, h = 3, ratio_start = ratio_start_lag,
                          Mstop = 50, tau_in = 0.95)
b_quantile_95_3 = boosting_reg_quantile(Y_or_lag, Y_lag, X_lag, v = v_in, h = 3, ratio_start = ratio_start_lag,
                          Mstop = 50, tau_in = 0.95)
b_quantile_05_3 = boosting_reg_quantile(Y_or_lag, Y_lag, X_lag, v = v_in, h = 3, ratio_start = ratio_start_lag,
                          Mstop = 50, tau_in = 0.05)
b_quantile_05_3 = boosting_reg_quantile(Y_or_lag, Y_lag, X_lag, v = v_in, h = 3, ratio_start = ratio_start_lag,
                          Mstop = 50, tau_in = 0.05)
b_quantile_05_3 = boosting_reg_quantile(Y_or_lag, Y_lag, X_lag, v = v_in, h = 3, ratio_start = ratio_start_lag,
                          Mstop = 50, tau_in = 0.05)

ratio_start_in = 0.8

b_sarima_1 = SARIMA_bench(Y_or, Y, h=1, ratio_start = ratio_start_in)
b_sarima_2 = SARIMA_bench(Y_or, Y, h=2, ratio_start = ratio_start_in)
b_sarima_3 = SARIMA_bench(Y_or, Y, h=3, ratio_start = ratio_start_in)

b_qsarima_5_1 = quantile_sarima(Y_or, Y, h=1, ratio_start = 0.79, tau_in = 0.5)
b_qsarima_5_2 = quantile_sarima(Y_or, Y, h=2, ratio_start = 0.79, tau_in = 0.5)
b_qsarima_5_3 = quantile_sarima(Y_or, Y, h=3, ratio_start = 0.79, tau_in = 0.5)

b_qsarima_95_1 = quantile_sarima(Y_or, Y, h=1, ratio_start = 0.79, tau_in = 0.95)
b_qsarima_95_2 = quantile_sarima(Y_or, Y, h=2, ratio_start = 0.79, tau_in = 0.95)
b_qsarima_95_3 = quantile_sarima(Y_or, Y, h=3, ratio_start = 0.79, tau_in = 0.95)

b_qsarima_05_1 = quantile_sarima(Y_or, Y, h=1, ratio_start = 0.79, tau_in = 0.05)
b_qsarima_05_2 = quantile_sarima(Y_or, Y, h=2, ratio_start = 0.79, tau_in = 0.05)
b_qsarima_05_3 = quantile_sarima(Y_or, Y, h=3, ratio_start = 0.79, tau_in = 0.05)

n_tot_lag <- length(Y_lag)
n_out_lag <- ceiling(n_tot_lag - ratio_start_lag*n_tot_lag)
ind_out_lag <- seq(to = n_tot_lag, by = 1, length = n_out_lag)

n_tot <- length(Y)
n_out <- ceiling(n_tot - ratio_start_in*n_tot)
ind_out <- seq(to = n_tot, by = 1, length = n_out)

# PERFORMANCE MEASURES
b_aic_e_1 = evaluation(b_aic_1$forecast, Y_or_lag, ind_out_lag, "boost_aic h=1")
b_aic_e_2 = evaluation(b_aic_2$forecast, Y_or_lag, ind_out_lag, "boost_aic h=2")
b_aic_e_3 = evaluation(b_aic_3$forecast, Y_or_lag, ind_out_lag, "boost_aic h=3")

b_kfold_e_1 = evaluation(b_kfold_50_1$forecast, Y_or_lag, ind_out_lag, "boost_kfold h=1")
b_kfold_e_2 = evaluation(b_kfold_50_2$forecast, Y_or_lag, ind_out_lag, "boost_kfold h=2")
b_kfold_e_3 = evaluation(b_kfold_50_3$forecast, Y_or_lag, ind_out_lag, "boost_kfold h=3")

# RESULTS FOR BENCH

b_sarima_e_1 = evaluation(b_sarima_1$benchmark, Y_or, ind_out, "SARIMA h=1")
b_sarima_e_2 = evaluation(b_sarima_2$benchmark, Y_or, ind_out, "SARIMA h=2")
b_sarima_e_3 = evaluation(b_sarima_3$benchmark, Y_or, ind_out, "SARIMA h=3")

# dm test

dm.test(exp(b_aic_1$forecast)[-1]-exp(tail(Y_or,38)),
        exp(b_sarima_1$benchmark)[-1] - exp(tail(Y_or,38)),
        h=1, power = 1, alternative = "less")

dm.test(exp(b_aic_1$forecast)[-1]-exp(tail(Y_or,38)),
        exp(b_kfold_50_1$forecast)[-1] - exp(tail(Y_or,38)),
        h=1, power = 1, alternative = "greater")

# Variable Importance
## Varimp

var_imp = rowSums(b_aic_1$varimp)/ncol(b_aic_1$varimp)

names(var_imp) = c("Intercept", names)

df_1 = data.frame(imp = var_imp[order(var_imp, decreasing = T)],
                name = names(var_imp[order(var_imp, decreasing = T)]))

## Frequency

n_prev = length(b_aic_3$selected)
freq_search = function(x) {
  y=0
  for (i in 1:n_prev) {
    y = y + (x %in% b_aic_2$selected[[i]])
  }
  return(y/n_prev)
}

vector_freq = mapply(freq_search, 1:(ncol(X_lag)+1))
df_freq_2 = data.frame(freq = vector_freq[order(vector_freq, decreasing = T)],
                     index = c("Intercept", names)[order(vector_freq, decreasing = T)])
View(df_freq_2)

# Plotting graph

library(ggplot2)
library(reshape2)

# Plot
# h=1

date = seq(as.Date("2003/1/1"), as.Date("2017/10/1"), "months")

intervalo = 48

Boosting = exp(tail(Y_or, 48))
SARIMA = exp(tail(Y_or, 48))
Original = exp(tail(Y_or, 48))
date = tail(date, 48)

n = length(tail(exp(b_aic_1$forecast),-1))

Boosting[c(rep(F,length(date)-n),rep(T,n))] = tail(exp(b_aic_1$forecast),-1)
SARIMA[c(rep(F,length(date)-n),rep(T,n))] = tail(exp(b_sarima_1$benchmark),-1)

df = data.frame(date = date, SARIMA = SARIMA, 
                Boosting  = Boosting, Original = Original)

meltdf = melt(df, id = "date")
colnames(meltdf)[2] = "Model"

meltdf_1 = meltdf

# h=2

date = seq(as.Date("2003/2/1"), as.Date("2017/10/1"), "months")

intervalo = 48

Boosting = exp(tail(Y_or, 48))
SARIMA = exp(tail(Y_or, 48))
Original = exp(tail(Y_or, 48))
date = tail(date, 48)

n = length(tail(exp(b_aic_2$forecast),-1))

Boosting[c(rep(F,length(date)-n),rep(T,n))] = tail(exp(b_aic_2$forecast),-1)
SARIMA[c(rep(F,length(date)-n),rep(T,n))] = tail(exp(b_sarima_2$benchmark),-1)

df = data.frame(date = date, SARIMA = SARIMA, 
                Boosting  = Boosting, Original = Original)

meltdf = melt(df, id = "date")
colnames(meltdf)[2] = "Model"

meltdf_2 = meltdf

# h=3

date = seq(as.Date("2003/2/1"), as.Date("2017/10/1"), "months")

intervalo = 48

Boosting = exp(tail(Y_or, 48))
SARIMA = exp(tail(Y_or, 48))
Original = exp(tail(Y_or, 48))
date = tail(date, 48)

n = length(tail(exp(b_aic_3$forecast),-1))

Boosting[c(rep(F,length(date)-n),rep(T,n))] = tail(exp(b_aic_3$forecast),-1)
SARIMA[c(rep(F,length(date)-n),rep(T,n))] = tail(exp(b_sarima_3$benchmark),-1)

df = data.frame(date = date, SARIMA = SARIMA, 
                Boosting  = Boosting, Original = Original)

meltdf = melt(df, id = "date")
colnames(meltdf)[2] = "Model"

meltdf_3 = meltdf

meltdf_1$date = seq(as.Date("2013/11/1"), as.Date("2017/10/1"), "months")
meltdf_2$date = seq(as.Date("2013/11/1"), as.Date("2017/10/1"), "months")
meltdf_3$date = seq(as.Date("2013/11/1"), as.Date("2017/10/1"), "months")

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
