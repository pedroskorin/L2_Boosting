######################################################################
## Developed by: Pedro Skorin and Guilherme Lindenmeyer             ##
## Find more: https://github.com/pedroskorin/L2_Boosting            ##
## Paper:                                                           ##
## Using boosting for forecasting electric energy consumption       ##
## during a recession: a case study for the Brazilian State Rio     ##
## Grande do Sul                                                    ##
######################################################################

rm(list=ls())
require(mboost)
require(R.matlab)
require(forecast)
require(sandwich)
require(readxl)
source("l2_boosting/helper_functions.R")

##############
## Data Imp ##
##############

Y_or <- read.csv("https://raw.githubusercontent.com/pedroskorin/L2_Boosting/master/l2_boosting/data/target.csv", encoding = "UTF-8")[,4]
Y <- read.csv("https://raw.githubusercontent.com/pedroskorin/L2_Boosting/master/l2_boosting/data/target.csv", encoding = "UTF-8")[-1,3]
X <- read.csv("https://raw.githubusercontent.com/pedroskorin/L2_Boosting/master/l2_boosting/data/predictors.csv", encoding = "UTF-8")[,-c(1,2)]

##############
## Settings ##
##############

save_directory <- "objects/"    # Add / to end or leave empty!

# Boosting parameters
v_in <- 0.1               # Def=12    #
h_set <- 1:3              # Def=1:3   #
Mstop_in <- 50            # Def=50    #
ratio_start_in <- 0.8     # Def=0.8   #
ratio_start_lag <- 0.79   # Def=0.79  #
store_objects <- 0        # Def=0     #

###################
## 2 Forecasting ##
###################

# Data Processing
X_lag <- add_lags(X,Y)
Y_lag <- tail(Y,-11)
Y_or_lag <- tail(Y_or,-11)
names <- colnames(X_lag)
colnames(X_lag) <- 1:(823*12)

# Ensure reproducibility
set.seed(12345)

for (h_in in h_set) {
  
  # Point Estimate  
  # Fit L2-Boosting AIC
  b_aic <- boosting_reg_aic(Y_or_lag, Y_lag, X_lag, v_in, h_in, ratio_start_lag, Mstop_in)
  # Fit L2-Boosting kfold
  b_kfold <- boosting_reg_kfold(Y_or_lag, Y_lag, X_lag, v_in, h_in, ratio_start_lag, Mstop_in)
  # Fit Sarima Benchmark
  b_sarima <- SARIMA_bench(Y_or, Y, h_in, ratio_start_in)

  # Interval Estimate
  # Fit Quantile Boosting
  b_quantile_50 <- boosting_reg_quantile(Y_or_lag, Y_lag, X_lag, v_in, h_in, ratio_start_lag, Mstop_in, tau_in = 0.5)
  b_quantile_95 <- boosting_reg_quantile(Y_or_lag, Y_lag, X_lag, v_in, h_in, ratio_start_lag, Mstop_in, tau_in = 0.95)
  b_quantile_05 <- boosting_reg_quantile(Y_or_lag, Y_lag, X_lag, v_in, h_in, ratio_start_lag, Mstop_in, tau_in = 0.05)

  # Interval Estimate
  # Fit Quantile Sarima
  b_qsarima_50 <- quantile_sarima(Y_or, Y, h_in, ratio_start_in, tau_in = 0.5)
  b_qsarima_95 <- quantile_sarima(Y_or, Y, h_in, ratio_start_in, tau_in = 0.95)
  b_qsarima_05 <- quantile_sarima(Y_or, Y, h_in, ratio_start_in, tau_in = 0.05)
  
  # Storing Objects
  if (store_objects) {
    saveRDS(b_aic, file = paste(save_directory, "b_aic_", as.character(h_in), sep=""))
    saveRDS(b_kfold, file = paste(save_directory, "b_kfold_", as.character(h_in), sep=""))
    saveRDS(b_sarima, file = paste(save_directory, "b_sarima_", as.character(h_in), sep=""))
    saveRDS(b_quantile_50, file = paste(save_directory, "b_aic_", as.character(h_in), sep=""))
    saveRDS(b_quantile_95, file = paste(save_directory, "b_aic_", as.character(h_in), sep=""))
    saveRDS(b_quantile_05, file = paste(save_directory, "b_aic_", as.character(h_in), sep=""))
    saveRDS(b_qsarima_50, file = paste(save_directory, "b_aic_", as.character(h_in), sep=""))
    saveRDS(b_qsarima_95, file = paste(save_directory, "b_aic_", as.character(h_in), sep=""))
    saveRDS(b_qsarima_05, file = paste(save_directory, "b_aic_", as.character(h_in), sep=""))
   }
}

##############################
## 3.1 Performance measures ##
##############################

n_tot_lag <- length(Y_lag)
n_out_lag <- ceiling(n_tot_lag - ratio_start_lag*n_tot_lag)
ind_out_lag <- seq(to = n_tot_lag, by = 1, length = n_out_lag)

n_tot <- length(Y)
n_out <- ceiling(n_tot - ratio_start_in*n_tot)
ind_out <- seq(to = n_tot, by = 1, length = n_out)

# Import Objects created in step 2

# Print MAPE, P95, P90, MPE, MNE, change the models for other results
b_aic_e_1 <- evaluation(b_aic_1$forecast, Y_or_lag, ind_out_lag, "boost_aic h=1")
b_aic_e_2 <- evaluation(b_aic_2$forecast, Y_or_lag, ind_out_lag, "boost_aic h=2")
b_aic_e_3 <- evaluation(b_aic_3$forecast, Y_or_lag, ind_out_lag, "boost_aic h=3")

b_sarima_e_1 <- evaluation(b_sarima_1$benchmark, Y_or, ind_out, "SARIMA h=1")
b_sarima_e_2 <- evaluation(b_sarima_2$benchmark, Y_or, ind_out, "SARIMA h=2")
b_sarima_e_3 <- evaluation(b_sarima_3$benchmark, Y_or, ind_out, "SARIMA h=3")

# Diebold–Mariano Test

dm.test(exp(b_aic_1$forecast)[-1]-exp(tail(Y_or,38)),
        exp(b_sarima_1$benchmark)[-1] - exp(tail(Y_or,38)),
        h=1, power = 1, alternative = "less")

dm.test(exp(b_aic_1$forecast)[-1]-exp(tail(Y_or,38)),
        exp(b_kfold_50_1$forecast)[-1] - exp(tail(Y_or,38)),
        h=1, power = 1, alternative = "greater")

###########################################
## 3.2 Variable importance and frequency ##
###########################################

n_prev <- length(b_aic_3$selected)

vector_freq <- mapply(freq_search, 1:(ncol(X_lag)+1))
df_freq_2 <- data.frame(freq = vector_freq[order(vector_freq, decreasing = T)],
                       index = c("Intercept", names)[order(vector_freq, decreasing = T)])
View(df_freq_2)

#####################
## Plotting Graphs ##
#####################

# Import required libraries
require(ggplot2)
require(dplyr)
require(hrbrthemes)
require(grid)
loadfonts()

# PLOT 1 - FIGURE 2

# Data Imp
Y_test <- read.csv("https://raw.githubusercontent.com/pedroskorin/L2_Boosting/master/l2_boosting/data/target.csv",
                  encoding = "UTF-8")
Y_test$exp <- exp(Y_test$Y_or)

# Data Processing
desemprego <- read.csv("https://raw.githubusercontent.com/pedroskorin/L2_Boosting/master/l2_boosting/data/history/ipeadata[27-12-2020-01-20].csv")
colnames(desemprego) <- c("data", "des", "nada")
des_2 <- desemprego$des[206:395]
date <- seq(as.Date("2002/1/1"), as.Date("2017/10/1"), "months")
data <- data.frame("Time"=date,"des"=des_2,"luz"=Y_test$exp)

# Value used to transform the data
coeff <- 100000

# A few constants
desColor <- "#69b3a2"
luzColor <- rgb(0.2, 0.6, 0.9, 1)

p_1 = ggplot(data, aes(x=Time)) +
  geom_line( aes(y=des), size=1, color=desColor) + 
  geom_line( aes(y=luz/coeff), size=1, color=luzColor) +
  scale_y_continuous(
    # Features of the first axis
    name = "Unemployment Rate (%)",
    # Add a second axis and specify its features
    sec.axis = sec_axis(~.*coeff, name="Elecricity Consumption (MWh)")
  ) +
  annotation_custom(grobTree(textGrob("Cor = -0.65", x=0.08,  y=0.87, hjust=0,
                                      gp=gpar(col="black", fontsize=13)))) +
  theme_ipsum() +
  theme(
    axis.title.y = element_text(color = desColor, size=17),
    axis.title.y.right = element_text(color = luzColor, size=17)
  )

ggsave("objects/des_luz.pdf", p_1, width = 12, height = 6)

# PLOT 2 - FIGURE 3

# Data Imp
metereologicos <- read.csv("https://raw.githubusercontent.com/pedroskorin/L2_Boosting/master/l2_boosting/data/history/base_bruta/metereologicos.csv")

# Data Processing
date <- seq(as.Date("2002/2/1"), as.Date("2017/10/1"), "months")
data2 <- data.frame("Time"=date,
                   "temp"=diff(log(head(metereologicos$torres.TempMaximaMedia,-2))),
                   "luz"=diff(log(Y_test$exp)))

# Value used to transform the data
coeff2 <- 1/2

p_2 = ggplot(data2, aes(x=Time)) +
  geom_line( aes(y=temp), size=0.75, color="darkorange2") + 
  geom_line( aes(y=luz/coeff2), size=0.5) +
  scale_y_continuous(
    # Features of the first axis
    name = "Log diff of AMT in Torres (Cº)",
    # Add a second axis and specify its features
    sec.axis = sec_axis(~.*coeff2, name="Log diff of Elecricity Consumption (MWh)")
  ) +
  annotation_custom(grobTree(textGrob("Cor = 0.48", x=0.08,  y=0.9, hjust=0,
                                      gp=gpar(col="black", fontsize=13)))) +
  theme_ipsum() +
  theme(
    axis.title.y = element_text(color = "darkorange2", size=17),
    axis.title.y.right = element_text(color = "black", size=17)
  )

ggsave("objects/temp_luz.pdf", p_2, width = 12, height = 6)
