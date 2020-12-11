
Y_or = read.csv("https://raw.githubusercontent.com/pedroskorin/L2_Boosting/master/l2_boosting/data/target.csv",
                               encoding = "UTF-8")[,4]
Y = read.csv("https://raw.githubusercontent.com/pedroskorin/L2_Boosting/master/l2_boosting/data/target.csv",
                                    encoding = "UTF-8")[-1,3]
X = read.csv("https://raw.githubusercontent.com/pedroskorin/L2_Boosting/master/l2_boosting/data/predictors.csv",
                               encoding = "UTF-8")[,-c(1,2)]

library(forecast)
library(mboost)

source("l2_boosting/functions.R")

# Example with v=0.1, h=1, M_max = 2500 and ratio_start = 80%

v_in = 0.1
h_in = 3
Mstop_in = 2500
ratio_start_lag = 0.8

b = prediciton_boost_2(Y_or_lag,
                       Y_lag,
                       X,
                       v = v_in,
                       h = h_in,
                       ratio_start = ratio_start_lag,
                       Mstop = Mstop_in)

b_2 = SARIMA_bench(Y_or, Y, h=3, ratio_start = 0.8)

n_tot_lag <- length(Y_lag)
n_out_lag <- ceiling(n_tot_lag - ratio_start_lag*n_tot_lag)
ind_out_lag <- seq(to = n_tot_lag, by = 1, length = n_out_lag)

n_tot <- length(Y)
n_out <- ceiling(n_tot - ratio_start_in*n_tot)
ind_out <- seq(to = n_tot, by = 1, length = n_out)

# RESULTS FOR BOOSTING

evaluation(b$forecast_uni, Y_or_lag, ind_out_lag, "boost")

# RESULTS FOR BENCH

evaluation(b_2$benchmark, Y_or, ind_out, "SARIMA")

# Variable Importance
## Varimp

var_imp = rowSums(b$varimp)/ncol(b$varimp)

names(var_imp) = c("Intercept", colnames(X))

df = data.frame(imp = var_imp[order(var_imp, decreasing = T)],
                name = names(var_imp[order(var_imp, decreasing = T)]))

## Frequency

n_prev = length(b$selected)
freq_search = function(x) {
  y=0
  for (i in 1:n_prev) {
    y = y + (x %in% b$selected[[i]])
  }
  return(y/n_prev)
}

vector_freq = mapply(freq_search, 1:(ncol(X)+1))
df_freq = data.frame(freq = vector_freq[order(vector_freq, decreasing = T)],
                     index = (1:ncol(X))[order(vector_freq, decreasing = T)])
View(df_freq)

# Plotting graph

library(ggplot2)
library(reshape2)

# Plot

date = seq(as.Date("2003/2/1"), as.Date("2017/10/1"), "months")

intervalo = 48

Boosting = exp(tail(Y_or, 48))
SARIMA = exp(tail(Y_or, 48))
Original = exp(tail(Y_or, 48))
date = tail(date, 48)

n = length(tail(exp(b$forecast_uni),-1))

Boosting[c(rep(F,length(date)-n),rep(T,n))] = tail(exp(b$forecast_uni),-1)
SARIMA[c(rep(F,length(date)-n),rep(T,n))] = tail(exp(b_2$benchmark),-1)

df = data.frame(date = date, SARIMA = SARIMA, 
                Boosting  = Boosting, Original = Original)

meltdf = melt(df, id = "date")
colnames(meltdf)[2] = "Model"

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