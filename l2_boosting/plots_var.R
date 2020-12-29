library(ggplot2)
library(dplyr)
library(patchwork) # To display 2 charts together
library(hrbrthemes)
library(reshape2)
library(grid)
library(extrafont)
loadfonts()

# PLOT 1 - FIGURE 2

Y_test = read.csv("https://raw.githubusercontent.com/pedroskorin/L2_Boosting/master/l2_boosting/data/target.csv",
                  encoding = "UTF-8")
Y_test$exp = exp(Y_test$Y_or)

desemprego <- read.csv("https://raw.githubusercontent.com/pedroskorin/L2_Boosting/master/l2_boosting/data/history/ipeadata[27-12-2020-01-20].csv")
colnames(desemprego) = c("data", "des", "nada")
des_2 = desemprego$des[206:395]

date = seq(as.Date("2002/1/1"), as.Date("2017/10/1"), "months")

data = data.frame("Time"=date,"des"=des_2,"luz"=Y_test$exp)

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

ggsave("des_luz.pdf", p_1, width = 12, height = 6)

# PLOT 2 - FIGURE 3

metereologicos = read.csv("https://raw.githubusercontent.com/pedroskorin/L2_Boosting/master/l2_boosting/data/history/base_bruta/metereologicos.csv")

date = seq(as.Date("2002/2/1"), as.Date("2017/10/1"), "months")
data2 = data.frame("Time"=date,
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

ggsave("temp_luz.pdf", p_2, width = 12, height = 6)
