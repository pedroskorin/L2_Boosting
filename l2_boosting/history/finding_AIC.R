# Outro

library(mboost)

Y_in = Y[2:152]
X_in = X[1:151,]

Fitting(Y_in, X_in, 1656, 0.3)[2]

X_in = as.matrix(X_in)
colnames(X_in) = colnames(X)

install.packages("caret")



library(caret)

df_imp = varimp(model_1)

plot(df_imp, percent = T)

df_imp_2 = as.data.frame(df_imp)

df_imp_2 = df_imp_2[order(df_imp_2$reduction, decreasing = T),]
View(df_imp_2)
