#L2 BOOSTING ALGORITHM
#
#Implementation of the L2 Boosting algorithm
#
#theoretical source: "Jing Zang, 2014, Forecasting Aggregates with Disaggregate Variables: 
#Does boosting help to select the most informative predictors?"
#
#Authors: guislinden, pedroskorin
#UFRGS
#
#2019-6-19

# Getting data from github without the first column
data = read.csv("https://raw.githubusercontent.com/pedroskorin/L2_Boosting/master/l2_boosting/dados.csv")[,-1]

# Selecting data
Y = data$Y # Response variable

X = data[,-ncol(data)] # Predictors variables

# Creating function L2_boost:
# INPUT: Y - response ; X - predictors ; M - Iterations ; v - shrinkage parameter (standart 0.1)
# OUTPUT : MSE of algorithm ; times each predictor was choosed

L2_boost = function(Y,X,M,v=0.1){

# Defining variables
ft = mean(Y) # First appearence of function ft
m = 0 # m index

u = 0 # error vector (Y-ft)

b = 0 # optimum coefficient for each predictor
teta = c() # vector of optimum coefficient for all possible predictors

sum_squared_resid = 0 # SSR for regression of ut for each predictor
SSR = c() # list of all SSR

g = c() # vector of selected predictor * optimum coefficient
g_optimum = list() # list of all g's
index = 0 # index of selected predictor in teta's vector

X_optimum = c() # Vector of selected predictor
Matrix = c() # Matrix of all selected predictors*coefficient*v
f_optimum = c() # final prediction vector

choosed_predictors = rep(0,ncol(X)) # controling choosed predictors
names(choosed_predictors) = names(X)

while (m < M) { # Iterate M times
u = Y - ft # Calculating error in m

# calculating optimum coefficients
teta = c()

for (i in X){
  b = solve(t(i)%*%i)%*%t(i)%*%Y
  
  teta = c(teta, b)
}

# calculating SSR

SSR = c()

for (i in 1:length(teta)){
  sum_squared_resid = sum((u - teta[i]*X[,i])^2)
  
  SSR = c(SSR, sum_squared_resid)
  }

# Minimizing SSR

index = which.min(SSR)

X_optimum = X[,index]

# Calculating g (vector of selected predictor * optimum coefficient)
g = teta[index]*X_optimum

choosed_predictors[index] = choosed_predictors[index] + 1

# Updating ft
ft = ft + v*g

# Iterating
m = m + 1

g_optimum[[m]] <- g

}

# Matrix of all g vectors
Matrix = do.call(cbind, g_optimum)

# Final function
f_optimum = mean(Y) + ft

# MSE of L2_boost
print(sum((f_optimum - Y)^2)/nrow(X))

print(choosed_predictors)
}

#Example
L2_boost(Y,X,20)

# Hudson, ficamos com dúvidas com a relação M x MSE
# Intuitivamente, parece fazer sentido que quanto mais interações, mais específico fica o modelo
# E menor o MSE, entretanto, como visto no gráfico, o MSE parece aumentar
# Não sabemos se isso se deve ao fato de que algo está errado no algoritmo nosso
# Ou se é uma característica normal do modelo

