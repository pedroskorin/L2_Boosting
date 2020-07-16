### Comentario sobre previsao ###

f.ols <- function(X,y,x0) {
  y <- as.matrix(y)
  if(ncol(x0) != ncol(X)){
    stop("x0 deve ser matriz")
  }
  b <- solve(a = t(X)%*%X, b = t(X)%*%y)
  y.for <- x0%*%b
  return(y.for)
}

## DGP ##

set.seed(1234)

n <- 1000

e <- rnorm(n,0,0.2)
ex <- rnorm(n,0,0.2)

x <- 0
y <- 0
b <- c(1,0.6,0.7)

for(i in 1:n){
  x[i+1] <- 0.55*x[i] + ex[i]
  y[i+1] <- b[1] + b[2]*y[i] + b[3]*x[i+1] + e[i]
}

y <- y[-1]
x <- x[-1]

plot.ts(y)

## Previsao fora da amostra ##

n.tot <- length(y) # numero de obs.
n.out <- 300 # numero de obs. fora da amostra
h <- 1 # horizonte de previsao

ind.out <- seq(to = n.tot, by = 1, length = n.out)
#ind.in <- seq(from = 1, to = ind.out[1] - h, by = 1)

forecast.y <- rep(NA, n.out)

for (i in 1:n.out) {
  
  ind.in <- seq(from = 1, to = ind.out[i] - h, by = 1) #expanding window
  
  y.dep <- y[tail(ind.in,-h)] # variavel dependente t = h+1, ..., T.in
  y.ind <- y[head(ind.in,-h)] # y independente t = 1, ..., T.in-h
  x.ind <- x[head(ind.in,-h)] # x independente t = 1, ..., T.in-h
  
  X.reg <- cbind(1,y.ind,x.ind) # modelo auto-regressivo em y e defasado em x
  y.reg <- as.matrix(y.dep)
  x0.reg <- matrix(c(1,y[tail(ind.in,1)],x[tail(ind.in,1)]), nrow = 1)
  
  forecast.y[i] <- f.ols(X.reg,y.reg,x0.reg)
  
}

plot.ts(y[ind.out])
lines(forecast.y, col = "red")
