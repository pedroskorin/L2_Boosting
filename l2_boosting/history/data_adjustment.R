# Rotina para dessazonalizar e estacionarizar variaveis 
#
# Autores: guislinden e pedroskorin
# UFRGS
# 2020-07-09

# Pacotes
library(readr)
library(urca)

# Retiradas de bases brutas
nacional_mensal = data.frame(read.csv("https://raw.githubusercontent.com/pedroskorin/L2_Boosting/master/l2_boosting/dados/history/base_bruta/nacional_mensal.csv",
                                      encoding = "latin1"))
metereologicos = data.frame(read.csv("https://raw.githubusercontent.com/pedroskorin/L2_Boosting/master/l2_boosting/dados/history/base_bruta/metereologicos.csv",
                                     encoding = "UTF-8"))
internacional = data.frame(read.csv("https://raw.githubusercontent.com/pedroskorin/L2_Boosting/master/l2_boosting/dados/history/base_bruta/internacional.csv",
                                    encoding = "UTF-8"))
regional = data.frame(read.csv("https://raw.githubusercontent.com/pedroskorin/L2_Boosting/master/l2_boosting/dados/history/base_bruta/regional.csv",
                               encoding = "UTF-8"))


nacional_mensal = nacional_mensal[,-1]
metereologicos = metereologicos[,-1:-2]
regional = regional[,-1]
internacional = internacional[,-1]


base_bruta = cbind(nacional_mensal,
                   metereologicos,
                   regional,
                   internacional)

base_bruta = base_bruta[,-167]

## Processo de estacionarizacao ####

base_ponto = base_bruta

## Deteccao de tipo de dado

tipo = c(4)

for (i in 2:ncol(base_ponto)){
  print(i)
  if (sum(base_ponto[,i]<0) > 0) {
    
    if (sum(base_ponto[,i]==0) > 0) {
      
      tipo = append(tipo, 3)
      
    } else {
    
    tipo = append(tipo, 1) 
    
    }
    
  } else {
  
  if (sum(base_ponto[,i]==0) > 0){
    tipo = append(tipo, 2)
  } else {
    
    tipo = append(tipo, 0)
  }
  
  }
  
}

# Deteccao de quantidade de lags necessarios ####

test = c(0)

preparacao = function(X, i) {
  
  if (tipo[i] == 0) {
    
    return(log(X))
  }
  
  if (tipo[i] == 2) {
    
    return(log(X+2))
    
  }
  
  if (tipo[i] == 1) {
    
    return(X)
  }
  
  if (tipo[i] == 3) {
    
    return(X+1)
  }
}

cresc_discreto = function(X) {
  
  Y = c()
  
  for (i in 2:length(X)) {
    
    y = X[i]/X[i-1]-1
    
    Y = append(Y, y)
    
  }
  
  return(Y)
}

for(i in 2:ncol(base_ponto)){
  print(i)
  X = base_ponto[,i]
  
  k = 0
  
  j = 0
  
  status = "non-stationary"
  
  while (status == "non-stationary"){
    
    adf_test = summary(ur.df(X, "none"))
    if(k == 2) {
      dale = colnames(base_ponto)[i]
    }
    #adf.test(X)[4] <= 0.05
    
    if (adf_test@teststat[1] <= adf_test@cval[1,1]){
      
      status = "stationary"
      
      if (j == 1){
        k = 0.5
      }
      
    } else {
      
      if (j == 0) {
        
        X = preparacao(X,i)
        
        j = 1
        
      } else {
        
      k = k + 1
      
      if (tipo[i] == 0 | tipo[i] == 2) {
        X = diff(X)
      }
      if (tipo[i] == 1 | tipo[i] == 3) {
        X = cresc_discreto(X)
      }
      
      j = j+1
      
      }
      
    }
    
    
    print(k)
  }
  
  test = append(test, k)
  
}

## Plot de casos necessidades de lags especificas 

for (i in which(test == 2)) {
  
  plot((base_ponto[,i]), type = "l", xlab = as.character(i),
       main = colnames(base_ponto)[i])
  
}

table(test)

# Aplicacao dos Lags ####

base_estacionaria = as.data.frame(base_ponto[-1,1])

for (i in 2:ncol(base_ponto)) {
  print(i)
  X = base_ponto[,i]
  
  if (test[i] == 0) {
    
    #X = X[-1:-2]
    X = X[-1]
  }
  
  if (test[i] == 0.5) {
    
    X = preparacao(X, i)
    
    #X = X[-1:-2]
    X = X[-1]
    
  }
  
  if (test[i] == 1) {
    
    X = preparacao(X, i)
    
    if (tipo[i] == 0 | tipo[i] == 2) {
      X = diff(X)
    }
    if (tipo[i] == 1 | tipo[i] == 3) {
      X = cresc_discreto(X)
    }
    
#    X = X[-1]
  }
  
  if (test[i] == 2) {
    
    X = preparacao(X, i)
    
    if (tipo[i] == 0 | tipo[i] == 2) {
      X = diff(diff(X))
    }
    if (tipo[i] == 1 | tipo[i] == 3) {
      X = crescimento_discreto(cresc_discreto(X))
    }
    
  }

  base_estacionaria = cbind(base_estacionaria, X)
  
}

# Nomeacao das colunas
colnames(base_estacionaria) = colnames(base_ponto)
colnames(base_estacionaria)[1] = "Data"
#base_estacionaria = base_estacionaria[-191,]


write.csv(base_estacionaria,"C:\\Users\\Pedro.Pedro-PC\\Documents\\GitHub\\L2_Boosting\\l2_boosting\\dados\\predictors.csv", fileEncoding = "UTF-8" )

write.csv(Y_dados, "C:\\Users\\Pedro.Pedro-PC\\Documents\\GitHub\\L2_Boosting\\l2_boosting\\dados\\target.csv", fileEncoding = "UTF-8")
