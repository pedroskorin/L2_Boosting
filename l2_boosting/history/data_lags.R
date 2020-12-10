X_sem_1 = head(X,-1)
X_sem_2 = head(X_sem_1,-1)
X_sem_3 = head(X_sem_2,-1)
X_sem_4 = head(X_sem_3,-1)
X_sem_5 = head(X_sem_4,-1)
X_sem_6 = head(X_sem_5,-1)
X_sem_7 = head(X_sem_6,-1)
X_sem_8 = head(X_sem_7,-1)
X_sem_9 = head(X_sem_8,-1)
X_sem_10 = head(X_sem_9,-1)
X_sem_11 = head(X_sem_10,-1)
X_sem_12 = head(X_sem_11,-1)

X_s_0 = tail(X,nrow(X_sem_12))
X_s_1 = tail(X_sem_1,nrow(X_sem_12))
X_s_2 = tail(X_sem_2,nrow(X_sem_12))
X_s_3 = tail(X_sem_3,nrow(X_sem_12))
X_s_4 = tail(X_sem_4,nrow(X_sem_12))
X_s_5 = tail(X_sem_5,nrow(X_sem_12))
X_s_6 = tail(X_sem_6,nrow(X_sem_12))
X_s_7 = tail(X_sem_7,nrow(X_sem_12))
X_s_8 = tail(X_sem_8,nrow(X_sem_12))
X_s_9 = tail(X_sem_9,nrow(X_sem_12))
X_s_10 = tail(X_sem_10,nrow(X_sem_12))
X_s_11 = tail(X_sem_11,nrow(X_sem_12))
X_s_12 = tail(X_sem_12,nrow(X_sem_12))

X_total = cbind(X_s_0,
                X_s_1,
                X_s_2,
                X_s_3,
                X_s_4,
                X_s_5,
                X_s_6,
                X_s_7,
                X_s_8,
                X_s_9,
                X_s_10,
                X_s_11,
                X_s_12)

Y_total = tail(Y,-12)
Y_or_total = tail(Y_or,-12)

troca_nome = function(x, numero = n) {return(paste("L",n," - ", x, sep = ""))}

nome_all = c()

for (i in 0:12) {
  n = i
  nome_in = unname(lapply(colnames(X), troca_nome))
  nome_all = c(nome_all, nome_in)

}

colnames(X_total) = unlist(nome_all)


X = X_total
Y = Y_total
Y_or = Y_or_total
