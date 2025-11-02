
# rho estimator (each cluster)
rho_optim <- function(u, var.u, D, X, corr_type, id){
  rho_hat <- numeric(length(unique(id)))
  for(i in unique(id)){
    h_i <- function(rho){
      R <- R_i(length(id[id==i]), rho, corr_type)
      W_i <- prod_2(R, (var.u[id==i])^(1/2), D[id==i])
      zs_i <- D[id==i]^(-1)*u[id==i]
      sum_i <- prod_1(t(X[id==i,]), prod_1(W_i, zs_i))
      return(sum(abs(sum_i)))  
    }
    rho_hat[i] <- optim(par=0, fn=h_i, lower=-1, upper=1, method="Brent")$par
  }
  return(mean(rho_hat))
}

# rho estimator (all clusters)
rho_optim <- function(u, var.u, D, X, corr_type, time, id){
  h <- function(rho){
    sum_2 <- 0
    for(i in unique(id)){
      R <- R_i(rho, corr_type, time[id==i])
      W_i <- prod_2(R, (var.u[id==i])^(1/2), D[id==i])
      zs_i <- D[id==i]^(-1)*u[id==i]
      sum_2 <-  sum_2 + prod_1(t(X[id==i,]), prod_1(W_i, zs_i))
    }
    return(sum(abs(sum_2)))  
  }
  return(optim(par=0, fn=h, lower=-1/(max(time)-1), upper=1, method="Brent")$par)
}
