# density 
dul <- function(y, mu){
  d <- (((1-mu)^2)/(mu*(1-y)^3))*exp(-y*(1-mu)/(mu*(1-y)))
  d[y<=0] <- 0; d[y>=1] <- 0
  return(d)
}

# distribution
pul <- function(q, mu){
  p <- numeric(length(q))
  p <- 1-((1-mu*q)/(1-q))*exp(-q*(1-mu)/(mu*(1-q)))
  p[q<=0] <- 0; p[q>=1] <- 1
  return(p)
}

# quantile
qul <- function(p, mu){
  pos <- mu > 0.002; q <- numeric(length(p))
  W1 <- lamW::lambertWm1((p[pos]-1)*exp(-1/mu[pos])/mu[pos])
  q[pos] <- (1+W1*mu[pos])/((W1+1)*mu[pos])
  q[!pos] <- Vectorize(function(p, mu) 
    uniroot((function(x) pul(x, mu) - p), c
            (0+.Machine$double.eps,1-.Machine$double.eps), 
            f.upper=1-.Machine$double.eps, tol=1e-10)$root)(p[!pos], mu[!pos])
  return(unlist(q))
}

# variance
ul.var <- function(mu){
  E1 <- expint::expint_En(1/mu-1, 1)
  return(unlist(((1-mu)^2)*(E1*exp(1/mu-1) - mu)/mu))
}

# random generate (indep.)
rul <- function(n, mu){
  return(qul(runif(n), mu))
}

# make Ri structure
R_i <- function(rho, corr_type, time_i){
  R <- matrix(1)
  if(length(time_i) >= 2){
    if(corr_type == "IDE"){
      R <- matrix(0, length(time_i), length(time_i))
      diag(R) <- 1
    }
    
    if(corr_type == "AR1"){
      R <- rho^abs(matrix(time_i-1, nrow=length(time_i), ncol=length(time_i),
                          byrow=T) - (time_i-1))
    }
    
    if(corr_type == "EXC"){
      R <- matrix(rho, length(time_i), length(time_i))
      diag(R) <- 1
    }
    
    if(corr_type == "UNS"){
      n_i <- (length(rho))^(1/2)
      R <- matrix(rho, n_i, n_i)[time_i-min(time_i)+1,time_i-min(time_i)+1]
    }
  }
  
  # return
  return(R)
}

# random generate (corr. + unbalanced data)
rul_corr <- function(mu, rho, corr_type, time, id, show.step=T){
  # uniform corr.
  if(show.step) pbapply::pboptions(txt.width=50, style=3, char="=")
  ids <- unique(id)
  if(show.step) pb <- pbapply::startpb(min=0, max=length(ids))
  unif_corr <- unlist(lapply(ids, function(i){
    pos <- id==i; if(show.step) pbapply::setpb(pb, i)
    pnorm(mvnfast::rmvn(1, rep(0, length(id[pos])), 
                        R_i(rho, corr_type, time[pos])))
    }))
  
  # return
  return(qul(unif_corr, mu))
}

# rho moments
rho_moments <- function(u, var.u, corr_type, time, id){
  if(corr_type == "IDE"){
    rho_hat <- 0
    return(rho_hat)
  }
  
  # prep. correlation
  r <- u/sqrt(var.u); rho_hat <- 0
  ids <- as.numeric(names(table(id))[table(id)>=2])
  
  if(corr_type == "AR1"){
    for(i in ids){
      pos <- id==i; time_i <- time[pos]  
      r_i <- numeric(max(time_i)-min(time_i)+1)
      r_i[time_i-min(time_i)+1] <- r[pos]
      aux <- outer(r_i, r_i, FUN="*")
      rho_hat <- rho_hat + mean(pracma::Diag(aux, 1))
    }
    rho_hat <- rho_hat/length(ids)
    rho_hat[rho_hat>1] <- 0.5; rho_hat[rho_hat<(-1)] <- -1/(2*length(r_i))
  }

  if(corr_type == "EXC"){
    for(i in ids){
      pos <- id==i; r_i <- r[pos]
      aux <- outer(r_i, r_i, FUN="*"); diag(aux) <- 0
      aux <- sum(aux)/(length(r_i)*(length(r_i)-1))
      rho_hat <- rho_hat + aux
    }
    rho_hat <- rho_hat/length(ids)
    rho_hat[rho_hat>1] <- 0.5; rho_hat[rho_hat<(-1)] <- -1/(2*length(r_i))
  }
  
  if(corr_type == "UNS"){
    n <- 0
    for(i in ids){
      pos <- id==i; time_i <- time[pos]
      r_i <- numeric(max(time)-min(time)+1)
      r_i[time_i] <- r[pos]
      aux <- outer(r_i, r_i, FUN="*")
      n <- n + 1*(aux!=0)
      rho_hat <- rho_hat + aux
    }
    rho_hat <- rho_hat/n; diag(rho_hat) <- 1
    rho_hat <- as.numeric(rho_hat)
    rho_hat[rho_hat>1] <- 0.5; rho_hat[rho_hat<(-1)] <- -1/(2*length(r_i))
  }
  
  # return
  return(rho_hat)
}

# newton score unit lindley
ulgee <- function(y, X, time, id, corr_type, link, epsilon_1, max.iter, 
                  show.step=T){
  if(link=="logit"){
    # link 
    g <- function(mu){
      return(gamlss.dist::qLO(p=mu)) 
    }
    
    # inv. link 
    inv_g <- function(eta){
      thresh <- -gamlss.dist::qLO(.Machine$double.eps)
      eta <- pmin(thresh, pmax(eta, -thresh))
      return(gamlss.dist::pLO(eta))
    }
    
    # der. link 
    dmu.deta <- function(eta){
      return(pmax(gamlss.dist::dLO(eta), .Machine$double.eps))
    }
  }
  
  if(link=="probit"){
    # link
    g <- function(mu){
      return(gamlss.dist::qNO(p=mu)) 
    }
    
    # inv. link 
    inv_g <- function(eta){
      thresh <- -gamlss.dist::qNO(.Machine$double.eps)
      eta <- pmin(thresh, pmax(eta, -thresh))
      return(gamlss.dist::pNO(eta))
    }
    
    # der. link 
    dmu.deta <- function(eta){
      return(pmax(gamlss.dist::dNO(eta), .Machine$double.eps))
    }
  }
  
  if(link=="cloglog"){
    # link 
    g <- function(mu){
      return(gamlss.dist::qGU(p=mu)) 
    }
    
    # inv. link 
    inv_g <- function(eta){
      thresh <- -gamlss.dist::qGU(.Machine$double.eps)
      eta <- pmin(thresh, pmax(eta, -thresh))
      return(gamlss.dist::pGU(eta))
    }
    
    # der. link
    dmu.deta <- function(eta){
      return(pmax(gamlss.dist::dGU(eta), .Machine$double.eps))
    }
  }
  
  # not appropriate data
  if(sum(y<=0)!=0 | sum(y>=1)!=0){
    # return
    output <- list(FALSE); names(output) <- c("converged"); return(output)
  }
  
  # iterative process
  beta_2 <- as.numeric(coef(lm(g(y)~X-1)))
  epsilon_2 <- 10; iter <- 1
  if(show.step) cat("Iteration", iter-1, "of", max.iter,
                    "- evaluated criterium =", epsilon_2,"\r")
  while(epsilon_2 > epsilon_1 & iter <= max.iter){
    beta_1 <- beta_2
    eta <- prod_1(X, beta_1)
    mu <- inv_g(eta)
    u <- -2/(1-mu) - 1/mu + (1/mu^2)*(y/(1-y)) 
    var.u <- (2 - (1 - mu)^2)/((mu*(1 - mu))^2)
    D <- -var.u*dmu.deta(eta)
    t <- eta - u/D
    rho <- rho_moments(u, var.u, corr_type, time, id)
    W <- list(); H_1i <- 0; sum_1 <- 0; sum_2 <- 0; sum_3 <- 0
    for(i in unique(id)){
      pos <- id==i; Di <- D[pos]; Xi <- rbind(X[pos,]); ti <- t[pos]
      if(ncol(X)==1) Xi <- cbind(X[pos,])
      Ri <- R_i(rho, corr_type, time[pos]) 
      Wi <- prod_2(Ri, sqrt(var.u[pos]), Di)
      W[[length(W)+1]] <- Wi
      rpi <- prod_1(Wi, eta[pos]-ti)
      H_1i <- H_1i + prod_3(Xi, diag(Di^2))
      sum_1 <- sum_1 + prod_3(Xi, Wi)
      sum_2 <- sum_2 + prod_1(t(Xi), prod_1(Wi, ti))  
      sum_3 <- sum_3 + prod_3(Xi, rpi%*%t(rpi))
    }
    beta_2 <- prod_4(sum_1, sum_2)
    if(sum(is.na(beta_2))!=0){iter <- max.iter; break}
    epsilon_2 <- max(abs((beta_2 - beta_1)/beta_1)) 
    if(show.step) cat("Iteration", iter, "of",
                      max.iter, "- evaluated criterium =", epsilon_2,"\r")
    iter <- iter + 1
  }

  # converged
  if(iter < max.iter & sum(pul(y, mu) >= 1) == 0 &
     sum(pul(y, mu) <= 0) == 0 & sum(rho[rho!=1]==0.5) == 0){
    # p-values
    var.beta <- prod_5(sum_1, sum_3)
    wald.stat <- (beta_1^2)/diag(var.beta) 
    pvalue <- 1-pchisq(wald.stat, 1)
    
    # qic
    qic <- -2*sum(log(dul(y, mu))) + 2*sum(diag(var.beta%*%H_1i))
    
    # sensitivity case-weight
    t.delta_c <- Matrix::Diagonal(x=u*(D^(-1)))%*%Matrix::.bdiag(W)%*%X
    B_c.tr <- sum(Matrix::diag((t.delta_c%*%solve(sum_1))%*%
                                 (Matrix::t(t.delta_c)%*%
                                    (t.delta_c%*%solve(sum_1)))%*%
                                 Matrix::t(t.delta_c)))
    Bi_c <- Matrix::diag((t.delta_c%*%solve(sum_1))%*%
                           Matrix::t(t.delta_c))/sqrt(B_c.tr)
    
    # sensitivity response
    f <- sqrt(ul.var(mu))/((1-y)*mu)^2
    t.delta_r <- Matrix::Diagonal(x=f*(D^(-1)))%*%Matrix::.bdiag(W)%*%X
    B_r.tr <- sum(Matrix::diag((t.delta_r%*%solve(sum_1))%*%
                                 (Matrix::t(t.delta_r)%*%
                                    (t.delta_r%*%solve(sum_1)))%*%
                                 Matrix::t(t.delta_r)))
    Bi_r <- Matrix::diag((t.delta_r%*%solve(sum_1))%*%
                           Matrix::t(t.delta_r))/sqrt(B_r.tr)
    
    # return 
    output <- list(beta_1, var.beta, pvalue, rho, mu, qnorm(pul(y, mu), 0, 1),
                   Bi_c, Bi_r, qic, corr_type, time, id, link, TRUE)
    names(output) <- c("mu.coefs", "vcov", "pvalues", "rho", "mu.hat", "rq", "Bi_c",
                       "Bi_r", "qic", "corr_type", "time", "id", "link", "converged")
    return(output)
  }
  
  # not converged
  else{
    # return
    output <- list(FALSE); names(output) <- c("converged"); return(output)
  }
}

# newton score unit lindley (just parameters)
ulgee.fast <- function(y, X, time, id, corr_type, link, epsilon_1, max.iter){
  if(link=="logit"){
    # link 
    g <- function(mu){
      return(gamlss.dist::qLO(p=mu)) 
    }
    
    # inv. link 
    inv_g <- function(eta){
      thresh <- -gamlss.dist::qLO(.Machine$double.eps)
      eta <- pmin(thresh, pmax(eta, -thresh))
      return(gamlss.dist::pLO(eta))
    }
    
    # der. link 
    dmu.deta <- function(eta){
      return(pmax(gamlss.dist::dLO(eta), .Machine$double.eps))
    }
  }
  
  if(link=="probit"){
    # link
    g <- function(mu){
      return(gamlss.dist::qNO(p=mu)) 
    }
    
    # inv. link 
    inv_g <- function(eta){
      thresh <- -gamlss.dist::qNO(.Machine$double.eps)
      eta <- pmin(thresh, pmax(eta, -thresh))
      return(gamlss.dist::pNO(eta))
    }
    
    # der. link 
    dmu.deta <- function(eta){
      return(pmax(gamlss.dist::dNO(eta), .Machine$double.eps))
    }
  }
  
  if(link=="cloglog"){
    # link 
    g <- function(mu){
      return(gamlss.dist::qGU(p=mu)) 
    }
    
    # inv. link 
    inv_g <- function(eta){
      thresh <- -gamlss.dist::qGU(.Machine$double.eps)
      eta <- pmin(thresh, pmax(eta, -thresh))
      return(gamlss.dist::pGU(eta))
    }
    
    # der. link
    dmu.deta <- function(eta){
      return(pmax(gamlss.dist::dGU(eta), .Machine$double.eps))
    }
  }
  
  # not appropriate data
  if(sum(y<=0)!=0 | sum(y>=1)!=0){
    # return
    output <- list(FALSE); names(output) <- c("converged"); return(output)
  }
  
  # iterative process
  beta_2 <- as.numeric(coef(lm(g(y)~X-1)))
  epsilon_2 <- 10; iter <- 1
  while(epsilon_2 > epsilon_1 & iter <= max.iter){
    beta_1 <- beta_2
    eta <- prod_1(X, beta_1)
    mu <- inv_g(eta)
    u <- -2/(1-mu) - 1/mu + (1/mu^2)*(y/(1-y)) 
    var.u <- (2 - (1 - mu)^2)/((mu*(1 - mu))^2)
    D <- -var.u*dmu.deta(eta)
    t <- eta - u/D
    rho <- rho_moments(u, var.u, corr_type, time, id)
    sum_1 <- 0; sum_2 <- 0
    for(i in unique(id)){
      pos <- id==i; Di <- D[pos]; Xi <- rbind(X[pos,]); ti <- t[pos]
      if(ncol(X)==1) Xi <- cbind(X[pos,])
      Ri <- R_i(rho, corr_type, time[pos]) 
      Wi <- prod_2(Ri, sqrt(var.u[pos]), Di)
      sum_1 <- sum_1 + prod_3(Xi, Wi)
      sum_2 <- sum_2 + prod_1(t(Xi), prod_1(Wi, ti))  
    }
    beta_2 <- prod_4(sum_1, sum_2)
    if(sum(is.na(beta_2))!=0){iter <- max.iter; break}
    epsilon_2 <- max(abs((beta_2 - beta_1)/beta_1)) 
    iter <- iter + 1
  }
  
  # converged
  if(iter < max.iter & sum(pul(y, mu) >= 1) == 0 & sum(pul(y, mu) <= 0) == 0 & 
     sum(rho[rho!=1]==0.5) == 0){
    output <- list(beta_1, rho, mu, qnorm(pul(y, mu), 0, 1), corr_type, 
                   time, id, link, TRUE)
    names(output) <- c("mu.coefs", "rho", "mu.hat", "rq", "corr_type", "time", 
                       "id", "link", "converged")
    return(output)
  }
  
  # not converged
  else{
    # return
    output <- list(FALSE); names(output) <- c("converged"); return(output)
  }
}

# auto simulation
ul_sim <- function(beta, R, rho, n, s, corr_type1, corr_type2, show.step=T){
  # prep. tables
  table1_sim <- matrix(NA, length(rho)*length(s)*length(n), 3*(length(beta)+2))
  table1_sim[,1] <- rep(rep(rho, each=length(s)), length(n)) 
  table1_sim[,2] <- rep(n, each=length(s)*length(rho))
  table1_sim[,3] <- rep(rep(s, length(rho)), length(n))
  table2_sim <- table1_sim
  if(show.step) pbapply::pboptions(txt.width=50, style=3, char="=")
  if(show.step) pb <- pbapply::startpb(min=0, max=nrow(table1_sim)*R)
  
  # simulation
  for(j in 1:nrow(table1_sim)){
    table1_coefs <- matrix(NA, R, length(beta)+1)
    table2_coefs <- matrix(NA, R, length(beta)+1)
    for(i in 1:R){
      # data generate
      rho_sim <- table1_sim[j,1]
      n_sim <- table1_sim[j,2]
      s_sim <- table1_sim[j,3]
      id <- rep(1:n_sim, each=s_sim)
      time <- rep(1:s_sim, n_sim)
      x <- runif(n_sim*s_sim, 0, 1)  
      data <- cbind(1, x, time, id)
      
      # mu simulation (probit link)
      eta <- data[,-c(3,4)]%*%beta; mu <- gamlss.dist::pNO(eta)
      
      # simulation 
      sim_ok <- 0
      while(sim_ok==0){
        y_sim <- rul_corr(mu, rho_sim, corr_type1, data[,3], data[,4], F)
        fit1_i <- ulgee.fast(y_sim, data[,-c(3,4)], data[,3], data[,4], 
                             corr_type1, "probit", 1e-06, 20)
        fit2_i <- ulgee.fast(y_sim, data[,-c(3,4)], data[,3], data[,4], 
                             corr_type2, "probit", 1e-06, 20)
        if(fit1_i$converged==TRUE & fit2_i$converged==TRUE){sim_ok <- 1}
      }
      
      # mu coefficients
      table1_coefs[i,1:length(beta)] <- fit1_i$mu.coefs
      table2_coefs[i,1:length(beta)] <- fit2_i$mu.coefs
      
      # correlation coefficient
      table1_coefs[i,length(beta)+1] <- fit1_i$rho
      table2_coefs[i,length(beta)+1] <- fit2_i$rho
      if(show.step) pbapply::setpb(pb, (j-1)*R + i)
    }
    
    # update tables 
    theta <- c(beta, rho_sim)
    for(i in 1:ncol(table1_coefs)){
      # fill column 1
      table1_sim[j,(3*i+1)] <- mean(table1_coefs[,i])
      table2_sim[j,(3*i+1)] <- mean(table2_coefs[,i])
      
      # fill column 2
      table1_sim[j,(3*i+2)] <- abs(mean(table1_coefs[,i])-theta[i])
      table2_sim[j,(3*i+2)] <- abs(mean(table2_coefs[,i])-theta[i])
      
      # fill column 3
      table1_sim[j,(3*i+3)] <- mean((table1_coefs[,i]-theta[i])^2)
      table2_sim[j,(3*i+3)] <- mean((table2_coefs[,i]-theta[i])^2)
    }
  }

  # return 
  return(rbind(table1_sim,table2_sim))
}

# diagnostic quantile residual 
diag_quant <- function(fit.model, X, nsim, show.step=T, random=F, n, label.id, 
                       label.time){
  # random selection
  if(random==T){
    pos <- rep(T, nrow(X))
    ids <- unique(fit.model$id)
    random_time <- numeric(nrow(X))
    for(i in ids){
      pos_i <- rep(0, nrow(X))
      time_i <- sample(fit.model$time[fit.model$id==i], 1)
      pos_i[fit.model$time==time_i & fit.model$id==i] <- 1
      random_time <- random_time + pos_i
    }
    pos <- as.logical(random_time)
    
    # plot quantile residual vs value adjusted
    par(mar=c(5.5,5.5,2,2), mfrow=c(1,2))
    plot(fit.model$mu.hat[pos], fit.model$rq[pos], xlab="Fitted value", 
         ylab="Quantile residual", pch=16, lwd=2, cex.lab=1.5, cex.axis=1.2)
    if(n > 0) identify(x=fit.model$mu.hat[pos], y=fit.model$rq[pos], 
                       label=paste(paste("(",paste(label.id[pos], label.time[pos],
                                                   sep=","), sep=""),")",sep=""),
                       n=n, cex=1.2)
    
    # plot qq-norm quantile residual 
    faixa <- range(fit.model$rq[pos]) 
    qqnorm(fit.model$rq[pos], xlab="Quantile of N(0,1)", ylab="Quantile residual",
           pch=16, lwd=2, cex.lab=1.5, cex.axis=1.2, ylim=faixa, main="")
    
    # plot reference line
    par(new=TRUE)
    abline(a=0, b=1, xlab="", ylab="", lty=2, lwd=1, main="")
  }
  
  # non random selection
  else{
    # plot quantile residual vs value adjusted
    par(mar=c(5.5,5.5,2,2), mfrow=c(1,2))
    plot(fit.model$mu.hat, fit.model$rq, xlab="Fitted value", 
         ylab="Quantile residual", pch=16, lwd=2, cex.lab=1.5, cex.axis=1.2)
    if(n > 0) identify(x=fit.model$mu.hat, y=fit.model$rq, 
                       label=paste(paste("(",paste(label.id, label.time, sep=","),
                                         sep=""),")",sep=""), n=n, cex=1.2)
    
    # create band
    if(nsim > 0){
      # prep. simulation
      if(show.step) pbapply::pboptions(txt.width=50, style=3, char="=")
      if(show.step) pb <- pbapply::startpb(min=0, max=nsim)
      e <- matrix(0, nrow(X), nsim); med <- numeric(nrow(X))
      e1 <- numeric(nrow(X)); e2 <- numeric(nrow(X))
      
      # simulation 
      for(j in 1:nsim){
        sim_ok <- 0
        while(sim_ok==0){
          y_sim <- rul_corr(fit.model$mu.hat, fit.model$rho, fit.model$corr_type, 
                            fit.model$time, fit.model$id, F)
          fit_j <- ulgee.fast(y_sim, X, fit.model$time, fit.model$id,
                              fit.model$corr_type, fit.model$link, 1e-6, 20)
          if(fit_j$converged==TRUE){sim_ok <- 1}
        }
        e[,j] <- sort(fit_j$rq)
        if(show.step) pbapply::setpb(pb, j)
      }
      
      # cut simulation
      for(i in 1:nrow(X)){
        e1[i] <- quantile(e[i,], probs=(0.995))
        e2[i] <- quantile(e[i,], probs=(0.005))
        med[i] <- quantile(e[i,], probs=(0.5))
      }
    
      # plot qq-norm quantile residual 
      faixa <- range(e1, e2) 
      qqnorm(fit.model$rq, xlab="Quantile of N(0,1)", ylab="Quantile residual",
             pch=16, lwd=2, cex.lab=1.5, cex.axis=1.2, ylim=faixa, main="")
      
      # plot qq-norm quantile residual (with band)
      par(new=TRUE)
      qqnorm(e1, axes=FALSE, xlab="", ylab="", type="l", ylim=faixa, lty=1,
             main="", lwd=1)
      par(new=TRUE)
      qqnorm(e2, axes=FALSE, xlab="", ylab="", type="l", ylim=faixa, lty=1,
             main="", lwd=1)
      par(new=TRUE)
      qqnorm(med, axes=FALSE, xlab="", ylab="", type="l", ylim=faixa, lty=2, 
             main="", lwd=1)
    }
  }
}

# sensitivity via conformal normal curvature 
sens_conf <- function(fit.model, c_c, n_c, c_r, n_r, label.id, label.time){
  # plot Bi_c vs index
  par(mar=c(5.5,5.5,2,2), mfrow=c(1,2))
  plot(fit.model$Bi_c, ylab=expression("B"["ij"]), xlab="Index", cex.lab=1.5,
       cex.axis=1.2, pch=16)
  if(n_c > 0) identify(x=1:length(fit.model$Bi_c), y=fit.model$Bi_c, 
                       label=paste(paste("(",paste(label.id, label.time, sep=","),
                                         sep=""),")",sep=""), n=n_c, cex=1.2)
  abline(a=mean(fit.model$Bi_c) + c_c*sd(fit.model$Bi_c), b=0, lty=2, lwd=2)
  
  # plot Bi_r vs index
  plot(fit.model$Bi_r, ylab=expression("B"["ij"]), xlab="Index", cex.lab=1.5,
       cex.axis=1.2, pch=16)
  if(n_r > 0) identify(x=1:length(fit.model$Bi_r), y=fit.model$Bi_r, 
                       label=paste(paste("(",paste(label.id, label.time, sep=","),
                                         sep=""),")",sep=""), n=n_r, cex=1.2)
  abline(a=mean(fit.model$Bi_r) + c_r*sd(fit.model$Bi_r), b=0, lty=2, lwd=2)
}

sens_mrc <- function(fit.model, y, X, n.sample, pos, show.step=T){
  # prep
  table <- matrix(NA, n.sample+1, 3)
  fit_i <- ulgee.fast(y[!pos], X[!pos,], fit.model$time[!pos], fit.model$id[!pos],
                      fit.model$corr_type, fit.model$link, 1e-6, 40)
  table[,1] <- 0:n.sample
  table[1,2] <- max(abs((fit.model$mu.coefs-fit_i$mu.coefs)/fit.model$mu.coefs))
  table[1,3] <- fit_i$rho
  
  # simulation
  if(show.step) pbapply::pboptions(txt.width=50, style=3, char="=")
  if(show.step) pb <- pbapply::startpb(min=0, max=n.sample)
  for(i in 1:n.sample){
    sim_ok <- 0
    while(sim_ok==0){
      pos_i <- (1:length(y)) %in% sample((1:length(y))[!pos], sum(pos))
      fit_i <- ulgee.fast(y[!pos_i], X[!pos_i,], fit.model$time[!pos_i],
                          fit.model$id[!pos_i], fit.model$corr_type, 
                          fit.model$link, 1e-6, 40)
      if(fit_i$converged==TRUE){sim_ok <- 1}
    }
    table[i+1,2] <- max(abs((fit.model$mu.coefs-fit_i$mu.coefs)/fit.model$mu.coefs))
    table[i+1,3] <- fit_i$rho
    if(show.step) pbapply::setpb(pb, i)
  }
  
  # return
  return(table)
}

sens_coef <- function(fit.model, X, nsim, show.step=T){
  # prep
  table <- matrix(NA, nsim, ncol(X))
  
  # simulation 
  if(show.step) pbapply::pboptions(txt.width=50, style=3, char="=")
  if(show.step) pb <- pbapply::startpb(min=0, max=nsim)
  for(j in 1:nsim){
    sim_ok <- 0
    while(sim_ok==0){
      y_sim <- rul_corr(fit.model$mu.hat, fit.model$rho, fit.model$corr_type, 
                        fit.model$time, fit.model$id, F)
      fit_j <- ulgee.fast(y_sim, X, fit.model$time, fit.model$id, 
                          fit.model$corr_type, fit.model$link, 1e-6, 20)
      if(fit_j$converged==TRUE){sim_ok <- 1}
    }
    table[j,] <- fit_j$mu.coefs
    if(show.step) pbapply::setpb(pb, j)
  }
  
  # return
  return(table)
}
