
# source ulgee
#source("R/ulgee.R")

# data generate
set.seed(4379)
n <- 500
s <- 4
id <- rep(1:n, each=s)
time <- rep(1:s, n)
beta <- c(-3, 6)
rho_sim <- 0.5
corr_type <- "AR1"
x <- runif(n*s, 0, 1)  
data <- cbind(1, x, time, id)

# y simulation (logit link)
eta <- data[,-c(3,4)]%*%beta; mu <- gamlss.dist::pLO(eta)
y_sim <- rul_corr(mu, rho_sim, corr_type, data[,3], data[,4])

# model
fit_1 <- ulgee(y_sim, data[,-c(3,4)], data[,3], data[,4], "AR1", 
               "logit", 1e-06, 20)

# coefficients 
round(fit_1$mu.coefs, 2)
round(diag(fit_1$vcov), 2)
round(fit_1$pvalues, 4)
round(fit_1$rho, 2)

# diagnostic
diag_quant(fit_1, data[,-c(3,4)], 100, T, T, n=0)

# correlation
r <- (y_sim-mu)/sqrt(ul.var(mu)); rho_hat <- 0
ids <- as.numeric(names(table(id))[table(id)>=2])
n <- 0
for(i in ids){
  pos <- id==i; time_i <- time[pos]
  r_i <- numeric(max(time)-min(time)+1)
  r_i[time_i] <- r[pos]
  aux <- outer(r_i, r_i, FUN="*")
  n <- n + 1*(aux!=0)
  rho_hat <- rho_hat + aux
}
rho_hat <- rho_hat/(n)
cov2cor(rho_hat)
