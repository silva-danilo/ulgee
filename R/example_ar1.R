
# source ulgee
#source("R/ulgee.R")

# data AR1 generate
set.seed(4379)
n <- 500
s <- 4
id <- rep(1:n, each=s)
time <- rep(1:s, n)
beta <- c(-3, 6)
rho_sim <- 0.7
corr_type <- "AR1"
x <- runif(n*s, 0, 1)  
data <- cbind(1, x, time, id)
data <- data[-c(1:3, 6:7, 11:14),]

# y simulation (logit link)
eta <- data[,-c(3,4)]%*%beta; mu <- gamlss.dist::pLO(eta)
y_sim <- rul_corr(mu, rho_sim, corr_type, data[,3], data[,4])
rm(n, s, time, id, beta, rho_sim, x, corr_type, mu)

# model
fit_1 <- ulgee(y_sim, data[,-c(3,4)], data[,3], data[,4], "AR1", 
               "logit", 1e-06, 20)

# coefficients 
round(fit_1$mu.coefs, 2)
round(diag(fit_1$vcov), 2)
round(fit_1$pvalues, 4)
round(fit_1$rho, 2)

# diagnostic
diag_quant(fit_1, data[,-c(3,4)], 1, n=0, label.id=data[,3], 
           label.time=data[,4])
diag_quant(fit_1, data[,-c(3,4)], 100, T, T, n=0, label.id=data[,3], 
           label.time=data[,4])
