
# pacotes
library(MASS)

# prep
n <- 1000                             
mu <- c(10, 5, 1)                                   
S <- matrix(c(1, 0.8, 0.1, 0.8, 1, 0.5, 0.1, 0.5, 1), ncol=3)
S <- matrix(c(1, 0.7, 0.7, 0.7, 1, 0.7, 0.7, 0.7, 1), ncol=3)

# sim
y <- mvrnorm(n=n, mu=mu, Sigma=S)

# test1
cor(y[,1], y[,2], method="pearson")
cor(y[,1], y[,2]^2, method="pearson")
cor(y[,1]^2, y[,2], method="pearson")
cor(y[,1]^2, y[,2]^2, method="pearson")

# prep
temp1 <- cor(y[,1], y[,2], method="pearson")
temp2 <- cor(y[,1], y[,3], method="pearson")
temp3 <- cor(y[,2], y[,3], method="pearson")

# test2
cor(y[,1], y[,2]*y[,3], method="pearson")
(temp1+temp2)*temp3/2
cor(y[,1], y[,2]*y[,3]^2, method="pearson")
(temp1+2*temp2)*temp3/3
cor(y[,1], y[,2]^2*y[,3], method="pearson")
(2*temp1+temp2)*temp3/3
