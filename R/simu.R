
# source ulgee
#source("R/ulgee.R")

# prep. simulation
beta <- c(-3, 6)
R <- 5000
rho <- c(-0.1, 0.3, 0.7)
n <- c(500, 50, 10)
s <- c(10, 5, 3)

# AR1
set.seed(4598)
table_1 <- ul_sim(beta, R, rho, n, s, "AR1", "EXC")
print(xtable::xtable(table_1[,-2], digits=4), type="latex", 
      include.rownames=F, include.colnames=F)

# EXC
set.seed(8878)
table_2 <- ul_sim(beta, R, rho, n, s, "EXC", "AR1")
print(xtable::xtable(table_2[,-2], digits=4), type="latex", 
      include.rownames=F, include.colnames=F)
