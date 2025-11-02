
# confidence band
par(mar=c(5.5,5.5,2,2), mfrow=c(1,3))
colors <- viridis::viridis(3, a=0.65)

# 1991
z <- seq(0.48, 0.68, 0.01)
z <- cbind(1, z, 0, 0, 0, 0, 0); c <- qchisq(0.95, 7)
var.z <- z%*%fit_1$vcov%*%t(z)
sd.z <- sqrt(c*diag(var.z))
eta <- z%*%fit_1$mu.coefs
eta.s <- eta + sd.z
eta.i <- eta - sd.z
prob <- exp(eta)/(1+exp(eta))
prob.s <- exp(eta.s)/(1+exp(eta.s))
prob.i <- exp(eta.i)/(1+exp(eta.i))
plot(z[,2], prob, type="n", ylim=c(0,0.7), xlab="Gini", ylab="Proporção ajustada", 
     cex.axis=1.5, cex.lab=1.8)
polygon(c(z[,2], rev(z[,2])), c(prob.i, rev(prob.s)), col=colors[1],
        border=NA)
lines(z[,2], prob, lty=2,lwd=2)
lines(z[,2], prob.s,lwd=2)
lines(z[,2], prob.i,lwd=2)
legend("top", col=colors[1], legend=c("1991"), cex=1.5, bty="n", 
       title="Ano:", horiz=T, lwd=2, inset=c(0,.03))
#points(data$gini[time==1], y[time==1], pch=16, cex=1.5)

# 2000
z <- seq(0.48, 0.68, 0.01)
z <- cbind(0, 0, 1, z, z^2, 0, 0); c <- qchisq(0.95, 7)
var.z <- z%*%fit_1$vcov%*%t(z)
sd.z <- sqrt(c*diag(var.z))
eta <- z%*%fit_1$mu.coefs
eta.s <- eta + sd.z
eta.i <- eta - sd.z
prob <- exp(eta)/(1+exp(eta))
prob.s <- exp(eta.s)/(1+exp(eta.s))
prob.i <- exp(eta.i)/(1+exp(eta.i))
plot(z[,4], prob, type="n", ylim=c(0,0.7), xlab="Gini", ylab="Proporção ajustada", 
     cex.axis=1.5, cex.lab=1.8)
polygon(c(z[,4], rev(z[,4])), c(prob.i, rev(prob.s)), col=colors[2],
        border=NA)
lines(z[,4], prob, lty=2,lwd=2)
lines(z[,4], prob.s,lwd=2)
lines(z[,4], prob.i,lwd=2)
legend("top", col=colors[2], legend=c("2000"), cex=1.5, bty="n", 
       title="Ano:", horiz=T, lwd=2, inset=c(0,.03))
#points(data$gini[time==2], y[time==2], pch=16, cex=1.5)

# 2010
z <- seq(0.48, 0.68, 0.01)
z <- cbind(0, 0, 0, 0, 0, 1, z); c <- qchisq(0.95, 7)
var.z <- z%*%fit_1$vcov%*%t(z)
sd.z <- sqrt(c*diag(var.z))
eta <- z%*%fit_1$mu.coefs
eta.s <- eta + sd.z
eta.i <- eta - sd.z
prob <- exp(eta)/(1+exp(eta))
prob.s <- exp(eta.s)/(1+exp(eta.s))
prob.i <- exp(eta.i)/(1+exp(eta.i))
plot(z[,7], prob, type="n", ylim=c(0,0.7), xlab="Gini", ylab="Proporção ajustada",
     cex.axis=1.5, cex.lab=1.8)
polygon(c(z[,7], rev(z[,7])), c(prob.i, rev(prob.s)), col=colors[3], 
        border=NA)
lines(z[,7], prob, lty=2,lwd=2)
lines(z[,7], prob.s,lwd=2)
lines(z[,7], prob.i,lwd=2)
legend("top", col=colors[3], legend=c("2010"), cex=1.5, bty="n", 
       title="Ano:", horiz=T, lwd=2, inset=c(0,.03))
#points(data$gini[time==3], y[time==3], pch=16, cex=1.5)
