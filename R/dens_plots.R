# source ulgee
source("R/ulgee.R")

# density function
par(mar=c(5.5,5.5,2,2), mfrow=c(1,2))
mu <- c(0.2, 0.4, 0.6, 0.8)
dul_y <- function(y) dul(y, mu=0.2)
curve(dul_y, lty=1, ylab="Densidade", xlab=expression(italic(y)), cex.lab=1.5, 
      cex.axis=1.2, lwd=2, ylim=c(0,6.5), xlim=c(0.0001,1))
for(i in 2:length(mu)){
  dul_y <- function(y) dul(y, mu=mu[i])
  curve(dul_y, add=TRUE, lty=i, lwd=2)
}
legend("top", lty=1:length(mu), legend=mu, cex=0.9, bty="n", 
       title=expression(paste("Valores de ", italic("\u03bc"), ":")), 
       horiz=T, inset=c(0,.03))

# variance
curve(ul.var, lty=1, ylab="Variância", xlab=expression(italic("\u03bc")), 
      cex.lab=1.5, cex.axis=1.2, lwd=2, from=0.01, to=0.99)
abline(v=0.5, lty=2, lwd=2)
