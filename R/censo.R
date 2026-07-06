# source ulgee
#source("R/ulgee.R")

# prep. 2010
data_2010 <- readxl::read_excel("data/censo/censo_2010.xlsx")
data_2010 <- data_2010[-c(1,29:31),]
names(data_2010) <- c("state", "gini", "pwater", "life", "peletro", "psewage", 
                      "hdi", "income")
data_2010$time <- 3
data_2010$year <- 2010
data_2010$id <- 1:nrow(data_2010)
data_2010$psewage <- data_2010$psewage/100
data_2010$state <- c("AC","AL","AP","AM","BA","CE","DF","ES","GO","MA","MT",
                     "MS","MG","PA","PB","PR","PE","PI","RJ","RN","RS","RO",
                     "RR","SC","SP","SE","TO")

# prep. 2000
data_2000 <- readxl::read_excel("data/censo/censo_2000.xlsx")
data_2000 <- data_2000[-c(1,29:31),]
names(data_2000) <- c("state", "gini", "pwater", "life", "peletro", "psewage",
                      "hdi", "income")
data_2000$time <- 2
data_2000$year <- 2000
data_2000$id <- 1:nrow(data_2000)
data_2000$psewage <- data_2000$psewage/100
data_2000$state <- c("AC","AL","AP","AM","BA","CE","DF","ES","GO","MA","MT",
                     "MS","MG","PA","PB","PR","PE","PI","RJ","RN","RS","RO",
                     "RR","SC","SP","SE","TO")

# prep. 1991
data_1991 <- readxl::read_excel("data/censo/censo_1991.xlsx")
data_1991 <- data_1991[-c(1,29:31),]
names(data_1991) <- c("state", "gini", "pwater", "life", "peletro", "psewage", 
                      "hdi", "income")
data_1991$time <- 1
data_1991$year <- 1991
data_1991$id <- 1:nrow(data_1991)
data_1991$psewage <- data_1991$psewage/100
data_1991$state <- c("AC","AL","AP","AM","BA","CE","DF","ES","GO","MA","MT",
                     "MS","MG","PA","PB","PR","PE","PI","RJ","RN","RS","RO",
                     "RR","SC","SP","SE","TO")

# prep. data
data <- rbind(data_2010, data_2000, data_1991)

# num. clusters
table(table(data$id))

# boxplot 
par(mar=c(5.5,5.5,2,2), mfrow=c(1,2))
colors <- viridis::viridis(3, a=0.65)
boxplot(gamlss.dist::qLO(data$psewage) ~ time, data, ylab="logit(psewage)", 
        xlab="tempo", pch=16, cex.lab=1.5, cex.axis=1.2, ylim=c(-5.3,0.8),
        names=c("1991", "2000", "2010"), col=colors)

# dispersion plot
colors <- viridis::viridis(3)
plot(data$gini, gamlss.dist::qLO(data$psewage), pch=16, 
     col=colors[data$time], ylab="logit(psewage)", xlab="Gini", cex.lab=1.5, 
     cex.axis=1.2, lwd=2, ylim=c(-5.3,0.8))
legend("top", col=colors, legend=c("1991", "2000", "2010"), cex=0.9, bty="n", 
       title="Ano:", horiz=T, lwd=2, inset=c(0,.03))

# add 1991 line
colors <- viridis::viridis(3, a=0.65)
newdata <- data.frame(gini=seq(min(data_1991$gini), 
                               max(data_1991$gini), length.out=100))
newdata$pred <- predict(loess(gamlss.dist::qLO(psewage)~gini, 
                              span=2, data_1991), newdata)
with(newdata, lines(x=gini, y=pred, col=colors[1], lwd=3))

# add 2000 line
newdata <- data.frame(gini=seq(min(data_2000$gini),
                               max(data_2000$gini), length.out=100))
newdata$pred <- predict(loess(gamlss.dist::qLO(psewage)~gini,
                              span=2, data_2000), newdata)
with(newdata, lines(x=gini, y=pred, col=colors[2], lwd=3))

# add 2010 line
newdata <- data.frame(gini=seq(min(data_2010$gini), 
                               max(data_2010$gini), length.out=100))
newdata$pred <- predict(loess(gamlss.dist::qLO(psewage)~gini, 
                              span=2, data_2010), newdata)
with(newdata, lines(x=gini, y=pred, col=colors[3], lwd=3))

# histogram
par(mar=c(5.5,5.5,2,2), mfrow=c(1,3))

# add 1991
set.seed(1546)
colors <- viridis::viridis(3, a=0.65)
p1a <- hist(rul(500, mu=mean(data$psewage[data$year==1991])), plot=F)
p1b <- hist(data$psewage[data$year==1991], plot=F)
plot(p1a, ylim=c(0,15), xlab="psewage", freq=F, main="", cex.axis=1.5, 
     cex.lab=1.8, xlim=c(0,0.5), ylab="Densidade")
plot(p1b, ylim=c(0,15), xlab="psewage", freq=F, main="", cex.axis=1.5, 
     cex.lab=1.8, xlim=c(0,0.5), col=colors[1], add=T)
legend("top", legend=c("   1991", "UL(0.13)"), cex=1.5, bty="n", 
       col=c(colors[1], "gray"), title="Densidades de probabilidade:", 
       horiz=F, lwd=c(2,2), inset=c(0,.03))

# add 2000
set.seed(1546)
p1a <- hist(rul(500, mu=mean(data$psewage[data$year==2000])), plot=F)
p1b <- hist(data$psewage[data$year==2000], plot=F)
plot(p1a, ylim=c(0,15), xlab="psewage", freq=F, main="", cex.axis=1.5, 
     cex.lab=1.8, xlim=c(0,0.5), ylab="Densidade")
plot(p1b, ylim=c(0,15), xlab="psewage", freq=F, main="", cex.axis=1.5, 
     cex.lab=1.8, xlim=c(0,0.5), col=colors[2], add=T)
legend("top", legend=c("   2000", "UL(0.15)"), cex=1.5, bty="n", 
       col=c(colors[2], "gray"), title="Densidades de probabilidade:", 
       horiz=F, lwd=c(2,2), inset=c(0,.03))

# add 2010
set.seed(1546)
p1a <- hist(rul(500, mu=mean(data$psewage[data$year==2010])), plot=F)
p1b <- hist(data$psewage[data$year==2010], plot=F)
plot(p1a, ylim=c(0,15), xlab="psewage", freq=F, main="", cex.axis=1.5, 
     cex.lab=1.8, xlim=c(0,0.5), ylab="Densidade")
plot(p1b, ylim=c(0,15), xlab="psewage", freq=F, main="", cex.axis=1.5, 
     cex.lab=1.8, xlim=c(0,0.5), col=colors[3], add=T)
legend("top", legend=c("   2010", "UL(0.09)"), cex=1.5, bty="n", 
       col=c(colors[3], "gray"), title="Densidades de probabilidade:", 
       horiz=F, lwd=c(2,2), inset=c(0,.03))

# prep. model
y <- data$psewage; time <- data$time; id <- data$id; state <- data$state
X <- cbind(as.numeric(time==1), as.numeric(time==1)*data$gini,
           as.numeric(time==2), as.numeric(time==2)*data$gini, 
           as.numeric(time==2)*data$gini^2,
           as.numeric(time==3), as.numeric(time==3)*data$gini)

# model
fit_1 <- ulgee(y, X, time, id, "EXC", "logit", 1e-6, 40)

# est. 
round(fit_1$mu.coefs, 2)
round(diag(fit_1$vcov), 2)
round(fit_1$pvalues, 4)
round(fit_1$rho, 2)
round(fit_1$qic, 2)

# qic comparison
table <- matrix(NA, 4, 1)
corr_type <- c("EXC", "AR1", "UNS", "IDE")
rownames(table) <- corr_type; colnames(table) <- "qic"
for(i in 1:length(corr_type)){
  fit_i <- ulgee(y, X, time, id, corr_type[i], "logit", 1e-6, 50, show.step=F)
  table[i,1] <- fit_i$qic
}
round(table,2)

# diagnostic 
set.seed(3489)
diag_quant(fit_1, X, 100, n=1, label.id=data$state, label.time=data$time)
diag_quant(fit_1, X, 100, T, T, n=0, label.id=data$state, label.time=data$time)

# sensitivity 
sens_conf(fit_1, 4, 1, 4, 1, label.id=data$state, label.time=data$time)

# show points highlighted
par(mar=c(5.5,5.5,2,2), mfrow=c(1,1))
colors <- viridis::viridis(3, a=1)
plot(data$gini, gamlss.dist::qLO(data$psewage), pch=16, col=colors[data$time],
     ylab="logit(psewage)", xlab="Gini", cex.lab=1.5, cex.axis=1.2, 
     lwd=2, ylim=c(-5.3,0.8))
legend("top", col=colors, legend=c("1991", "2000", "2010"), cex=1.1, bty="n", 
       title="Ano:", horiz=T, lwd=2, inset=c(0,.03))
identify(x=data$gini, y=gamlss.dist::qLO(data$psewage), 
         label=paste(paste("(",paste(data$state, data$time, sep=","), 
                           sep=""),")",sep=""), n=3, cex=1.2)

# inference remove points
fit_i <- ulgee(y[!(id==24&time==2)], X[!(id==24&time==2),], 
               time[!(id==24&time==2)], id[!(id==24&time==2)], 
               "EXC", "logit", 1e-6, 40)

# est. 
round(fit_1$mu.coefs, 2)
round(diag(fit_1$vcov), 2)
round(fit_1$pvalues, 4)
round(fit_1$rho, 2)
round(fit_1$qic, 2)

# prep position 
pos <- as.numeric(data$id==3 & data$time==1)
pos <- pos + as.numeric(data$id==24 & data$time==2)
pos <- pos + as.numeric(data$id==10 & data$time==1)
pos <- pos==1

# sensitivity dropp
set.seed(9792)
table <- sens_mrc(fit_1, y, X, 10 , pos)
print(xtable::xtable(table, digits=2), type="latex",
      include.rownames=F, include.colnames=F)

# asymptotic check
set.seed(1111)
table <- sens_coef(fit_1, X, 100)
titles <- c(expression(paste("Quantis empíricos de ",hat(italic("\u03b1"))[1])),
            expression(paste("Quantis empíricos de ",hat(italic("\u03b2"))[1])),
            expression(paste("Quantis empíricos de ",hat(italic("\u03b1"))[2])),
            expression(paste("Quantis empíricos de ",hat(italic("\u03b2"))[2])),
            expression(paste("Quantis empíricos de ",hat(italic("\u03c4"))[2])),
            expression(paste("Quantis empíricos de ",hat(italic("\u03b1"))[3])),
            expression(paste("Quantis empíricos de ",hat(italic("\u03b2"))[3])))
par(mar=c(5.5,5.5,2,2), mfrow=c(1,3))
for(i in 1:7){
  if(i==7) plot(1, type="n", xaxt="n", yaxt="n", xlab="", ylab="", bty="n",
                xlim=c(-3,3), ylim=c(-3,3))
  pcoef_i <- (table[,i]-mean(table[,i]))/sd(table[,i])
  qqnorm(pcoef_i, xlab="Quantis teóricos", ylab=titles[i],
         pch=16, cex.lab=1.8, cex.axis=1.5, main="", xlim=c(-3,3), ylim=c(-3,3))
  abline(a=0, b=1, xlab="", ylab="", lty=2, lwd=1)
}

