##### HW 3, Problem 6
sig = sqrt(4)
n = 500

### Generating Data
set.seed(1)
ar = 0.95
#ar = 0.30
x1 = arima.sim(list(ar=ar),n=n, sd=sig)

### Certain portions of the data vector
yn = x1[-1]
yn_1 = x1[-n]

### Sampling for sigma^2
B = 10000
alpha = (n-2)/2
beta = (sum(yn^2)-sum(yn*yn_1)^2/sum(yn_1^2))/2
po.sig = 1/rgamma(B, alpha, beta)

### Sampling for rho
mu = sum(yn*yn_1)/sum(yn_1^2)
si = po.sig/sqrt(sum(yn_1^2))
po.ro = rnorm(B, mu, si)

### Plots
#pdf("./fig_1.pdf", width = 8, height = 8)
par(mfrow=c(1,1), mar = c(4.1, 4.1, 3.1, 1.1))
layout(matrix(c(1,2,1,3),2,2))
plot(x1, type='o', pch=20, main="Simulated data", cex.main = 1.5, ylab = expression(y[t]))
plot(density(po.ro), main = expression("Posterior for" ~ rho), xlab ="", ylab="", cex.main = 1.5, lwd = 2)
abline(v = mean(po.ro), col = 'red', lwd = 2)
plot(density(po.sig), main = expression("Posterior for" ~ sigma^2), xlab ="", ylab="", cex.main = 1.5, lwd = 2)
abline(v = mean(po.sig), col = 'red', lwd = 2)
#dev.off()

### Summary statistics
c("mean"=mean(po.ro), "var"=var(po.ro), quantile(po.ro, c(0, 0.025, 0.5, 0.975, 1)))
c("mean"=mean(po.sig), "var"=var(po.sig), quantile(po.sig, c(0, 0.025, 0.5, 0.975, 1)))

### Scatterplot for the joint distribution
# par(mfrow=c(1,1), mar = c(5.1,4.1,4.1,2.1))
# plot(po.ro, po.sig, pch = 20, xlab="", ylab="", cex=0.5)
# mtext(expression(rho), side=1, line=3, cex=1.5)
# mtext(expression(sigma^2), side=2, line=2, cex=1.5)
