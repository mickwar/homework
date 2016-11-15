source("./sampling.R")
library(mwBASE)
# library(KMsurv)
# library(survival)

dat = as.matrix(read.table("./kidney-dat.txt", header = TRUE))
#dat[,4] = 2 - dat[,4]

set.seed(3)
out = get_samples(nburn =  800000, nmcmc = 200000, window = 1000, k = 2)

colMeans(out$accept)
par(mfrow = c(2,3), mar = c(0,0,0,0))
for (i in 1:5)
    plot(out$params[,i], type='l', axes=FALSE)
par(mfrow = c(1,1), mar = c(5.1,4.1,4.1,2.1))

pairs(out$params[,1:5], pch = 16, col = rgb(seq(0, 1, length = 5000), 0, 0))
pairs(out$params[,6:10], pch = 16, col = rgb(seq(0, 1, length = 5000), 0, 0))
pairs(out$params[,11:15], pch = 16, col = rgb(seq(0, 1, length = 5000), 0, 0))

par(mfrow = c(2,3))
plot_hpd(out$params[,1], col1 = 'dodgerblue', bty='n')
for (i in 2:5)
    plot_hpd(out$params[,i], col1 = 'dodgerblue', bty='n')
par(mfrow = c(1,1), mar = c(5.1,4.1,4.1,2.1))

# eta, beta0, beta1, gamma, alpha
mean(out$params[,1])  # kappa
mean(out$params[,2])  # beta1
mean(out$params[,3])  # beta2
mean(out$params[,4])  # gamma
mean(out$params[,5])  # alpha

plot(density(apply(out$params[,6:43], 2, mean)))

# Survival for person 1, cluster 1
tt = seq(0, max(dat[,2]), length = 50)
vv = apply(out$params, 1, function(p)
    exp(-p[4]*tt^p[5]*p[6]*exp(dat[1,4:5] %*% p[2:3])))
mm = apply(vv, 1, mean)
qq = apply(vv, 1, quantile, c(0.025, 0.975))

plot(tt,  mm, type='l', lwd = 3)
lines(tt, qq[1,])
lines(tt, qq[2,])
