source("./sampling.R")
library(mwBASE)
#library(KMsurv)
library(survival)

trailing = function(x, digits = 4)
    formatC(x, digits=digits, format="f")

dat = as.matrix(read.table("./kidney-dat.txt", header = TRUE))
dat[,4] = dat[,4] - 1   # male 0, female 1
#dat[,4] = 2 - dat[,4]   # male 1, female 0


### Frequentist model
coxph(Surv(time, nu) ~ as.factor(Sex) + Age + frailty(cluster, theta = 1), data = data.frame(dat))


### Bayesian model
set.seed(1)
mod = get_samples(nburn =  800000, nmcmc = 200000, window = 5000, k = 2)


# colMeans(mod$accept)
# par(mfrow = c(2,3), mar = c(0,0,0,0))
# for (i in 1:5)
#     plot(mod$params[,i], type='l', axes=FALSE)
# par(mfrow = c(1,1), mar = c(5.1,4.1,4.1,2.1))
# 
# pairs(mod$params[,1:5], pch = 16, col = rgb(seq(0, 1, length = 5000), 0, 0))
# pairs(mod$params[,6:10], pch = 16, col = rgb(seq(0, 1, length = 5000), 0, 0))
# pairs(mod$params[,11:15], pch = 16, col = rgb(seq(0, 1, length = 5000), 0, 0))


pdf("./figs/posterior.pdf", width = 12, height = 12)
par(mfrow = c(2,3), oma = c(0,0,2,0))
plot_hpd(mod$params[,1], col1 = 'dodgerblue', bty='n', xlab ="",
    main = expression(kappa), cex.main = 2, cex.lab = 1.3)
plot_hpd(mod$params[,2], col1 = 'dodgerblue', bty='n', xlab ="",
    main = expression(beta[1]), cex.main = 2, cex.lab = 1.3)
plot_hpd(mod$params[,3], col1 = 'dodgerblue', bty='n', xlab ="",
    main = expression(beta[2]), cex.main = 2, cex.lab = 1.3)
plot_hpd(mod$params[,4], col1 = 'dodgerblue', bty='n', xlab ="",
    main = expression(gamma), cex.main = 2, cex.lab = 1.3)
plot_hpd(mod$params[,5], col1 = 'dodgerblue', bty='n', xlab ="",
    main = expression(alpha), cex.main = 2, cex.lab = 1.3)
# plot(0, type='n', axes=FALSE, xlab="", ylab="", ylim = c(-2, 1.1))
# text(c(0.7, 0.9, 1.1, 1.3), rep(1.1, 4), c("Mean", "SD", "2.5%", "97.5%"), cex = 1.8)
# text(0.7,  rev(seq(-.2, 0.8, length = 5)), trailing(colMeans(mod$params[,1:5]), 3), cex = 1.5)
# text(0.9,  rev(seq(-.2, 0.8, length = 5)), trailing(apply(mod$params[,1:5], 2, sd), 3), cex = 1.5)
# text(1.1,  rev(seq(-.2, 0.8, length = 5)), trailing(apply(mod$params[,1:5], 2, quantile, 0.025), 3), cex = 1.5)
# text(1.3,  rev(seq(-.2, 0.8, length = 5)), trailing(apply(mod$params[,1:5], 2, quantile, 0.975), 3), cex = 1.5)
# mtext(expression(kappa), 2, at = seq(-.2, 0.8, length = 5)[5], las = 1, line = 1)
# mtext(expression(beta[1]), 2, at = seq(-.2, 0.8, length = 5)[4], las = 1, line = 1)
# mtext(expression(beta[2]), 2, at = seq(-.2, 0.8, length = 5)[3], las = 1, line = 1)
# mtext(expression(gamma), 2, at = seq(-.2, 0.8, length = 5)[2], las = 1, line = 1)
# mtext(expression(alpha), 2, at = seq(-.2, 0.8, length = 5)[1], las = 1, line = 1)
# title("Posteriors", outer = TRUE, cex.main = 2)
par(mfrow = c(1,1), mar = c(5.1,4.1,4.1,2.1), oma = c(0,0,0,0))
dev.off()


# # Survival for person 1, cluster 1
# tt = seq(0, max(dat[,2]), length = 50)
# vv = apply(mod$params, 1, function(p)
#     exp(-p[4]*tt^p[5]*p[6]*exp(dat[1,4:5] %*% p[2:3])))
# mm = apply(vv, 1, mean)
# qq = apply(vv, 1, quantile, c(0.025, 0.975))
# 
# plot(tt,  mm, type='l', lwd = 3)
# lines(tt, qq[1,])
# lines(tt, qq[2,])


pdf("./figs/frailty.pdf", width = 12, height = 12)
boxplot(mod$params[seq(1, NROW(mod$params), by = 50),-(1:5)], axes = FALSE, horizontal = TRUE,
    main = "Posterior frailties", xlab = expression(w[i]), at = 38:1,
    ylab = "Observation", col = 'dodgerblue', cex.main = 2, cex.lab = 1.5,
    xlim = c(1, 38), ylim = c(0, max(mod$params[,-(1:5)])))
axis(1)
axis(2, at = unique(c(seq(1, 38, by = 3), 38)), labels = unique(c(seq(38, 1, by = -3), 1)))
text(c(12.5, 14), c(39, 39), c("Mean", "SD"), cex = 1.5)
text(12.5, 38:1, trailing(colMeans(mod$params[,-(1:5)]), 3))
text(14, 38:1, trailing(apply(mod$params[,-(1:5)], 2, sd), 3))
dev.off()


# Stuff for Latex
tmp = cbind(colMeans(mod$params[,1:5]),
    apply(mod$params[,1:5], 2, sd),
    apply(mod$params[,1:5], 2, quantile, 0.025),
    apply(mod$params[,1:5], 2, quantile, 0.975))
colnames(tmp) = c("mean", "sd", "0.025%", "0.975%")
rownames(tmp) = c("kappa", "beta0", "beta1", "gamma", "alpha")

mean(exp(mod$params[,2]))
