#load("./workspaces/toy_data.RData")
#
#load("./workspaces/toy_para.RData")
#load("./workspaces/toy_cpo_para.RData")
#load("./workspaces/toy_para_preds.RData")
#
#para.theta = thin.theta
#para.mu = thin.mu
#para.sigma = thin.sigma
#para.tau = thin.tau
#
#para.t0 = pred.theta_0
#para.y0 = pred.y_0
#para.qhpd = qhpd
#
#
#save(para.theta, para.mu, para.sigma, para.tau, para.t0, para.y0,
#    para.qhpd, para.cpo, pred.x, x, y, groups,
#    file = "./workspaces/toy_para_everything.RData")
#
#load("./workspaces/toy_semi.RData")
#load("./workspaces/toy_cpo_semi.RData")
#load("./workspaces/toy_dpm_preds.RData")
#
#dpm.theta = thin.theta
#dpm.alpha = thin.alpha
#dpm.mu = thin.mu
#dpm.sigma = thin.sigma
#dpm.tau = thin.tau
#
#dpm.t0 = pred.theta_0
#dpm.y0 = pred.y_0
#dpm.qhpd = qhpd
#
#save(dpm.theta, dpm.alpha, dpm.mu, dpm.sigma, dpm.tau, dpm.t0, dpm.y0,
#    dpm.qhpd, dpm.cpo, x, x, y, groups,
#    file = "./workspaces/toy_dpm_everything.RData")
#
#load("./workspaces/toy_groups.RData")

load("./workspaces/toy_para_everything.RData")
load("./workspaces/toy_dpm_everything.RData")

# Red: linear         Green: quadratic
cols = ifelse(groups == 1, "red", "green")

para.nmcmc = NROW(para.theta)
dpm.nmcmc = NROW(dpm.theta)


pdf("./figs/toy_data.pdf", width = 8, height = 8)
matplot(x, t(y), type='l', col = "black", lty = 1, lwd = 0.5, cex.lab = 1.5, ylab = "y",
    main = "Simulated 'toy' data, n = 100", cex.main = 2)
dev.off()


# Black: parametric,     Blue: DPM

### Pairs plot for mu
#pairs(dpm.mu, pch = 20, cex = 0.5)
#pairs(para.mu, pch = 20, cex = 0.5)
pdf("./figs/toy_mu.pdf", width = 8, height = 8)
pairs(rbind(dpm.mu, para.mu), pch = 20, cex = 0.5, col = c(rep(rgb(0,0,1,0.2), para.nmcmc),
    rep(rgb(0,0,0,0.2), dpm.nmcmc)), labels = expression(mu[0], mu[1], mu[2]),
    cex.labels = 5, main = expression("Posterior for" ~ mu),
    cex.main = 2)
dev.off()


### Densities for log(det(Sigma))
#par(mfrow = c(2,1), mar = c(4.1, 4.1, 3.1, 2.1))
#plot(sapply(para.sigma, function(x) determinant(x)$modulus[1]), type='l', ylab = "log(det(Sigma)",
#    main = expression("Parametric -- Posterior for" ~ Sigma), cex.main = 2)
#plot(sapply(dpm.sigma, function(x) determinant(x)$modulus[1]), type='l', ylab = "log(det(Sigma)",
#    main = expression("DPM -- Posterior for" ~ Sigma), col = 'blue', cex.main = 2)
#par(mfrow = c(1,1), mar = c(5.1, 4.1, 4.1, 2.1))

pdf("./figs/toy_sigma.pdf", width = 8, height = 8)
plot(density(sapply(para.sigma, function(x) determinant(x)$modulus[1])), xlim = c(-13, -4),
    main = expression("Posterior for" ~ Sigma), cex.main = 2, xlab = "log(det(Sigma)", lwd = 2)
lines(density(sapply(dpm.sigma, function(x) determinant(x)$modulus[1])), col = 'blue', lwd = 2)
dev.off()

#par(mfrow = c(2,1), mar = c(4.1, 4.1, 3.1, 2.1))
#plot(para.tau, type='l', ylab = "tau",
#    main = expression("Parametric -- Posterior for" ~ tau), cex.main = 2)
#plot(dpm.tau, type='l', ylab = "tau",
#    main = expression("DPM -- Posterior for" ~ tau), col = 'blue', cex.main = 2)
#par(mfrow = c(1,1), mar = c(5.1, 4.1, 4.1, 2.1))


### densities for tau
pdf("./figs/toy_tau.pdf", width = 8, height = 8)
plot(density(para.tau), xlim = c(0.004, 0.012), yli = c(0, 1200),
    main = expression("Posterior for" ~ tau), cex.main = 2, xlab = "tau", lwd = 2)
lines(density(dpm.tau), col = 'blue', lwd = 2)
dev.off()


### density for alpha
pdf("./figs/toy_alpha.pdf", width = 8, height = 8)
plot(density(dpm.alpha), main = expression("Posterior for" ~ alpha),
    cex.main = 2, xlab = "alpha", lwd = 2, col = 'blue')
dev.off()

### summary for dpm.alpha
mean(dpm.alpha)
quantile(dpm.alpha, c(0.025, 0.975))


### cpo plot
pdf("./figs/toy_cpo.pdf", width = 8, height = 8)
plot(log(dpm.cpo) - log(para.cpo), pch = 20, main = "log(CPO.DPM / CPO.Parametric)",
    cex.main = 2)
abline(h = 0, lty = 2)
legend("topleft", box.lty = 0, legend = paste0(mean(dpm.cpo > para.cpo), "% favor the DPM"),
    cex = 2)
dev.off()


### posterior predictions for a new observation y_0 based a new theta_0
### draws for theta_0
pdf("./figs/toy_theta_0.pdf", width = 8, height = 8)
pairs(rbind(dpm.t0, para.t0), pch = 20, cex = 0.5, col = c(rep(rgb(0,0,1,0.2), para.nmcmc),
    rep(rgb(0,0,0,0.2), dpm.nmcmc)), labels = expression(beta[0], beta[1], beta[2]),
    cex.labels = 5, main = expression("Predictions for a new" ~ theta ~ "= (" ~ beta[0] ~~
    beta[1] ~~ beta[2] ~ ")"), cex.main = 2)
dev.off()


### draws for a new y_0 (dpm)
pdf("./figs/toy_y_0.pdf", width = 16, height = 8)
par(mfrow = c(1,2), mar = c(4.1, 4.1, 3.1, 2.1))
plot(0, type='n', xlim = range(x), ylim = c(.5, 7.5), cex.lab = 1.5,
    main = "DP mixture model posterior predictions", cex.main = 2, xlab = "x", ylab = "y")
matplot(x, t(y), type='l', lty = 1, lwd = 0.5, add = TRUE, col = 'gray20')
for (i in 1:4)
    lines(pred.x, sapply(dpm.qhpd, function(x) x[i]), col = 'darkblue', lwd = 3)
legend("topleft", box.lty = 0, legend = c("Simulated data", "95% Prediction bounds"),
    col = c("gray20", "darkblue"), lwd = c(0.5, 3), cex = 2)


### y_0 (parametric)
plot(0, type='n', xlim = range(x), ylim = c(.5, 7.5), cex.lab = 1.5,
    main = "Parametric model posterior predictions", cex.main = 2, xlab = "x", ylab = "y")
matplot(x, t(y), type='l', lty = 1, lwd = 0.5, add = TRUE, col = 'gray20')
lines(pred.x, para.qhpd[1,], col = 'darkblue', lwd = 3)
lines(pred.x, para.qhpd[2,], col = 'darkblue', lwd = 3)
legend("topleft", box.lty = 0, legend = c("Simulated data", "95% Prediction bounds"),
    col = c("gray20", "darkblue"), lwd = c(0.5, 3), cex = 2)
par(mfrow = c(1,1), mar = c(5.1, 4.1, 4.1, 2.1))
dev.off()

