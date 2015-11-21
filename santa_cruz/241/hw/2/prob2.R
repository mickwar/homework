library(DPpackage)
source("~/files/R/mcmc/bayes_functions.R")

y = as.numeric(unlist(read.table("./data.txt")))

n = length(y)
ord = order(y)

par(mfrow = c(1,1), mar = c(5.1, 4.1, 4.1, 2.1))
plot(density(y))
curve(0.2*dnorm(x, -5, 1) + 0.5*dnorm(x, 0, 1) + 0.3*dnorm(x, 3.5, 1), add = TRUE, col = 'red')
legend("topleft", legend = c("data", "truth"), col = c("black", "red"), lwd = 1, box.lwd = 0, cex = 1.5)

nburn = 2000
nmcmc = 7000


prior1 = list(m2 = 0, s2 = sqrt(3),    # prior mean and sd for mean of normal part of G0
              tau1 = 2, tau2 = 2,      # twice the shape and rate for kappa
              nu1 = 3,                 # fixed c, shape in inverse gamma part of G0
              nu2 = 5, psiinv2 = 1/10, # rate in inverse gamma part of G0 (an inverse gamma itself)
              a0 = 1, b0 = 1)          # priors on alpha (same as prob 1)
mcmc1 = list(nburn = nburn, nsave = nmcmc, ndisplay = 1)
                

fit1 = DPdensity(y = y, prior = prior1, mcmc = mcmc1, state = NULL, status = TRUE)

str(fit1$save.state)

plot((fit1$save.state$thetasave[,1]), main="mean of normal part", type='l')
plot((fit1$save.state$thetasave[,2]), main=expression(kappa), type='l')
plot((fit1$save.state$thetasave[,3]), main="rate of ig part", type='l')
plot((fit1$save.state$thetasave[,4]), type='l')
plot((fit1$save.state$thetasave[,5]), main=expression(alpha), type='l')

plot(density(fit1$save.state$thetasave[,1]), main="mean of normal part")
plot(density(fit1$save.state$thetasave[,2]), main=expression(kappa))
plot(density(fit1$save.state$thetasave[,3]), main="rate of ig part")
plot(table(fit1$save.state$thetasave[,4])/nmcmc)
plot(density(fit1$save.state$thetasave[,5]), main=expression(alpha))


plot(density(fit1$save.state$randsave[,501]))
plot(density(fit1$save.state$randsave[,502]))
fit1$save.state$randsave[,seq(1, 2*n, by = 2)]
fit1$save.state$randsave[,seq(2, 2*n, by = 2)], 2, mean)

z1 = apply(fit1$save.state$randsave[,c(2*245-1,2*245+0)], 1,
    function(x) length(unique(x)))
z2 = apply(fit1$save.state$randsave[,c(2*140-1,2*140+0)], 1,
    function(x) length(unique(x)))
z3 = apply(fit1$save.state$randsave[,c(2*245-1,2*245+0,2*140-1,2*140+0)], 1,
    function(x) length(unique(x)))
table(z1); table(z2); table(z3)

### Plot of posterior means for location and scale of eath observation
cluster_loc = apply(fit1$save.state$randsave[,seq(1, 2*n, by = 2)], 2, mean)
cluster_sca = apply(fit1$save.state$randsave[,seq(2, 2*n, by = 2)], 2, mean)
plot(cluster_loc, cluster_sca, pch = 20, cex = 0.5)

plot(0, type = 'n', xlim = range(fit1$save.state$randsave[,seq(1, 2*n, by = 2)]),
    ylim = range(fit1$save.state$randsave[,seq(2, 2*n, by = 2)]))
points(cluster_loc, cluster_sca, pch = 20, cex = 0.6, col = 'red')
    

### Predictions on mixing parameters
theta = fit1$save.state$randsave[,c(501, 502)]
d = kde2d(theta[,1], theta[,2], n = 200)
contour(d, xlim = c(-7, 7), ylim = c(0, 3), nlevels = 100,
    xlab = "location", ylab = "scale", main = "clusters")
points(theta, pch = 20, cex = 0.5)

points(fit1$save.state$randsave[,c(501, 502)], pch = 20, cex = 0.5,
    xlab = "location", ylab = "scale", main = "clusters")

hpd.theta = apply(theta, 2, function(x) hpd.mult(x, density(x)))
par(mfrow = c(2,1), mar = c(3.1, 2.1, 2.1, 1.1))
hpd.plot(density(theta[,1]), hpd.theta[[1]], main = "new location")
hpd.plot(density(theta[,2]), hpd.theta[[2]], main = "new scale")

### Predictions on new observation
y0 = fit1$save.state$randsave[,503]
hpd.y = hpd.mult(y0, density(y0))
par(mfrow = c(1,1), mar = c(3.1, 2.1, 2.1, 1.1))
hpd.plot(density(y0), hpd.y, lwd = 3, main=expression("Posterior predictions of new" ~ y[0]))
lines(density(y), lwd = 3)
legend("topleft", box.lty = 0, legend = "Data", lwd = 3, cex = 1.5)

