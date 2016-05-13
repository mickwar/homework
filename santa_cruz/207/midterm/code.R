source("~/files/R/mcmc/bayes_functions.R")
dat = read.table("~/files/data/volcano.txt", header=TRUE)
y = log(dat[1:62,3])
n = length(y)

pdf("./figs/data.pdf", width = 6, height = 6)
#plot(density(y), lwd = 3, xlab = "Log Interevent Time", main = "Etna volcano data", cex.lab = 1.3)
hist(y, xlab = "Log Interevent Time", main = "Etna volcano data", cex.lab = 1.3,
    col = 'gray50', border = 'white', freq = FALSE, breaks = 20)
dev.off()

set.seed(1)
### Model 1
# Hyperpriors
m = 6   # mu mean
s2 = 10 # mu variance
a = 3   # sig2 alpha
b = 4   # sig2 beta
c = 3   # tau2 alpha
d = 4   # tau2 beta
df = 5  # df for the t, fixed

nburn = 10000
nmcmc = 20000
window = 500

params.lambda.i = matrix(0, nburn + nmcmc, n)
params.mu.i = matrix(0, nburn + nmcmc, n)
params.mu = double(nburn + nmcmc)
params.tau2 = double(nburn + nmcmc)
params.sig2 = double(nburn + nmcmc)

params.tau2[1] = 1
params.sig2[1] = 1

for (i in 2:(nburn + nmcmc)){
    if (floor(i / window) == i/window)
        cat(i, "/", nburn+nmcmc, "\r")
    # Update lambda's
    params.lambda.i[i,] = rgamma(n, (df + 1)/2, df/2 + 1/(2*params.sig2[i-1])*
        (y - params.mu.i[i-1,])^2)

    # Update mu.i's
    params.mu.i[i,] = rnorm(n,
        (params.tau2[i-1]*y*params.lambda.i[i,] + params.sig2[i-1]*params.mu[i-1]) / 
            (params.tau2[i-1]*params.lambda.i[i,] + params.sig2[i-1]),
        sqrt((params.tau2[i-1]*params.sig2[i-1]) /
            (params.tau2[i-1]*params.lambda.i[i,] + params.sig2[i-1])))

    # Update mu
    params.mu[i] = rnorm(1, 
        (n*s2*mean(params.mu.i[i,]) + params.tau2[i-1]*m) /
            (n*s2+params.tau2[i-1]),
        sqrt((s2 * params.tau2[i-1]) /
            (n*s2+params.tau2[i-1])))

    # Update tau2
    params.tau2[i] = 1/rgamma(1, c + n/2, d + 1/2*sum((params.mu.i[i,] - params.mu[i])^2))

    # Update sig2
    params.sig2[i] = 1/rgamma(1, a + n/2, b + 1/2*sum(params.lambda.i[i,] *
        (y - params.mu.i[i,])^2))

    if (i == (nburn + nmcmc))
        cat("\n")
    }

params.lambda.i = tail(params.lambda.i, nmcmc)
params.mu.i = tail(params.mu.i, nmcmc)
params.mu = tail(params.mu, nmcmc)
params.tau2 = tail(params.tau2, nmcmc)
params.sig2 = tail(params.sig2, nmcmc)


hpd.1.mu = hpd.uni(params.mu)
hpd.1.sig2 = hpd.uni(params.sig2)
hpd.1.tau2 = hpd.uni(params.tau2)

pdf("./figs/m_post.pdf", width = 12, height = 8)
par(mfrow = c(3,2), mar = c(4.1, 2.1, 2.1, 1.1))
hpd.plot(density(params.mu), hpd.1.mu, main = expression("M1: Posterior for"~ mu),
    xlab = expression(mu), cex.lab = 1.3, col1 = 'firebrick')
hpd.plot(density(m2.mu), hpd.2.mu, main = expression("M2: Posterior for"~ mu),
    xlab = expression(mu), cex.lab = 1.3)
hpd.plot(density(params.sig2), hpd.1.sig2, main = expression("M1: Posterior for"~ sigma^2),
    xlab = expression(sigma^2), cex.lab = 1.3, col1 = 'firebrick')
hpd.plot(density(m2.sig2), hpd.2.sig2, main = expression("M2: Posterior for"~ sigma^2),
    xlab = expression(sigma^2), cex.lab = 1.3)
hpd.plot(density(params.tau2), hpd.1.tau2, main = expression("M1: Posterior for"~ tau^2),
    xlab = expression(tau^2), cex.lab = 1.3, col1 = 'firebrick')
par(mfrow = c(1,1), mar = c(5.1, 4.1, 4.1, 2.1))
dev.off()

#plot(params.sig2, type='l')
#plot(params.tau2, type='l')
#plot(params.mu, type='l')
#matplot(tail(params.mu.i, 1000), type='l')
#matplot(tail(params.lambda.i, 1000), type='l')

### Posterior predictive distribution (for each observation)
pred.y = matrix(0, nmcmc, n)
for (j in 1:n)
    pred.y[,j] = rnorm(nmcmc, params.mu.i[,j], sqrt(params.sig2 / params.lambda.i[,j]))
qq = apply(pred.y, 2, quantile, c(0.025, 0.975))

pdf("./figs/m1_obs_fit.pdf", height = 6, width = 6)
plot(y, apply(pred.y, 2, mean), pch = 20, ylim = range(qq), xlab = "Observed", ylab = "Fitted",
    cex.lab = 1.3, main = "M1: Predictions based on random effects", col = 'darkred')
segments(x0 = y, x1 = y, y0 = qq[1,], y1 = qq[2,], col = 'firebrick')
abline(0, 1)
dev.off()

### Posterior predictive for a new observation
pred.mu.0 = rnorm(nmcmc, params.mu, sqrt(params.tau2))
pred.lambda.0 = rgamma(nmcmc, df/2, df/2)
pred.y.0 = rnorm(nmcmc, pred.mu.0, sqrt(params.sig2 / pred.lambda.0))
#plot(density(y))
#lines(density(pred.y.0), col = 'darkgreen', lwd = 3)
#for (j in 1:n){
#    dens = density(pred.y[,j])
#    lines(dens$x, 0.1*dens$y, col = rgb(0, 1, 0, 0.25))
#    }
#points(y, rep(0, n), pch = 20)

# Bayesian goodness-of-fit
m1.all = cbind(params.lambda.i, params.mu.i, params.mu, params.sig2, params.tau2)
m1.pvals = bayes.gof(y, m1.all, function(y, x) pnorm(y, x[(n+1):(2*n)], sqrt(x[2*n+2] / x[1:n])))
mean(m1.pvals > 0.05)
mean(m1.pvals > 0.5)





### Model 2
m2.mu = double(nburn + nmcmc)
m2.sig2 = double(nburn + nmcmc)
m2.sig2[1] = 1

for (i in 2:(nburn + nmcmc)){
    # Update mu
    m2.mu[i] = rnorm(1, (n*s2*mean(y) + m2.sig2[i-1]*m) / (n*s2 + m2.sig2[i-1]),
        sqrt((s2*m2.sig2[i-1])/(n*s2+m2.sig2[i-1])))

    # Update sig2
    m2.sig2[i] = 1/rgamma(1, a + n/2, b + 0.5*sum((y - m2.mu[i])^2))
    }

m2.mu = tail(m2.mu, nmcmc)
m2.sig2 = tail(m2.sig2, nmcmc)

hpd.2.mu = hpd.uni(m2.mu)
hpd.2.sig2 = hpd.uni(m2.sig2)

pdf("./figs/m2_post.pdf", width = 6, height = 8.5)
par(mfrow = c(2,1), mar = c(5.1, 4.1, 4.1, 2.1))
par(mfrow = c(1,1), mar = c(5.1, 4.1, 4.1, 2.1))
dev.off()

m2.pred.y = rnorm(nmcmc, m2.mu, sqrt(m2.sig2))

pdf("./figs/post_pred.pdf", width = 6, height = 6)
hist(y, xlab = "Log Interevent Time", cex.lab = 1.3,
    col = 'gray50', border = 'white', freq = FALSE, breaks = 20,
    main = "Posterior predictive distributions")
lines(density(pred.y.0), col = 'firebrick', lwd = 3)
lines(density(m2.pred.y), col = 'dodgerblue', lwd = 3)
legend("topleft", box.lty = 0, col = c("firebrick", "dodgerblue", 1),
    legend = c("M1", "M2", "Data"), lty=1, lwd = c(2,2,3), cex = 1.5)
dev.off()

# Bayesian goodness-of-fit
m1.all = cbind(params.lambda.i, params.mu.i, params.mu, params.sig2, params.tau2)
m2.pvals = bayes.gof(y, cbind(m2.mu, m2.sig2), function(y, x) pnorm(y, x[1], sqrt(x[2])))
mean(m2.pvals)
mean(m2.pvals > 0.05)
mean(m2.pvals > 0.5)


pdf("./figs/gof.pdf", width = 6, height = 6)
plot(density(m2.pvals), col = 'dodgerblue', lwd = 2, xlab = "p-value",
    main = "Bayesian goodness-of-fit p-values", cex.lab = 1.3)
lines(density(m1.pvals), col = 'firebrick', lwd = 2)
legend("topleft", box.lty = 0, col = c("firebrick", "dodgerblue"),
    legend = c("M1", "M2"), lty=1, lwd = 2, cex = 1.5)
dev.off()





##### Model Comparison
### Posterior predictive loss
#c(sum((y - apply(pred.y, 2, mean))^2), sum(apply(pred.y, 2, var)))
c(sum((y - mean(pred.y.0))^2), n*var(pred.y.0))     # Model 1
sum(c(sum((y - mean(m2.pred.y))^2), n*var(m2.pred.y)))   # Model 2

### DIC
m1.dtheta = matrix(0, nmcmc, length(y))
for (i in 1:length(y))
    m1.dtheta[,i] = dnorm(y[i], params.mu.i[,i], sqrt(params.sig2 / params.lambda.i[,i]), log = TRUE)
m1.dtheta = -2*apply(m1.dtheta, 1, sum)
m1.DIC = mean(m1.dtheta) + var(m1.dtheta)/2

m2.dtheta = matrix(0, nmcmc, length(y))
for (i in 1:length(y))
    m2.dtheta[,i] = dnorm(y[i], m2.mu, sqrt(m2.sig2), log = TRUE)
m2.dtheta = -2*apply(m2.dtheta, 1, sum)
m2.DIC = mean(m2.dtheta) + var(m2.dtheta)/2

m1.DIC
m2.DIC


### Bayes factor (from Monte Carlo estimates)
B = 1000000 # Increase B!!!!
prior.means = rnorm(B, rnorm(B, m, sqrt(s2)), sqrt(1/rgamma(B, c, d)))
prior.sds = sqrt(1/rgamma(B, a, b) / rgamma(B, df/2, df/2))
bf1 = matrix(0, B, length(y))
for (i in 1:length(y))
    bf1[,i] = dnorm(y[i], prior.means, prior.sds, log = TRUE)
bf1 = apply(bf1, 1, sum)

prior.means = rnorm(B, m, sqrt(s2))
prior.sds = sqrt(1/rgamma(B, a, b))
bf2 = matrix(0, B, length(y))
for (i in 1:length(y))
    bf2[,i] = dnorm(y[i], prior.means, prior.sds, log = TRUE)
bf2 = apply(bf2, 1, sum)

BF12 = mean(exp(bf1)) / mean(exp(bf2))
BF12

