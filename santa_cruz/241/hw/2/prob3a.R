### Simple Bayesian model
library(MASS)
source("~/files/R/mcmc/bayes_functions.R")

### Posterior predictive loss criterion
pplc = function(y, ypred, k = Inf){
    n = length(y)
    vars = apply(ypred, 2, var)
    means = apply(ypred, 2, mean)

    factor = k / (k + 1)
    if (k == Inf)
        factor = 1

    return (sum(vars) + factor*sum((y - means)^2))
    }

dat = read.table("~/files/data/fabric.txt", header = TRUE)
x = dat$length / 100
y = dat$faults
n = length(y)
ord = order(y)

set.seed(1)
jity = jitter(y)

plot(x, y, pch = 20)

calc.post = function(param){
    # likelihood
    out = sum(dpois(y, param[1]*exp(param[2]*x), log = TRUE))
    # priors
    out = out + dgamma(param[1], prior.theta.a, prior.theta.b, log = TRUE)
    out = out + dnorm(param[2], prior.beta.a, prior.beta.b, log = TRUE)
    return (out)
    }

prior.theta.a = 1
prior.theta.b = 1
prior.beta.a = 0
prior.beta.b = 1


nburn = 5000
nmcmc = 20000

nparam = 2
params = matrix(0, nburn + nmcmc, nparam)
accept = double(nburn + nmcmc)
cand.sig = diag(0.1, nparam)

params[1,] = c(rgamma(1, prior.theta.a, prior.theta.b), rnorm(1, prior.beta.a, prior.beta.b))

lower = c(0, -Inf)
upper = c(Inf, Inf)
window = 100

post = calc.post(params[1,])


for (i in 2:(nburn + nmcmc)){
    cat("\r", i, "/", nburn+nmcmc)
    params[i,] = params[i-1,]
    cand = mvrnorm(1, params[i-1,], cand.sig)
    if (all(cand > lower) && all(cand < upper)){
        cand.post = calc.post(cand)
        if (log(runif(1)) <= cand.post - post){
            post = cand.post
            params[i,] = cand
            accept[i] = 1
            }
        }
    if ((floor(i/window) == i/window) && (i <= nburn))
        cand.sig = autotune(mean(accept[(i-window+1):i]), target = 0.234, k = window/50) *
            (cand.sig + window * var(params[(i-window+1):i,]) / i)
    }

params = tail(params, nmcmc)
accept = tail(accept, nmcmc)


### Trace plots
plot(params, type='l')
plot(params[,1], type='l')
plot(params[,2], type='l')

mean(accept)

### predictions
y0 = matrix(0, nmcmc, n)
for (i in 1:n)
    y0[,i] = rpois(nmcmc, params[,1]*exp(params[,2]*x[i]))
q0 = apply(y0, 2, quantile, c(0.025, 0.975))
m0 = apply(y0, 2, mean)

### Predictions at each x_i
plot(x, y, pch = 1, ylim = range(c(y, q0)), lwd = 1.5)
segments(x0 = x, y0 = q0[1,], y1 = q0[2,], col = 'forestgreen')
points(x, m0, col = 'darkgreen', pch = 20)


### Observed vs. Fitted plot
plot(0, type='n', xlim = range(y), ylim = range(c(y, y0)),
    xlab = "Observed", ylab = "Fitted", main = "Predictions", cex.lab = 1.3, cex.main = 2)
segments(x0 = jity, y0 = q0[1,], y1 = q0[2,], col = 'forestgreen')
points(jity, m0, col = 'darkgreen', pch = 20)
abline(0, 1, lwd = 3, lty = 2)

pplc(y, y0, Inf)
# 950.46
