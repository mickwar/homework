### Hierarchical Bayesian model
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
#x = as.numeric(scale(dat$length))
y = dat$faults
n = length(y)
ord = order(y)

set.seed(1)
jity = jitter(y)

par(mfrow = c(1,1), mar = c(5.1, 4.1, 4.1, 2.1))
plot(x, y, pch = 20)

calc.post = function(param){
    theta = param[1:n]
     beta = param[n+1]
       mu = param[n+2]
     zeta = param[n+3]

    # likelihood
    out = sum(dpois(y, theta*exp(beta*x), log = TRUE))

    # priors
    out = out + sum(dgamma(theta, zeta, zeta/mu, log = TRUE))
    out = out + dnorm(beta, prior.beta.a, prior.beta.b, log = TRUE)
    out = out + dgamma(mu, prior.mu.a, prior.mu.b, log = TRUE)
    out = out + dgamma(zeta, prior.zeta.a, prior.zeta.b, log = TRUE)

    return (out)
    }


prior.beta.a = 0
prior.beta.b = 1
prior.mu.a = 1
prior.mu.b = 1/2
prior.zeta.a = 1
prior.zeta.b = 1/2


nburn = 100000
nmcmc = 100000

#last = tail(params, 1)
nparam = n + 3 # n thetas, 1 for beta, mu, and zeta
params = matrix(0, nburn + nmcmc, nparam)
accept = double(nburn + nmcmc)
cand.sig = diag(0.01, nparam)
#params[1,] = last

params[1, c(n+2, n+3)] = rgamma(2, c(prior.mu.a, prior.zeta.a), c(prior.mu.b, prior.zeta.b))
params[1,1:(n+1)] = c(rgamma(n, params[1,n+3], params[1,n+3]/params[1,n+2]),
    rnorm(1, prior.beta.a, prior.beta.b))

lower = c(rep(0, n), -Inf, 0, 0)
upper = rep(Inf, nparam)
window = 400

post = calc.post(params[1,])


for (i in 2:(nburn + nmcmc)){
    if (floor(i/window) == i/window)
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



apply(params[,n+c(1,3,2)], 2, mean)
t(apply(params[,n+c(1,3,2)], 2, quantile, c(0, 0.025, 0.5, 0.975, 1)))


### Trace plots
par(mfrow = c(8,4), mar = c(0,0,0,0))
for (i in 1:n)
    plot(params[,i], type='l')

hpds = apply(params[,n+(1:3)], 2, hpd.uni)
par(mfrow = c(3,2), mar = c(5.1,4.1,4.1,2.1))
for (i in 1:3){
    plot(params[,n+i], type='l')
    hpd.plot(density(params[,n+i]), hpds[,i])
    }

mean(accept)

### predictions
theta0 = rgamma(nmcmc, params[,n+3], params[,n+3] / params[,n+2])
y0 = matrix(0, nmcmc, n)
for (i in 1:n)
    y0[,i] = rpois(nmcmc, theta0*exp(params[,n+1]*x[i]))
q0 = apply(y0, 2, quantile, c(0.025, 0.975))
m0 = apply(y0, 2, mean)


#y0 = matrix(0, nmcmc, n)
#for (i in 1:n)
#    y0[,i] = rpois(nmcmc, params[,i]*exp(params[,n+1]*x[i]))
#q0 = apply(y0, 2, quantile, c(0.025, 0.975))
#m0 = apply(y0, 2, mean)

### Predictions at each x_i
pdf("./figs/pred_3b.pdf", width = 12, height = 6)
par(mfrow = c(1,2), mar = c(4.1,4.1,2.1, 1.1))
plot(x, y, pch = 1, ylim = range(c(y, q0)), lwd = 1.5, main = "Predictions")
segments(x0 = x, y0 = q0[1,], y1 = q0[2,], col = 'forestgreen')
points(x, m0, col = 'darkgreen', pch = 20)
legend("topleft", box.lty = 0, legend = "Data", pch = 1, cex = 1.5)

### Observed vs. Fitted plot
plot(0, type='n', xlim = range(y), ylim = range(c(y, q0)),
    xlab = "Observed", ylab = "Fitted", main = "Observed vs. Fitted", cex.lab = 1.3)
segments(x0 = jity, y0 = q0[1,], y1 = q0[2,], col = 'forestgreen')
points(jity, m0, col = 'darkgreen', pch = 20)
abline(0, 1, lwd = 3, lty = 2)
dev.off()

pplc(y, y0, 0)                      # 888.66
pplc(y, y0, Inf)                    # 1538.68
pplc(y, y0, Inf) - pplc(y, y0, 0)   # 650.01

