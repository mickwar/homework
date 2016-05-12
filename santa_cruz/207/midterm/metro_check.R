source("~/files/R/mcmc/bayes_functions.R")
library(MASS)
dat = read.table("~/volcano.txt", header=TRUE)
y = log(dat[1:62,3])
n = length(y)

#n = 50
#set.seed(1)
#y = rt(n, 5)

plot(density(y))

digamma = function(x, alpha, beta)
    alpha*log(beta) - lgamma(alpha) - (alpha+1)*log(x) - beta/x
dt.new = function(x, df, mu, sig2)
    lgamma((df+1)/2) - lgamma(df/2) - 0.5*log(pi*df*sig2) -(df + 1)/2 *
        log(1 + 1/(df*sig2)*(x - mu)^2)

#sum(dt((y - mu.i)/sqrt(sig2), 5, log = TRUE))
#sum(dt.new(y, 5, mu.i[1,], sig2[1]))

calc.post = function(param){
    mu.i = head(param, n) 
    mu = param[n+1]
    sig2 = param[n+2]
    tau2 = param[n+3]
    # likelihood
#   out = sum(dt((y - mu.i)/sqrt(sig2), 5, log = TRUE))
    out = sum(dt.new(y, 5, mu.i, sig2))
    # priors
    out = out + sum(dnorm(mu.i, mu, sqrt(tau2), log = TRUE))
    out = out + dnorm(mu, m, sqrt(s2), log = TRUE)
    out = out + digamma(sig2, a, b)
    out = out + digamma(tau2, c, d)
    return (out)
    }

# Hyperpriors
m = 0   # mu mean
s2 = 10 # mu variance
a = 3   # sig2 alpha
b = 3   # sig2 beta
c = 3   # tau2 alpha
d = 3   # tau2 beta

nburn = 20000
nmcmc = 10000
window = 500

nparam = n + 3
params = matrix(0, nburn + nmcmc, nparam)
accept = matrix(0, nburn+nmcmc, nparam)
sigs = rep(1, nparam)

lower = c(rep(-Inf, n), -Inf, 0, 0)
upper = c(rep(Inf, n), Inf, Inf, Inf)

#params[1,1:n] = tail(mu.i, 1)
#params[1,n+1] = tail(mu, 1)
#params[1,n+2] = tail(sig2, 1)
#params[1,n+3] = tail(tau2, 1)

params[1,] = 1
post = calc.post(params[1,])

cand.param = params[1,]

for (i in 2:(nburn+nmcmc)){
    if (floor(i/window) == i/window)
        cat("\r", i, "/", nburn+nmcmc)
    params[i,] = params[i-1,]
    for (j in 1:nparam){
        cand = rnorm(1, params[i,j], sigs[j])
        if (cand >= lower[j] && cand <= upper[j]){
            cand.param[j] = cand
            cand.post = calc.post(cand.param)
            if (log(runif(1)) < cand.post - post){
                post = cand.post
                params[i,j] = cand
                accept[i,j] = 1
            } else {
                cand.param[j] = params[i,j]
                }
        } else {
            cand.param[j] = params[i,j]
            }
        }
    if (floor(i/window) == i/window && i <= nburn)
        sigs = sigs*autotune(apply(accept[(i-window+1):i,], 2,
            mean), k = max(window/50, 1.1))
    if (i == (nburn+nmcmc))
        cat("\n")
    }


#for (i in 2:(nburn + nmcmc)){
#    if (floor(i/window) == i/window)
#        cat("\r", i, "/", nburn+nmcmc)
#    params[i,] = params[i-1,]
#    cand = mvrnorm(1, params[i-1,], cand.sig)
#    if (all(cand > lower) && all(cand < upper)){
#        cand.post = calc.post(cand)
#        if (log(runif(1)) <= cand.post - post){
#            post = cand.post
#            params[i,] = cand
#            accept[i] = 1
#            }
#        }
#    if ((floor(i/window) == i/window) && (i <= nburn))
#        cand.sig = autotune(mean(accept[(i-window+1):i]), target = 0.234, k = window/50) *
#            (cand.sig + window * var(params[(i-window+1):i,]) / i)
#    }

params = tail(params, nmcmc)
accept = tail(accept, nmcmc)

mu.i = params[,1:n]
mu = params[,n+1]
sig2 = params[,n+2]
tau2 = params[,n+3]

apply(accept, 2, mean)

plot(mu, type='l')
plot(sig2, type='l')
plot(tau2, type='l')


###
pred.y = matrix(0, nmcmc, n)
for (j in 1:n)
    pred.y[,j] = mu.i[,j] + sqrt(sig2)*rt(nmcmc, 5)
qq = apply(pred.y, 2, quantile, c(0.025, 0.975))

plot(y, apply(pred.y, 2, mean), pch = 20, ylim = range(qq))
segments(x0 = y, x1 = y, y0 = qq[1,], y1 = qq[2,])
abline(0, 1)


# pred.y = matrix(0, nmcmc, n)
# for (j in 1:n)
#     pred.y[,j] = mu.i[,j] + sqrt(sig2)*rt(nmcmc, 5)
# plot(y, apply(pred.y, 2, mean), pch = 20)
# abline(0, 1)

pred.mu.0 = rnorm(nmcmc, mu, sqrt(tau2))
pred.y.0 = pred.mu.0 + sqrt(sig2)*rt(nmcmc, 5)

#plot(density(pred.mu.0))
plot(density(y))
#lines(density(c(pred.y)), col = 'blue', lwd = 2)
lines(density(pred.y.0), col = 'darkgreen', lwd = 2)
# curve(exp(dt.new(x, 5, 0, 0.5)), add = TRUE, col = 'red')
# 
# xx = seq(-5, 5, length = 100)
# plot(xx, dt((xx-0)/sqrt(1.5), 5), type='l')
# lines(xx, exp(dt.new(xx, 5, 0, sqrt(1.5))), col = 'red', type='l')
# 
# 
# for (j in 1:n){
#     plot(density(pred.y[,j])); abline(v = y[j])
#     readline()
#     }
