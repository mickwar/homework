source("./sampling.R")
library(mwBASE)
# library(KMsurv)
# library(survival)

dat = as.matrix(read.table("./kidney-dat.txt", header = TRUE))
dat[,4] = 2 - dat[,4]

calc_post = function(dat, param){
    cluster = dat[,1]
    time = dat[,2]
    nu = dat[,3]
    x = dat[,4:5]

    eta = param[1]
    beta = param[2:3]
    gamma = param[4]
    alpha = param[5]
    w = tail(param, 38)

    if (any(param[-(2:3)] <= 0))
        return (-Inf)

    # Likelihood
    out = sum(log(alpha) + log(gamma) + log(w[cluster[nu == 1]]) +
        (alpha-1)*log(time[nu == 1]) + x[nu == 1,] %*% beta) -
        sum(w[cluster]*gamma*(time^alpha)*exp(x %*% beta))

    # Prior
    out = out + dgamma(eta, 0.001, 0.001, log = TRUE)
    out = out + sum(dnorm(beta, 0, sqrt(10^3), log = TRUE))
    out = out + dgamma(gamma, 0.001, 0.001, log = TRUE)
    out = out + dgamma(alpha, 0.001, 0.001, log = TRUE)
    out = out + sum(dgamma(w, eta, eta, log = TRUE))

    return (out)
    }

out = get_samples(nburn =  200000, nmcmc = 50000)

chain_init = apply(tail(out$params, 500), 2, mean)
sig_init = out$cand_sig

out = get_samples(nburn =  10000, nmcmc = 10000, chain_init = chain_init,
    cand_sig = sig_init)


mean(out$accept)

par(mfrow = c(9,5), mar = c(0,0,0,0))
for (i in 1:43)
    plot(out$params[,i], type='l', axes=FALSE)
par(mfrow = c(1,1), mar = c(5.1,4.1,4.1,2.1))

pairs(out$params[,1:5], pch = 16, col = rgb(seq(0, 1, length = 5000), 0, 0))
pairs(out$params[,6:10], pch = 16, col = rgb(seq(0, 1, length = 5000), 0, 0))
pairs(out$params[,11:15], pch = 16, col = rgb(seq(0, 1, length = 5000), 0, 0))

