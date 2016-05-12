library(MASS)
source("~/files/R/mcmc/bayes_functions.R")
x = c(14, 30, 49, 70, 33, 15)
ints = c(-Inf, 66, 68, 70, 72, 74, Inf)
n = sum(x)

calc.post = function(param){
    mu = param[1]
    sig = sqrt(param[2])
    sum(x*log(pnorm(tail(ints, length(x)), mu, sig) -
        pnorm(head(ints, length(x)), mu, sig)))
    }

nburn = 10000
nmcmc = 100000

nparam = 2
params = matrix(0, nburn + nmcmc, nparam)
accept = double(nburn + nmcmc)
cand.sig = diag(0.1, nparam)

params[1,] = c(70, 3)

lower = c(-Inf, 0)
upper = c(Inf, Inf)
window = 200

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
mean(accept)


params.1 = params

apply(params, 2, mean)
apply(params, 2, var)


par(mfrow=c(2,2), mar = c(4.1,2.1,2.1,1.1))
plot(density(params[,1]))
plot(params[,1], type='l')
plot(density(params[,2]))
plot(params[,2], type='l')
par(mfrow=c(1,1), mar = c(5.1,4.1,4.1,2.1))


### Latent variables
library(truncnorm)

calc.post = function(param, Z){
    mu = param[1]
    sig = sqrt(param[2])
    sum(dnorm(Z, mu, sig, log = TRUE)) + 0
#   sum(x*log(pnorm(tail(ints, length(x)), mu, sig) -
#       pnorm(head(ints, length(x)), mu, sig)))
    }

nburn = 10000
nmcmc = 100000

nparam = 2
params = matrix(0, nburn + nmcmc, nparam)
accept = double(nburn + nmcmc)
cand.sig = diag(0.1, nparam)
latent.Z = matrix(0, nburn + nmcmc, n)


lower = c(-Inf, 0)
upper = c(Inf, Inf)
window = 200

index = list(
    (          0+1):sum(x[1:1]),
    (sum(x[1:1])+1):sum(x[1:2]),
    (sum(x[1:2])+1):sum(x[1:3]),
    (sum(x[1:3])+1):sum(x[1:4]),
    (sum(x[1:4])+1):sum(x[1:5]),
    (sum(x[1:5])+1):sum(x[1:6]))
    
params[1,] = c(70, 3)
for (i in 1:length(x))
    latent.Z[1, index[[i]]] = rtruncnorm(x[i], a = ints[i], b = ints[i+1],
        mean = params[1,1], sd = sqrt(params[1,2]))

post = calc.post(params[1,], latent.Z[1,])
set.seed(1)

for (i in 2:(nburn + nmcmc)){
    if (floor(i/window) == i/window)
        cat("\r", i, "/", nburn+nmcmc)

    # Update Z
    for (j in 1:length(x))
        latent.Z[i, index[[j]]] = rtruncnorm(x[j], a = ints[j], b = ints[j+1],
            mean = params[i-1,1], sd = sqrt(params[i-1,2]))

    # Update mu, sigma^2
    params[i,] = params[i-1,]
    cand = mvrnorm(1, params[i-1,], cand.sig)
    if (all(cand > lower) && all(cand < upper)){
        post = calc.post(params[i,], latent.Z[i,])
        cand.post = calc.post(cand, latent.Z[i,])
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
mean(accept)

params.2 = params

apply(params, 2, mean)
apply(params, 2, var)

par(mfrow=c(2,2), mar = c(4.1,2.1,2.1,1.1))
plot(density(params[,1]))
plot(params[,1], type='l')
plot(density(params[,2]))
plot(params[,2], type='l')
par(mfrow=c(1,1), mar = c(5.1,4.1,4.1,2.1))



apply(latent.Z, 2, mean)
apply(latent.Z, 2, var)


apply(params.1, 2, mean)
apply(params.2, 2, mean)
apply(params.1, 2, var)
apply(params.2, 2, var)
par(mfrow=c(2,1), mar = c(4.1,2.1,2.1,1.1))
plot(density(params.1[,1]))
lines(density(params.2[,1]), col = 'red')
plot(density(params.1[,2]))
lines(density(params.2[,2]), col = 'red')
par(mfrow=c(1,1), mar = c(5.1,4.1,4.1,2.1))

