library(coda)
set.seed(1)

x = c(5.25644,6.630203,1.548099,1.978492,2.154677,3.18829,0.8640084,2.298181,2.023496,
    2.956788,3.788934,0.8952826,1.179618,2.539608,5.115685,2.870664,2.752188,
    0.5418666,0.7017399,1.163273,3.130179,3.910224,1.39601,1.728189,0.988848)
n = length(x)
plot(density(x))

nburn = 10000
nmcmc = 100000

params = matrix(0, nburn + nmcmc, 2)
accept = double(nburn + nmcmc)
cand.sig = 0.6

# 1: theta, 2: log nu
params[1,1] = rgamma(1, 2 + n * exp(0), 2 + sum(x))
params[1,2] = 0

for (i in 2:(nburn + nmcmc)){
    params[i,] = params[i-1,]

    # Update theta
    params[i,1] = rgamma(1, 2 + n*exp(params[i,2]), 2 + sum(x))

    # Update log nu
    cand = rnorm(1, params[i,2], cand.sig) # propose log nu

    theta = params[i,1]
    lognu = params[i,2]
    cand.post = n*exp(cand)*log(theta) + 2*cand - n*lgamma(exp(cand)) + 
        (exp(cand)-1) * sum(log(x)) - exp(cand)
    curr.post = n*exp(lognu)*log(theta) + 2*lognu - n*lgamma(exp(lognu)) +
        (exp(lognu)-1) * sum(log(x)) - exp(lognu)

    if (log(runif(1)) < cand.post - curr.post){
        params[i,2] = cand
        accept[i] = 1
        }
    }

params = tail(params, nmcmc)
params[,2] = exp(params[,2])
accept = tail(accept, nmcmc)

mean(accept)

apply(params, 2, mean)
apply(params, 2, var)
apply(params, 2, effectiveSize)

plot(params, pch = 20, type='n')
segments(x0 = params[-nmcmc,1], y0 = params[-nmcmc,2],
    x1 = params[-1,1], y1 = params[-1,2],
    col = rgb(seq(0, 1, length = nmcmc-1), 0, 0))

pred = rgamma(nmcmc, params[,2], params[,1])
plot(density(x))
lines(density(pred), col = 'seagreen', lwd = 2)
