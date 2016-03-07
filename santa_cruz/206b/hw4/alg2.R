library(coda)
library(MASS)

autotune = function(accept, target = 0.25, k = 2.5)
    (1+(cosh(accept-target)-1)*(k-1)/(cosh(target-
        ceiling(accept-target))-1))^sign(accept-target)

x = c(5.25644,6.630203,1.548099,1.978492,2.154677,3.18829,0.8640084,2.298181,2.023496,
    2.956788,3.788934,0.8952826,1.179618,2.539608,5.115685,2.870664,2.752188,
    0.5418666,0.7017399,1.163273,3.130179,3.910224,1.39601,1.728189,0.988848)
n = length(x)
plot(density(x))

nburn = 20000
nmcmc = 20000

params = matrix(0, nburn + nmcmc, 2)
accept = double(nburn + nmcmc)
window = 500
#cand.sig = matrix(c(0.0818, 0.1649, 0.1649, 0.4201), 2, 2)
#cand.sig = matrix(c(0.0748, 0.1502, 0.1502, 0.3884), 2, 2)
cand.sig = matrix(c(0.0753, 0.0597, 0.0597, 0.0594), 2, 2)

keep.sig = NULL
keep.acc = NULL
keep.sig[length(keep.sig)+1] = determinant(cand.sig)$modulus

for (i in 2:(nburn + nmcmc)){
    params[i,] = params[i-1,]

    # Joint update
    cand = mvrnorm(1, params[i,], cand.sig)

    tau = params[i,1]
    eta = params[i,2]

    cand.post = cand[1]*(n*exp(cand[2])+2) + 3*cand[2] - n*lgamma(exp(cand[2])) +
        (exp(cand[2])-1)*sum(log(x)) - exp(cand[1])*(2+sum(x)) - exp(cand[2])
    curr.post = tau*(n*exp(eta)+2) + 3*eta - n*lgamma(exp(eta)) +
        (exp(eta)-1)*sum(log(x)) - exp(tau)*(2+sum(x)) - exp(eta)

    if (log(runif(1)) < cand.post - curr.post){
        params[i,] = cand
        accept[i] = 1
        }

    if ((floor(i/window) == i/window) && (i <= nburn)){
#       cand.sig = cov(params[(i-window+1):i,])
        cand.sig = (cand.sig + autotune(mean(accept[(i-window+1):i]),
            k = max(window / 50, 1.5))^2 * cov(params[(i-window+1):i,]))/2
        keep.sig[length(keep.sig)+1] = determinant(cand.sig)$modulus
        keep.acc[length(keep.acc)+1] = mean(accept[(i-window+1):i])
        }
    }

plot(keep.sig, type='l')
keep.acc

cand.sig
var(params)

params = tail(params, nmcmc)
accept = tail(accept, nmcmc)

mean(accept)

apply(params, 2, mean)
var(params)
apply(params, 2, effectiveSize)

plot(params, pch = 20, type='n')
segments(x0 = params[-nmcmc,1], y0 = params[-nmcmc,2],
    x1 = params[-1,1], y1 = params[-1,2],
    col = rgb(seq(0, 1, length = nmcmc-1), 0, 0))
y = mvrnorm(1000, apply(params, 2, mean), cand.sig)
points(y, col = 'green', pch = 20, cex = 0.5)

pred = rgamma(nmcmc, params[,2], params[,1])
plot(density(x))
lines(density(pred), col = 'seagreen', lwd = 2)
