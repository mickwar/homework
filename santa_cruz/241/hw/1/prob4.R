library(MCMCpack)

cdf = function(x, t){
    out = length(t)
    for (i in 1:length(t))
        out[i] = mean(x <= t[i])
    return (out)
    }
gen.y1 = function(n)
    rpois(n, 5)
gen.y2 = function(n){
    out = double(n)
    mix = sample(1:2, n, replace = TRUE, prob = c(0.7, 0.3))
    out = ifelse(mix == 1, rpois(sum(mix == 1), 3), out)
    out = ifelse(mix == 2, rpois(sum(mix == 2), 11), out)
    return (out)
    }

n = 300

set.seed(1)
y = gen.y2(n)

nj = as.vector(table(y))
ystar = as.numeric(names(table(y))) # sort(unique(y))
nstar = length(ystar)

calc.like = function(alpha, lambda){
    logzpar = function(z, m)
        lgamma(z + m) - lgamma(z)

    nstar * log(alpha) - logzpar(alpha, n) +
        sum(dpois(ystar, lambda, log = TRUE)) + 
        sum(logzpar(alpha * dpois(ystar, lambda, log = FALSE) + 1, nj - 1))
    }

F0.dist.post = function(t, alpha, lambda)
    alpha/(alpha+n)*ppois(t, lambda) + 1/(alpha+n)*sapply(t, function(t) sum(y <= t))

tvec = seq(0, 30, by = 1)

nburn = 1000
nmcmc = 20000

alpha = double(nburn + nmcmc)
lambda = double(nburn + nmcmc)
fparam = matrix(0, nburn + nmcmc, length(tvec))
al.accept = double(nburn + nmcmc)
#sig = matrix(c(29.1, -0.059, -0.059, 0.222), 2, 2) * 2.4^2/2
sig = matrix(c(326, -1.644, -1.644, 0.116), 2, 2) * 2.4^2

a.alpha = 1/10
b.alpha = 1/10
a.lambda = 5/2
b.lambda = 1/2

### Initial values
alpha[1] = rgamma(1, a.alpha, b.alpha)
lambda[1] = rgamma(1, a.lambda, b.lambda)
fparam[1,] = F0.dist.post(c(tvec[-1], Inf), alpha[1], lambda[1]) - 
    F0.dist.post(tvec, alpha[1], lambda[1])
fparam[1,] = F0.dist.post(c(tvec[-1], Inf), alpha[1], lambda[1]) - 
    F0.dist.post(tvec, alpha[1], lambda[1])
#fparam[1,] = cumsum(rdirichlet(1, (alpha[1] + n)*fparam[1,]))

for (i in 2:(nburn + nmcmc)){
    cat("\r", i, "/", nburn+nmcmc)
    # update alpha, lambda
    alpha[i] = alpha[i-1]
    lambda[i] = lambda[i-1]
    cand = mvrnorm(1, c(alpha[i-1], lambda[i-1]), sig)
    if (all(cand > 0)){
        post = calc.like(alpha[i-1], lambda[i-1]) +
            dgamma(alpha[i-1], a.alpha, b.alpha, log = TRUE) +
            dgamma(lambda[i-1], a.lambda, b.lambda, log = TRUE)
        cand.post = calc.like(cand[1], cand[2]) +
            dgamma(cand[1], a.alpha, b.alpha, log = TRUE) +
            dgamma(cand[2], a.lambda, b.lambda, log = TRUE)
        if (log(runif(1)) <= cand.post - post){
            al.accept[i] = 1
            alpha[i] = cand[1]
            lambda[i] = cand[2]
            }
        }

    # update F (fparam)
    fparam[i,] = F0.dist.post(c(tvec[-1]-0.5, Inf), alpha[i], lambda[i]) -
        F0.dist.post(tvec-0.5, alpha[i], lambda[i])
    fparam[i,] = F0.dist.post(c(tvec[-1]-0.5, Inf), alpha[i], lambda[i]) -
        F0.dist.post(tvec-0.5, alpha[i], lambda[i])
    fparam[i,] = rdirichlet(1, (alpha[i] + n)*fparam[i,])
    if (i == (nburn + nmcmc))
        cat("\n")
    }

alpha = tail(alpha, nmcmc)
lambda = tail(lambda, nmcmc)
fparam = tail(fparam, nmcmc)
al.accept = tail(al.accept, nmcmc)

mean(al.accept)
var(cbind(alpha, lambda))

par(mfrow = c(2,2))
plot(alpha, type = 'l')
plot(density(alpha), xlab = expression(alpha), main="")
plot(lambda, type = 'l')
plot(density(lambda), xlab = expression(lambda), main = "")

qlines = apply(fparam, 2, quantile, c(0, 0.025, 0.5, 0.975, 1))

par(mfrow = c(1,1))
plot(tvec, qlines[4,], pch = "_", col = 'darkgreen', cex = 2)
lines(tvec, qlines[3,], type='h', col = 'green', lwd = 3)
points(tvec, qlines[2,], pch = "_", col = 'darkgreen', cex = 2)
#lines(as.numeric(names(table(y)))+0.3, table(y)/n, col = 'blue', type='h', lwd = 2)
points(as.numeric(names(table(y))), table(y)/n, col = 'black', pch = 20, lwd = 2)

