### Problem 3
library(KMsurv)
library(survival)
library(mwBASE)
library(xtable)

data(tongue)
tongue$time = tongue$time / 10

aneuploid = tongue[tongue$type == 1,]
diploid = tongue[tongue$type == 2,]

calc.post = function(y, p){
    alpha = p[1]
    gamma = p[2]
    beta = p[3]
    if (alpha <= 0 || gamma <= 0)
        return (-Inf)

    x = 2-y[,1] # 1 is anueuploid, 0 is diploid
    t = y[,2]
    A = which(y[,3] == 1)   # the non-censored (i.e. observed)

    # Likelihood
    out = length(A)*log(alpha*gamma) + (alpha - 1)*sum(log(t[A])) + 
        sum(x[A]*beta) - gamma*sum((t^alpha) * exp(x*beta))

    # Priors
    out = out - log(alpha)
    out = out - log(gamma)
    out = out + dnorm(beta, 0, 10, log = TRUE)
    
    return (out)
    }

mod1 = mcmc_sampler(tongue, calc.post, 3, nburn = 20000, nmcmc = 30000)
mean(mod1$accept)

plot(mod1$params[,c(1,2)])
plot(mod1$params[,c(1,3)])
plot(mod1$params[,c(2,3)])
colMeans(mod1$params)

beta0 = -log(mod1$params[,2])/mod1$params[,1]
beta1 = -mod1$params[,3]/mod1$params[,1]
sigma = 1/mod1$params[,1]

mean(beta0)
mean(beta1)
mean(sigma)

var(beta0)
var(beta1)
var(sigma)

mean(mod1$params[,2]/exp(mod1$params[,3]))
