### Proportional hazards model
library(KMsurv)
library(survival)
library(mwBASE)
library(xtable)

data(tongue)
tongue$time = tongue$time / 10

### Weibull baseline hazard
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



### Piece-wise constant hazards
ss = c("0%"=0, quantile(tongue[,2], seq(0.05, 1.00, by = 0.05)))
#ss = c("0%"=0, quantile(tongue[,2], seq(0.10, 0.90, by = 0.10)), max(tongue[,2]+0.01))
#ss = seq(0, 40, length = 50)
J = length(ss)-1

# which interval each observation (failed or censored) is in
ints = sapply(tongue[,2], function(x) which.max(x <= ss[-1]))


calc.post = function(y, p){
    beta0 = p[1]
    beta1 = p[2]
    lambda = p[-c(1,2)]

    if (any(lambda <= 0))
        return (-Inf)

    x = 2 - y[,1]
    t = y[,2]
    A = which(y[,3] == 1)   # set of observations where not censored
    H = c(0, cumsum(lambda*(ss[-1] - ss[-(J+1)])))    # cumulative hazard

    # Likelihood
    out = sum(log(lambda[ints[A]]) + beta0 + x[A]*beta1) -
        sum((lambda[ints]*(t - ss[ints]) + H[ints])*exp(beta0 + x * beta1))

    # Priors
    out = out + dnorm(beta0, 0, 10, log = TRUE)
    out = out + dnorm(beta1, 0, 10, log = TRUE)
    out = out + sum(dgamma(lambda, 1, 1/10, log = TRUE))

    return (out)
    }

mod2 = mcmc_sampler(tongue, calc.post, nparam = (2 + J), nburn = 50000, nmcmc = 20000)

mean(mod2$accept)
plot(mod2$params[,1], type='l')
plot(mod2$params[,2], type='l')
plot(mod2$params[,3], type='l')
plot(mod2$params[,4], type='l')

colMeans(mod2$params)


# for (i in 1:NCOL(mod2$params)){
#     plot(mod2$params[,i], type='l')
#     readline()
#     }

dip.survival = t(apply(mod2$params,1,function(p) exp(-exp(p[1])*cumsum(p[-c(1,2)]*(ss[-1]-ss[-(J+1)])))))
ane.survival = t(apply(mod2$params,1,function(p) exp(-exp(p[1]+p[2])*cumsum(p[-c(1,2)]*(ss[-1]-ss[-(J+1)])))))

dip.mm = apply(dip.survival, 2, mean)
dip.vv = apply(dip.survival, 2, quantile, c(0.025, 0.975))
ane.mm = apply(ane.survival, 2, mean)
ane.vv = apply(ane.survival, 2, quantile, c(0.025, 0.975))

# Kaplan-Meier curves
plot(survfit(Surv(time[type == 1], delta[type == 1]) ~ 1, data = tongue), col = 'dodgerblue', lwd=1)
lines(survfit(Surv(time[type == 2], delta[type == 2]) ~ 1, data = tongue), col = 'firebrick1', lwd=1)

# Posterior survival curves
plot(0, type='n', xlim = range(ss), ylim = c(0, 1), bty = 'n')
lines(ss, c(1, dip.mm), col = 'firebrick1', lwd = 3)
lines(ss, c(1, dip.vv[1,]), col = 'firebrick1',)
lines(ss, c(1, dip.vv[2,]), col = 'firebrick1',)
lines(ss, c(1, ane.mm), col = 'dodgerblue', lwd = 3)
lines(ss, c(1, ane.vv[1,]), col = 'dodgerblue',)
lines(ss, c(1, ane.vv[2,]), col = 'dodgerblue',)

# Together
plot(survfit(Surv(time[type == 1], delta[type == 1]) ~ 1, data = tongue), col = 'dodgerblue')
lines(survfit(Surv(time[type == 2], delta[type == 2]) ~ 1, data = tongue), col = 'firebrick1')
lines(ss, c(1, dip.mm), col = 'firebrick1', lwd = 3)
lines(ss, c(1, dip.vv[1,]), col = 'firebrick1',)
lines(ss, c(1, dip.vv[2,]), col = 'firebrick1',)
lines(ss, c(1, ane.mm), col = 'dodgerblue', lwd = 3)
lines(ss, c(1, ane.vv[1,]), col = 'dodgerblue',)
lines(ss, c(1, ane.vv[2,]), col = 'dodgerblue',)

