### Problem 3
library(KMsurv)
library(survival)
library(mwBASE)

data(tongue)
tongue$time = tongue$time / 10

aneuploid = tongue[tongue$type == 1,]
diploid = tongue[tongue$type == 2,]


calc.post = function(x, p){
    t = x[,2]   # time to death or on-study time, weeks
    nu = x[,3]  # 0=alive (censored), 1=dead
    obs = which(nu == 1)
    alpha = p[1]
    lambda = p[2]
    if (alpha <= 0 || lambda <= 0)
        return (-Inf)

    # Likelihood
    out = length(obs)*log(alpha) + length(obs)*log(lambda) +
        (alpha - 1)*sum(log(t[obs])) - lambda*sum(t^alpha)

    # Priors
#   out = out - log(alpha)
#   out = out - log(lambda)
    return (out)
    }

out1 = mcmc_sampler(aneuploid, calc.post, nparam = 2, nburn = 40000, window = 500, nmcmc = 20000)
out2 = mcmc_sampler(diploid, calc.post, nparam = 2, nburn = 40000, window = 500, nmcmc = 20000)

mean(out1$accept)
plot(out1$params[,1], type='l')
plot(out1$params[,2], type='l')
plot(out1$params[,c(1,2)], pch = 16)

mean(out2$accept)
plot(out2$params[,1], type='l')
plot(out2$params[,2], type='l')
plot(out2$params[,c(1,2)], pch = 16)

plot(out1$params[,c(1,2)], pch = 16, xlim = c(0.3, 1.5), ylim = c(0, 0.7), col = rgb(1,0,0,0.1))
points(out2$params[,c(1,2)], pch = 16, col = rgb(0,1,0,0.1))

# survreg(Surv(time, delta) ~ 1, data = aneuploid, dist = "weibull")
colMeans(out1$params)

# survreg(Surv(time, delta) ~ 1, data = diploid, dist = "weibull")
colMeans(out2$params)

plot_hpd(density(out1$params[,1]), hpd_mult(out1$params[,1]), "dodgerblue", fade = 0.7)
plot_hpd(density(out2$params[,1]), hpd_mult(out2$params[,1]), "yellow", add = TRUE, fade = 0.7)
