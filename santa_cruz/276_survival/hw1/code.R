### Problem 3
library(KMsurv)
library(survival)
library(mwBASE)
library(xtable)

data(tongue)
tongue$time = tongue$time / 10

aneuploid = tongue[tongue$type == 1,]
diploid = tongue[tongue$type == 2,]


### part a

# survreg(Surv(time, delta) ~ 1, data = aneuploid, dist = "weibull")
# survreg(Surv(time, delta) ~ 1, data = diploid, dist = "weibull")

calc.post = function(x, p, priors = 1){
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
    if (priors == 1){
        out = out - log(alpha)
        out = out - log(lambda)
        }
    if (priors == 2){
        out = out + dgamma(alpha, 1, 1, log = TRUE)
        out = out + dgamma(lambda, 1, 1, log = TRUE)
        }
    if (priors == 3){
        out = out + dgamma(alpha, 1, 0.1, log = TRUE)
        out = out + dgamma(lambda, 1, 0.1, log = TRUE)
        }
    if (priors == 4){
        out = out + dgamma(alpha, 50, 50, log = TRUE)
        out = out + dgamma(lambda, 50, 50, log = TRUE)
        }
    return (out)
    }

ane = NULL
dip = NULL
for (i in 1:4){
    ane[[i]] = mcmc_sampler(aneuploid, function(x, p) calc.post(x, p, priors = i),
        nparam = 2, nburn = 5000, window = 500, nmcmc = 20000)
    dip[[i]] = mcmc_sampler(diploid, function(x, p) calc.post(x, p, priors = i),
        nparam = 2, nburn = 5000, window = 500, nmcmc = 20000)
    }
sapply(ane, function(x) mean(x$accept))
sapply(dip, function(x) mean(x$accept))

par(mfrow = c(2, 2))
lapply(ane, function(x) plot(x$params, pch = 16, bty = 'n', xlab = "alpha", ylab = "lambda"))
lapply(dip, function(x) plot(x$params, pch = 16, bty = 'n', xlab = "alpha", ylab = "lambda"))
par(mfrow = c(1, 1))

for (i in 1:4){
    ane[[i]]$hpds = apply(ane[[i]]$params, 2, hpd_mult)
    dip[[i]]$hpds = apply(dip[[i]]$params, 2, hpd_mult)
    }

pdf("a_sens_alpha.pdf", height = 9, width = 9)
par(mfrow = c(2, 2), mar = c(4.1, 2.1, 2.1, 1.1), oma = c(0, 0, 2, 0))
for (i in 1:4){
    tmp.x = density(ane[[i]]$params[,1])
    tmp.y = density(dip[[i]]$params[,1])
    plot_hpd(tmp.x, ane[[i]]$hpds[,1],
        xlim = range(tmp.x$x, tmp.y$x), ylim = range(tmp.x$y, tmp.y$y),
        "dodgerblue", fade = 0.7, bty = 'n', main = paste0("Prior set ", i),
        xlab = expression(alpha), cex.lab = 1.5)
    plot_hpd(tmp.y, dip[[i]]$hpds[,1],
        "firebrick1", add = TRUE, fade = 0.7)
    }
make_phantom(c("Aneuplooid", "Diploid"), 1, c("dodgerblue", "firebrick1"),
    outer = TRUE, line = 0.5, cex.main = 1.5, sep = "  ")
make_phantom(c("Aneuplooid", "Diploid"), 2, c("dodgerblue", "firebrick1"),
    outer = TRUE, line = 0.5, cex.main = 1.5, sep = "  ")
par(mfrow = c(1, 1), mar = c(5.1, 4.1, 4.1, 2.1), oma = c(0, 0, 0, 0))
dev.off()

pdf("a_sens_lambda.pdf", height = 9, width = 9)
par(mfrow = c(2, 2), mar = c(4.1, 2.1, 2.1, 1.1), oma = c(0, 0, 2, 0))
for (i in 1:4){
    tmp.x = density(ane[[i]]$params[,2])
    tmp.y = density(dip[[i]]$params[,2])
    plot_hpd(tmp.x, ane[[i]]$hpds[,2],
        xlim = range(tmp.x$x, tmp.y$x), ylim = range(tmp.x$y, tmp.y$y),
        "dodgerblue", fade = 0.7, bty = 'n', main = paste0("Prior set ", i),
        xlab = expression(lambda), cex.lab = 1.5)
    plot_hpd(tmp.y, dip[[i]]$hpds[,2],
        "firebrick1", add = TRUE, fade = 0.7)
    }
make_phantom(c("Aneuplooid", "Diploid"), 1, c("dodgerblue", "firebrick1"),
    outer = TRUE, line = 0.5, cex.main = 1.5, sep = "  ")
make_phantom(c("Aneuplooid", "Diploid"), 2, c("dodgerblue", "firebrick1"),
    outer = TRUE, line = 0.5, cex.main = 1.5, sep = "  ")
par(mfrow = c(1, 1), mar = c(5.1, 4.1, 4.1, 2.1), oma = c(0, 0, 0, 0))
dev.off()

t(rbind(apply(ane[[1]]$params, 2, function(x) c(mean(x), var(x))), ane[[1]]$hpds))
t(rbind(apply(dip[[1]]$params, 2, function(x) c(mean(x), var(x))), dip[[1]]$hpds))


### part b

calc.post = function(x, p, likelihood = "ev"){
    t = x[,2]   # time to death or on-study time, weeks
    nu = x[,3]  # 0=alive (censored), 1=dead
    obs = which(nu == 1)
    beta0 = p[1]
    beta1 = p[2]
    sig = p[3]
    if (sig <= 0)
        return (-Inf)

    like.base = switch(likelihood,
        "ev"        = function(x) sum(x - exp(x)),
        "normal"    = function(x) sum(dnorm(x, 0, 1, log = TRUE)),
        "logistic"  = function(x) sum(-x - 2*log(1 + exp(-x))))
    surv.base = switch(likelihood,
        "ev"        = function(x) sum(-exp(x)),
        "normal"    = function(x) sum(pnorm(x, 0, 1, log = TRUE, lower.tail = FALSE)),
        "logistic"  = function(x) sum(-x - log(1 + exp(-x))))


    # Likelihood
    z = (log(t) - beta0 - (x[,1] == 1)*beta1) / sig

    out = like.base(z[obs]) - length(obs)*log(sig) + surv.base(z[-obs])

    # Priors
    out = out + dnorm(beta0, 0, 10, log = TRUE)
    out = out + dnorm(beta1, 0, 10, log = TRUE)
    out = out - log(sig)
    return (out)
    }

aft.weibull     = mcmc_sampler(tongue, function(x, p) calc.post(x, p, "ev"),
    nparam = 3, nburn = 5000, window = 500, nmcmc = 20000, chain_init = c(0, 0, 1))
aft.lognormal   = mcmc_sampler(tongue, function(x, p) calc.post(x, p, "normal"),
    nparam = 3, nburn = 5000, window = 500, nmcmc = 20000, chain_init = c(0, 0, 1))
aft.loglogistic = mcmc_sampler(tongue, function(x, p) calc.post(x, p, "logistic"),
    nparam = 3, nburn = 5000, window = 500, nmcmc = 20000, chain_init = c(0, 0, 1))

# Sanity check
survreg(Surv(time, delta) ~ as.factor(ifelse(type == 1, 2, 1)), data = tongue, dist = "weibull")
apply(aft.weibull$params, 2, mean)

survreg(Surv(time, delta) ~ as.factor(ifelse(type == 1, 2, 1)), data = tongue, dist = "lognormal")
apply(aft.lognormal$params, 2, mean)

survreg(Surv(time, delta) ~ as.factor(ifelse(type == 1, 2, 1)), data = tongue, dist = "loglogistic")
apply(aft.loglogistic$params, 2, mean)

# some MCMC diagnostics
plot(aft.weibull$params[,1], type='l')
plot(aft.weibull$params[,2], type='l')
plot(aft.weibull$params[,3], type='l')
mean(aft.weibull$accept)

plot(aft.lognormal$params[,1], type='l')
plot(aft.lognormal$params[,2], type='l')
plot(aft.lognormal$params[,3], type='l')
mean(aft.lognormal$accept)

plot(aft.loglogistic$params[,1], type='l')
plot(aft.loglogistic$params[,2], type='l')
plot(aft.loglogistic$params[,3], type='l')
mean(aft.loglogistic$accept)

aft.weibull$hpds = apply(aft.weibull$params, 2, hpd_mult)
aft.lognormal$hpds = apply(aft.lognormal$params, 2, hpd_mult)
aft.loglogistic$hpds = apply(aft.loglogistic$params, 2, hpd_mult)

pdf("bc_posterior.pdf", width = 9, height = 12)
lab = c("beta0", "beta1", "sig")
par(mfrow = c(3, 1), mar = c(4.1, 2.1, 2.1, 1.1), oma = c(0, 0, 2, 0))
for (i in 1:3){
    tmp.x = density(aft.weibull$params[,i])
    tmp.y = density(aft.lognormal$params[,i])
    tmp.z = density(aft.loglogistic$params[,i])
    plot_hpd(tmp.x, aft.weibull$hpds[,i],
        xlim = range(tmp.x$x, tmp.y$x, tmp.z$x), ylim = range(tmp.x$y, tmp.y$y, tmp.z$y),
        "dodgerblue", fade = 0.7, bty = 'n', xlab = lab[i], main = "", cex.lab = 1.5)
    plot_hpd(tmp.y, aft.lognormal$hpds[,i],
        "firebrick1", fade = 0.7, add = TRUE)
    plot_hpd(tmp.z, aft.loglogistic$hpds[,i],
        "forestgreen", fade = 0.7, add = TRUE)
    }
make_phantom(c("Weibull", "Log-normal", "Log-logistic"), c(1,2,3),
    c("dodgerblue", "firebrick1", "forestgreen"), outer = TRUE, line = 0.5,
    cex.main = 1.5, sep = "   ")
par(mfrow = c(1, 1), mar = c(5.1, 4.1, 4.1, 2.1), oma = c(0, 0, 0, 0))
dev.off()



weib.dtheta = apply(aft.weibull$params, 1, function(p) {
    -2*(calc.post(tongue, p, "ev") - dnorm(p[1], 0, 10, log = TRUE) -
        dnorm(p[2], 0, 10, log = TRUE) + log(p[3]))
    })
logn.dtheta = apply(aft.lognormal$params, 1, function(p) {
    -2*(calc.post(tongue, p, "normal") - dnorm(p[1], 0, 10, log = TRUE) -
        dnorm(p[2], 0, 10, log = TRUE) + log(p[3]))
    })
logl.dtheta = apply(aft.loglogistic$params, 1, function(p) {
    -2*(calc.post(tongue, p, "logistic") - dnorm(p[1], 0, 10, log = TRUE) -
        dnorm(p[2], 0, 10, log = TRUE) + log(p[3]))
    })

# DICs
mean(weib.dtheta) + 0.5*var(weib.dtheta)
mean(logn.dtheta) + 0.5*var(logn.dtheta)
mean(logl.dtheta) + 0.5*var(logl.dtheta)

t(rbind(apply(aft.weibull$params, 2, function(x) c(mean(x), var(x))), aft.weibull$hpds))
t(rbind(apply(aft.lognormal$params, 2, function(x) c(mean(x), var(x))), aft.lognormal$hpds))
t(rbind(apply(aft.loglogistic$params, 2, function(x) c(mean(x), var(x))), aft.loglogistic$hpds))

mean(exp(aft.weibull$params[,2]))
