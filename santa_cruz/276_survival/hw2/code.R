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

#beta0 = -log(mod1$params[,2])/mod1$params[,1]
#beta1 = -mod1$params[,3]/mod1$params[,1]
#sigma = 1/mod1$params[,1]

boxplot(mod1$params)

tt = seq(0, max(tongue[,2]), length = 100)
dip.survival = t(apply(mod1$params, 1, function(p) exp(-p[2]*(tt^p[1]))))
ane.survival = t(apply(mod1$params, 1, function(p) exp(-p[2]*(tt^p[1])*exp(p[3]))))

dip.mm = apply(dip.survival, 2, mean)
dip.vv = apply(dip.survival, 2, quantile, c(0.025, 0.975))
ane.mm = apply(ane.survival, 2, mean)
ane.vv = apply(ane.survival, 2, quantile, c(0.025, 0.975))

plot(0, type='n', xlim = c(0, max(tongue[,2])), ylim = c(0, 1), bty = 'n',
    main = "Weibull baseline model", xlab = "Time (t)", ylab = "Survival S(t)")
lines(survfit(Surv(time[type == 1], delta[type == 1]) ~ 1, data = tongue),
    col = col_fade('dodgerblue', 0.5), lty = 2, lwd = c(2, 1, 1))
lines(survfit(Surv(time[type == 2], delta[type == 2]) ~ 1, data = tongue),
    col = col_fade('firebrick1', 0.5), lty = 2, lwd = c(2, 1, 1))
lines(tt, dip.mm, col = 'firebrick1', lwd = 3)
lines(tt, dip.vv[1,], col = 'firebrick1',)
lines(tt, dip.vv[2,], col = 'firebrick1',)
lines(tt, ane.mm, col = 'dodgerblue', lwd = 3)
lines(tt, ane.vv[1,], col = 'dodgerblue',)
lines(tt, ane.vv[2,], col = 'dodgerblue',)

par(mfrow = c(2,2))
plot_hpd(mod1$params[,1], col1 = 'orange', bty = 'n', xlab = expression(alpha),
    main = expression("Posterior for"~ alpha), cex.main = 2, cex.lab = 1.5, sub = "M1")
plot_hpd(mod1$params[,2], col1 = 'orange', bty = 'n', xlab = expression(gamma),
    main = expression("Posterior for"~ gamma), cex.main = 2, cex.lab = 1.5, sub = "M1")
plot_hpd(mod1$params[,3], col1 = 'orange', bty = 'n', xlab = expression(beta),
    main = expression("Posterior for"~ beta), cex.main = 2, cex.lab = 1.5, sub = "M1")
par(mfrow = c(1,1), mar = c(5.1, 4.1, 4.1, 2.1))

# dip.hazard = t(apply(mod1$params, 1, function(p) p[1]*p[2]*ss^(p[1]-1)))
# ane.hazard = t(apply(mod1$params, 1, function(p) p[1]*p[2]*ss^(p[1]-1)*exp(p[3])))
# 
# plot(ss, colMeans(dip.hazard))
# plot(ss, colMeans(ane.hazard))



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
# calc2 = function(y, p){
#     beta0 = p[1]
#     beta1 = p[2]
#     lambda = p[-c(1,2)]
# 
#     if (any(lambda <= 0))
#         return (-Inf)
# 
#     x = 2 - y[,1]
#     t = y[,2]
#     nu = y[,3]
# 
#     out = 0
#     for (i in 1:length(t)){
#         tmp = 1
#         for (j in 1:J){
#             tmp = tmp * (lambda[j]^(ints[i]==j)*exp(beta0 + beta1*(x[i] == 1)))^(nu[i] == 1)
#             if (j == 1){
#                 tmp = tmp * exp(-(ints[i] == j)*(lambda[j]*(t[i] - ss[j]))*
#                     exp(beta0 + beta1*(x[i] == 1)))
#             } else {
#                 tmp = tmp * exp(-(ints[i] == j)*(lambda[j]*(t[i] - ss[j]) +
#                     sum(lambda[1:(j-1)]*(ss[2:j] - ss[1:(j-1)])))*
#                     exp(beta0 + beta1*(x[i] == 1)))
#                 }
#             }
#         out = out + log(tmp)
#         }
# 
#     out = out + dnorm(beta0, 0, 10, log = TRUE)
#     out = out + dnorm(beta1, 0, 10, log = TRUE)
#     out = out + sum(dgamma(lambda, 1, 1/10, log = TRUE))
#     return (out)
#     }
# 
# rr = runif(J+1)
# calc.post(tongue, c(0, rr))
# calc2(tongue, c(0, rr))


mod2 = mcmc_sampler(tongue, calc.post, nparam = (2 + J), nburn = 50000, nmcmc = 20000)

mean(mod2$accept)
colMeans(mod2$params)

#colMeans(mod2$params[,-c(1,2)]*exp(mod2$params[,1]))


dip.survival = t(apply(mod2$params,1,function(p) exp(-exp(p[1])*cumsum(p[-c(1,2)]*(ss[-1]-ss[-(J+1)])))))
ane.survival = t(apply(mod2$params,1,function(p) exp(-exp(p[1]+p[2])*cumsum(p[-c(1,2)]*(ss[-1]-ss[-(J+1)])))))

dip.mm = apply(dip.survival, 2, mean)
dip.vv = apply(dip.survival, 2, quantile, c(0.025, 0.975))
ane.mm = apply(ane.survival, 2, mean)
ane.vv = apply(ane.survival, 2, quantile, c(0.025, 0.975))

# Kaplan-Meier curves with posterior estimates of survival function
plot(0, type='n', xlim = c(0, max(tongue[,2])), ylim = c(0, 1), bty = 'n',
    main = "Piece-wise consant hazard model", xlab = "Time (t)", ylab = "Survival S(t)")
lines(survfit(Surv(time[type == 1], delta[type == 1]) ~ 1, data = tongue),
    col = col_fade('dodgerblue', 0.5), lty = 2, lwd = c(2, 1, 1))
lines(survfit(Surv(time[type == 2], delta[type == 2]) ~ 1, data = tongue),
    col = col_fade('firebrick1', 0.5), lty = 2, lwd = c(2, 1, 1))
lines(ss, c(1, dip.mm), col = 'firebrick1', lwd = 3)
lines(ss, c(1, dip.vv[1,]), col = 'firebrick1')
lines(ss, c(1, dip.vv[2,]), col = 'firebrick1')
lines(ss, c(1, ane.mm), col = 'dodgerblue', lwd = 3)
lines(ss, c(1, ane.vv[1,]), col = 'dodgerblue')
lines(ss, c(1, ane.vv[2,]), col = 'dodgerblue')



plot_hpd(mod2$params[,1], col1 = 'forestgreen', bty = 'n', xlab = expression(beta),
    main = expression("Posterior for"~ beta), cex.main = 2, cex.lab = 1.5, sub = "M2")

boxplot(mod2$params[,-c(1,2)]*exp(mod2$params[,1]), axes = FALSE, horizontal = TRUE,
    main = "Boxplots for posterior hazards", xlab = "Hazard", at = J:1,
    ylab = "Index", col = 'forestgreen')
axis(1);
axis(2, at = c(1, 6, 11, 16, 20), labels = c("20", "15", "10", "5", "1"))



### Gamma process
# Same intervals as used in the piece-wise model
ss = c("0%"=0, quantile(tongue[,2], seq(0.05, 1.00, by = 0.05)))
J = length(ss) - 1

# Unique death times interval
ss = c(0, unique(sort(tongue[tongue[,3] == 1,2])), max(tongue[,2]) + 0.01)
J = length(ss) - 1

risk = rep(list(NULL), J)
fail = rep(list(NULL), J)
risk_not_fail = rep(list(NULL), J)

for (j in 1:J){
    risk[[j]] = which(tongue[,2] > ss[j])
    fail[[j]] = which((tongue[,3] == 1) & (tongue[,2] > ss[j]) & (tongue[,2] <= ss[j+1]))
    risk_not_fail[[j]] = risk[[j]][which(!(risk[[j]] %in% fail[[j]]))]
    }

calc.post = function(y, p){
    beta0 = p[1]
    beta1 = p[2]
    hj = p[-c(1,2)]

    if (any(hj <= 0))
        return (-Inf)

    x = 2 - y[,1]
    t = y[,2]

    out = 0
    for (j in 1:J){
        z1 = x[risk_not_fail[[j]]]
        z2 = x[fail[[j]]]
        if (length(z1) == 0)
            z1 = 0
        if (length(z2) == 0)
            z2 = 0
        out = out -hj[j] * sum(exp(beta0 + z1 * beta1)) + 
            sum(log(1 - exp(-hj[j]*exp(beta0 + z2 * beta1)))) +
            dgamma(hj[j], 1, 1/10, log = TRUE)
        }

    out = out + dnorm(beta0, 0, 10, log = TRUE)
    out = out + dnorm(beta1, 0, 10, log = TRUE)
    return (out)
    }

mod3 = mcmc_sampler(tongue, calc.post, nparam = (J + 2), nburn = 50000, nmcmc = 20000)

mean(mod3$accept)
plot(mod3$params[,c(1,2)], pch = 16)
pairs(mod3$params[,1:6], pch = 16)

colMeans(mod3$params)
colMeans(mod3$params[,-(1:2)]*exp(mod3$params[,1]))
