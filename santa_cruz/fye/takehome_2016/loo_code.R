library(splines)
source("~/files/R/mcmc/bayes_functions.R")
dat = read.table("~/files/data/fye_2016_data.txt", header = TRUE)

###
z = dat[,1]
rad = dat[,2]
temp = dat[,3]
wind = dat[,4]
m = length(z)

funpar = function(k){
    design.matrices = list(
        matrix(rep(1, m), m, 1, dimnames = list(NULL, "Intercept")),
        cbind("Intercept"=1, rad, temp, wind),
        cbind("Intercept"=1, temp, wind, temp*wind),
        cbind("Intercept"=1, temp, wind, temp^2, wind^2, temp*wind),
        cbind("Intercept"=1, rad, temp, wind, temp*wind),
        cbind("Intercept"=1, rad, ns(temp, knots = c(75)), ns(wind, knots = c(10)),
            temp*wind),
        cbind("Intercept"=1, rad, ns(temp, knots = c(70, 80)), ns(wind, knots = c(5, 12)),
            temp*wind),
        cbind("Intercept"=1, rad, ns(temp, knots = c(70, 80)), ns(wind, knots = c(5, 12)),
            ns(temp, knots = c(70, 80))*ns(wind, knots = c(5, 12))))#,
    qnum = length(design.matrices)

    emp.cdf = double(qnum)

    for (model.n in 1:qnum){
        y = z[-k]
        n = length(y)
        X = matrix(design.matrices[[model.n]][-k,], nrow = n)
        xstar = design.matrices[[model.n]][k,]

        nburn = 200
        nmcmc = 10000

        p = NCOL(X)
        prior.m = rep(0, p)
        prior.g = n
        prior.a = 0
        prior.b = 0

        A = t(X) %*% X
        max(eigen(A)$values) / min(eigen(A)$values)
        chol.A = t(chol(solve(A)))
        post.mean = 1/(prior.g + 1) * (prior.g * solve(t(X) %*% X) %*% t(X) %*% y + prior.m)

        param.beta = matrix(0, nburn + nmcmc, p)
        param.sig2 = double(nburn + nmcmc)
        param.sig2[1] = 1

        for (i in 2:(nburn + nmcmc)){
            param.beta[i,] = post.mean +
                sqrt(prior.g * param.sig2[i-1] / (prior.g + 1)) * chol.A %*% rnorm(p, 0, 1)
            param.sig2[i] = 1/rgamma(1,
                prior.a + n/2 + p/2, 
                prior.b + 0.5*sum((y - X %*% param.beta[i,])^2) +
                    0.5/prior.g * t(param.beta[i,] - prior.m) %*% A %*% (param.beta[i,] - prior.m))
            }

        # Truncate
        param.beta = tail(param.beta, nmcmc)
        param.sig2 = tail(param.sig2, nmcmc)

        ### Posterior predictions
        pred.y.0 = rnorm(nmcmc, param.beta %*% xstar, sqrt(param.sig2))

        emp.cdf[model.n] = mean(pred.y.0 <= z[k])
        }
    return (emp.cdf)
    }

library(foreach)
library(doMC)
registerDoMC(4)

probs = foreach(k = 1:length(z), .combine=rbind) %dopar% funpar(k)

qnum = 8
xlim = c(Inf, -Inf)
ylim = c(Inf, -Inf)
for (j in 1:qnum){
    dens = density(probs[,j])
    xlim[1] = min(xlim[1], dens$x)
    xlim[2] = max(xlim[2], dens$x)
    ylim[1] = min(ylim[1], dens$y)
    ylim[2] = max(ylim[2], dens$y)
    }

pdf("figs/loo.pdf", width = 6, height = 6)
plot(0,type='n',col='firebrick',lwd=2,xlab="Posterior predictive probability",
    main="Leave-one-out Analysis",cex.lab=1.3,xlim=xlim,ylim=ylim,ylab="Density")
for (j in 1:qnum)
    lines(density(probs[,j]), col = j, lwd = 2)
legend("topleft", box.lty = 0, col = 1:8,
    legend = paste0("M", 1:8), lty=1, lwd = 2, cex = 1.2)
dev.off()


### K-S test
t(t(apply(probs, 2, function(x) ks.test(jitter(x), 'punif')$p.value)))
