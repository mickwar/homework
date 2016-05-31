library(splines)
source("~/files/R/mcmc/bayes_functions.R")

dat = read.table("~/files/data/eu_poll.txt", header = TRUE)
n = NROW(dat)

#key_events = c("5-Jan-2016", "2-Feb-2016", "20-Feb-2016", "21-Feb-2016",
#    "20-Mar-2016", "15-Apr-2016", "18-Apr-2016", "22-Apr-2016")
#key_time = as.numeric(as.Date(key_events, "%d-%b-%Y")) - 16678
#key_time = key_time[c(1, 2, 5, 6, 8)]

logit = function(p) log(p/(1-p))
ilogit = function(x) 1/(1+exp(-x))

table(dat[,5:6])

### Re-order (so earliest poll is first) 
dat = dat[n:1,]
dat[,5] = as.character(dat[,5])
dat[,6] = as.character(dat[,6])

### Combine the companies which did few polls
dat[dat[,5] == "BMG" | dat[,5] == "Opinium" | dat[,5] == "Panelbase", 5] = "Other"
#dat[dat[,5] == "BMG" | dat[,5] == "Opinium" | dat[,5] == "Panelbase" |
#    dat[,5] == "IpsosMori" | dat[,5] == "TNS", 5] = "Other"

### Orginaize data
# Observations
z = logit(dat[,2] / (dat[,2] + dat[,3]))

# Time vector
time = as.numeric(as.Date(dat[,1], "%d-%b-%Y")) - 16678

# Design matrices
# Treat Other and Online as baselines
company = matrix(model.matrix(lm(z ~ 0 + factor(dat[,5]))), nrow=n)
colnames(company) = c(sort(unique(dat[,5])))
company = company[,-which(colnames(company) == "Other")]
type = matrix(model.matrix(lm(z ~ 0 + factor(dat[,6]))), nrow=n)
colnames(type) = c(sort(unique(dat[,6])))
type = type[,-which(colnames(type) == "online")]

plot(time, z, type='l')

#####
nburn = 200
nmcmc = 10000

funpar = function(k){
    design.matrices = list(
        matrix(rep(1, n), n, 1, dimnames = list(NULL, "Intercept")),
        cbind("Intercept"=1, time),
        cbind("Intercept"=1, company),
        cbind("Intercept"=1, type),
        cbind("Intercept"=1, time, company),
        cbind("Intercept"=1, time, type),
        cbind("Intercept"=1, company, type),
        cbind("Intercept"=1, time, company, type))

    emp.cdf = double(length(design.matrices))

    for (model.n in 1:length(design.matrices)){
        y = z[-k]
        n = length(y)
        X = matrix(design.matrices[[model.n]][-k,], nrow = n)
        xstar = design.matrices[[model.n]][k,]

        ### Gibbs sampler
        p = NCOL(X)
        prior.m = rep(0, p)
        prior.g = n
        prior.a = 0
        prior.b = 0

        A = t(X) %*% X
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

xlim = c(Inf, -Inf)
ylim = c(Inf, -Inf)
for (j in 1:8){
    dens = density(probs[,j])
    xlim[1] = min(xlim[1], dens$x)
    xlim[2] = max(xlim[2], dens$x)
    ylim[1] = min(ylim[1], dens$y)
    ylim[2] = max(ylim[2], dens$y)
    }

pdf("./figs/loo.pdf", width = 6, height = 6)
plot(0,type='n',col='firebrick',lwd=2,xlab="Posterior predictive probability",
    main="Leave-one-out Analysis",cex.lab=1.3,xlim=xlim,ylim=ylim,ylab="Density")
for (j in 1:8)
    lines(density(probs[,j]), col = j, lwd = 2)
legend("topleft", box.lty = 0, col = 1:8,
    legend = paste0("M", 1:8), lty=1, lwd = 2, cex = 1.2)
dev.off()

### K-S test
t(t(apply(probs, 2, function(x) ks.test(jitter(x), 'punif')$p.value)))
