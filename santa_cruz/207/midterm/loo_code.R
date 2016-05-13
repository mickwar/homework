### Leave-one-out crossvalidation
source("~/files/R/mcmc/bayes_functions.R")
dat = read.table("~/files/data/volcano.txt", header=TRUE)
z = log(dat[1:62,3])
n = length(z)

set.seed(1)
funpar = function(k){
    y = z[-k]
    n = length(y)

    ### Model 1
    # Hyperpriors
    m = 6   # mu mean
    s2 = 10 # mu variance
    a = 3   # sig2 alpha
    b = 4   # sig2 beta
    c = 3   # tau2 alpha
    d = 4   # tau2 beta
    df = 5  # df for the t, fixed

    nburn = 1000
    nmcmc = 20000

    params.lambda.i = matrix(0, nburn + nmcmc, n)
    params.mu.i = matrix(0, nburn + nmcmc, n)
    params.mu = double(nburn + nmcmc)
    params.tau2 = double(nburn + nmcmc)
    params.sig2 = double(nburn + nmcmc)

    params.tau2[1] = 1
    params.sig2[1] = 1

    for (i in 2:(nburn + nmcmc)){
        # Update lambda's
        params.lambda.i[i,] = rgamma(n, (df + 1)/2, df/2 + 1/(2*params.sig2[i-1])*
            (y - params.mu.i[i-1,])^2)

        # Update mu.i's
        params.mu.i[i,] = rnorm(n,
            (params.tau2[i-1]*y*params.lambda.i[i,] + params.sig2[i-1]*params.mu[i-1]) / 
                (params.tau2[i-1]*params.lambda.i[i,] + params.sig2[i-1]),
            sqrt((params.tau2[i-1]*params.sig2[i-1]) /
                (params.tau2[i-1]*params.lambda.i[i,] + params.sig2[i-1])))

        # Update mu
        params.mu[i] = rnorm(1, 
            (n*s2*mean(params.mu.i[i,]) + params.tau2[i-1]*m) /
                (n*s2+params.tau2[i-1]),
            sqrt((s2 * params.tau2[i-1]) /
                (n*s2+params.tau2[i-1])))

        # Update tau2
        params.tau2[i] = 1/rgamma(1, c + n/2, d + 1/2*sum((params.mu.i[i,] - params.mu[i])^2))

        # Update sig2
        params.sig2[i] = 1/rgamma(1, a + n/2, b + 1/2*sum(params.lambda.i[i,] *
            (y - params.mu.i[i,])^2))

        }

    params.lambda.i = tail(params.lambda.i, nmcmc)
    params.mu.i = tail(params.mu.i, nmcmc)
    params.mu = tail(params.mu, nmcmc)
    params.tau2 = tail(params.tau2, nmcmc)
    params.sig2 = tail(params.sig2, nmcmc)

    ### Posterior predictive for a new observation
    pred.mu.0 = rnorm(nmcmc, params.mu, sqrt(params.tau2))
    pred.lambda.0 = rgamma(nmcmc, df/2, df/2)
    pred.y.0 = rnorm(nmcmc, pred.mu.0, sqrt(params.sig2 / pred.lambda.0))


    ### Model 2
    m2.mu = double(nburn + nmcmc)
    m2.sig2 = double(nburn + nmcmc)
    m2.sig2[1] = 1

    for (i in 2:(nburn + nmcmc)){
        # Update mu
        m2.mu[i] = rnorm(1, (n*s2*mean(y) + m2.sig2[i-1]*m) / (n*s2 + m2.sig2[i-1]),
            sqrt((s2*m2.sig2[i-1])/(n*s2+m2.sig2[i-1])))

        # Update sig2
        m2.sig2[i] = 1/rgamma(1, a + n/2, b + 0.5*sum((y - m2.mu[i])^2))
        }

    m2.mu = tail(m2.mu, nmcmc)
    m2.sig2 = tail(m2.sig2, nmcmc)

    m2.pred.y = rnorm(nmcmc, m2.mu, sqrt(m2.sig2))

    c(mean(pred.y.0 <= z[k]), mean(m2.pred.y <= z[k]))
    }

library(foreach)
library(doMC)
registerDoMC(4)

probs = foreach(k = 1:length(z), .combine=rbind) %dopar% funpar(k)
m1.prob = probs[,1]
m2.prob = probs[,2]

pdf("./figs/loo.pdf", width = 6, height = 6)
plot(density(m1.prob), col = 'firebrick', lwd = 2, xlab = "Posterior predictive probability",
    main = "Leave-one-out Analysis", cex.lab = 1.3)
lines(density(m2.prob), col = 'dodgerblue', lwd = 2)
legend("topleft", box.lty = 0, col = c("firebrick", "dodgerblue"),
    legend = c("M1", "M2"), lty=1, lwd = 2, cex = 1.5)
dev.off()

### K-S test
ks.test(m1.prob, 'punif')
ks.test(m2.prob, 'punif')

# Do a little jitter if there are ties
#ks.test(jitter(m1.prob), 'punif')
#ks.test(jitter(m2.prob), 'punif')

