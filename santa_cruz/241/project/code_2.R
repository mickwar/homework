source("~/files/R/mcmc/bayes_functions.R")
library(MASS)
library(tmvtnorm)

### Load the data (use load("./chen_gray_fake.RData") when real data not available)
source("./read_data.R")
dat = read_data()
#dat = c(head(dat, 7), tail(dat, 1))
n = length(dat)

ni = sapply(dat, function(x) nrow(x$xy))

temp = sapply(dat, function(x) x$xy[,2])
dat.mean = mean(unlist(temp))
dat.sd = sd(unlist(temp))

# Plot the data
cols = rainbow(n)
darkcols = character(n)
for (i in 1:n)
    darkcols[i] = col.mult(cols[i], "gray50")
plot(0, type='n', xlim = c(0, 0.27), ylim = c(-3, 3), xlab = "Plastic Strain",
    ylab = "Stress")
for (i in 1:n){
    dat[[i]]$xy[,2] = (dat[[i]]$xy[,2] - dat.mean) / dat.sd
    dat[[i]]$fixed = c(dat[[i]]$strain_rate, dat[[i]]$temperature, dat[[i]]$xy[,1])
    lines(dat[[i]]$xy[,1], dat[[i]]$xy[,2], type='l', col = cols[i], lwd = 2)
    # Combine some data
    }


### Simplified Johnson-Cook model
jc = function(x, theta){
    A = theta[1]
    B = theta[2]
    n = theta[3]
    C = theta[4]
    m = theta[5]
    strain_rate = x[1]
    temperature = x[2] / 3250 # ~melting (kelvin)
    eps = x[-(1:2)]
    return ((A + B * eps^n) * (1 + C * log(strain_rate)) * (1 - temperature^m))
    }

### Multivariate normal density (up to a constant)
dmvnorm = function(x, mu, sigma, log = TRUE){
    p = length(x)
    if (NROW(sigma) == 1){
#       out = 0.5*p*log(sigma) - 0.5/sigma * t(x-mu) %*% (x-mu)
        out = -p/2*log(2*pi) - 0.5*p*log(sigma) - 0.5/sigma * sum((x-mu)^2)
    } else {
        inv = solve(sigma)
        out = -p/2*log(2*pi) + 0.5*determinant(inv)$modulus[1] - 0.5*t(x-mu) %*% inv %*% (x-mu)
        }
    if (log) return (out)
    return (exp(out))
    }

### Draws from an inverse gamma
rinvgamma = function(n, shape, rate)
    1/rgamma(n, shape = shape, rate = rate)

### Draws from a Wishart distribution (several implementations)
# Bartlett decomposition (for random draws)
# get one random draw
rwishart = function(V, df){
    p = nrow(V)
    # function to randomize lower triangle of a matrix with independent N(0,1)s
    lower.rand = function(x){
        n = nrow(x)
        z = sequence(n)
        index = cbind(
            row = unlist(lapply(2:n, function(x) x:n), use.names = FALSE),
            col = rep(z[-length(z)], times = rev(tail(z, -1))-1))
        x[index] = rnorm(nrow(index))
        return (x)
        }
    L = t(chol(V))
    A = diag(sqrt(rchisq(p, df-(1:p)+1)))
    A = lower.rand(A)
    X = L %*% A
    return (X %*% t(X))
#   return (L %*% A %*% t(A) %*% t(L))
    }

### Priors
# Measurement variance
prior.tau.a = 10
prior.tau.b = 1

# DP precision
prior.alpha.a = 2
prior.alpha.b = 1

# G0 (baseline) mean
prior.mu.a = c(0.04, 0.07, 0.7, 0.05, 0.04)
prior.mu.B = 1*diag(5)

#theta.lower = c(0.00, 0.00, 0.05, 0.00, 0.01)
#theta.upper = c(0.15, 0.30, 0.95, 0.08, 0.09)
theta.lower = rep(-Inf, 5)
theta.upper = rep(Inf, 5)

# G0 (basline) covariance
prior.sigma.df = 5
prior.sigma.V = diag(1, 5)




### MCMC
nburn = 20000
nmcmc = 10000

nparam = 5 # Parameters in JC model
param.theta = matrix(0, nburn + nmcmc, nparam * n)                  # JC parameters
param.alpha = double(nburn + nmcmc)                                 # DP precision
param.mu = matrix(0, nburn + nmcmc, nparam)                         # G0 mean
param.sigma = rep(list(matrix(0, nparam, nparam)), nburn + nmcmc)   # G0 variance
param.tau = double(nburn + nmcmc)                                   # Measurement variance

# cand.sig = rep(list(0.1*diag(nparam)), n)
# accept = matrix(0, nburn + nmcmc, n)
# window = 200
cand.sig = 0.005 * diag(nparam)


param.alpha[1] = rgamma(1, prior.alpha.a, prior.alpha.b)
param.mu[1,] = prior.mu.a
param.sigma[[1]] = diag(0.001, nparam)
param.tau[1] = rinvgamma(1, prior.tau.a, prior.tau.b)
param.theta[1,] = rep(param.mu[1,], n)
clus.w = rep(1, n)
#param.theta[1,] = runif(nparam * n)
#clus.w = 1:n



for (iter in 2:(nburn + nmcmc)){
    cat("\r", iter, "/", nburn + nmcmc)

    # Current values
    c.theta = matrix(param.theta[iter-1,], n, nparam, byrow = TRUE)
    c.alpha = param.alpha[iter-1]
    c.mu = param.mu[iter-1,]
    c.sigma = param.sigma[[iter-1]]
    c.tau = param.tau[iter-1]

    # Update the thetas (algorithm 6)
    for (i in 1:n){
#       n.j.minus = table(clus.w)
#       temp.label = as.numeric(names(n.j.minus))
#       temp.ind = which(temp.label == clus.w[i])
#       n.j.minus[temp.ind] = n.j.minus[temp.ind] - 1

#       n.star.minus = length(temp.label) # but not really
#       temp.star = double(n.star.minus)
#       for (j in 1:n.star.minus)
#           temp.star[j] = which(clus.w == temp.label[j])[1]
#       theta.star.minus = matrix(c.theta[temp.star,], ncol = nparam)

        temp.prob = rep(1, n) / (n - 1 + c.alpha)
        temp.prob[i] = c.alpha / (n - 1 + c.alpha)

        # Draw a candidate value from prior conditional
        temp.samp = sample(n, 1, prob = temp.prob)
        if (temp.samp == i){ # Make a draw from G0
#           cand = mvrnorm(1, c.mu, c.sigma)
            cand = rtmvnorm(1, c.mu, c.sigma, lower = theta.lower, upper = theta.upper,
                algorithm = "gibbs", burn.in.sample = 50, start.value = c.mu)
        } else { # Make a draw from the existing groups
            cand = c.theta[temp.samp,]
            }

        # Decide whether to accept the candidate
        temp.cand = dmvnorm(dat[[i]]$xy[,2], jc(dat[[i]]$fixed, cand), c.tau)
        temp.curr = dmvnorm(dat[[i]]$xy[,2], jc(dat[[i]]$fixed, c.theta[i,]), c.tau)
        temp.cand - temp.curr
        if (log(runif(1)) < temp.cand - temp.curr){
            c.theta[i,] = cand
            # Acceptance rate? Neal doesn't seem to mention it
#           accept[iter,i] = 1

            # Change clus.w
            if (temp.samp == i){ # a new cluster was formed
                clus.w[i] = which.min(1:n %in% clus.w[-i])
            } else {            # part of existing cluster
                clus.w[i] = clus.w[temp.samp]
                }

            # Re-arrange clus.w for empty clusters
            # (the n+1 accounts for the possibility that all clusters are unique)
            a = which.min(1:(n+1) %in% clus.w)
            b = (n+1) - which.max((n+1):1 %in% clus.w) + 1
            if (b > a)
                clus.w[clus.w == b] = a

            }

#       if ((floor(iter/window) == iter/window) && iter <= nburn){
#           cand.sig[[i]] = (cand.sig[[i]] +
#               autotune(mean(accept[(iter-window+1):iter,i]), k = max(window/50,1.5)) *
#               cov(param.theta[(iter-window+1):iter, ((i-1)*nparam+1):(i*nparam)]))/2
#           }
        }
    n.j = table(clus.w)
    temp.label = as.numeric(names(n.j))
    n.star = length(temp.label)
    temp.star = double(n.star)
    for (j in 1:n.star)
        temp.star[j] = which(clus.w == temp.label[j])[1]
    theta.star = matrix(c.theta[temp.star,], ncol = nparam)

    # Update tau (measurement variance)
    sumsq = double(n)
    for (i in 1:n)
        sumsq[i] = sum((dat[[i]]$xy[,2]-jc(dat[[i]]$fixed, c.theta[i,]))^2)
    c.tau = rinvgamma(1, prior.tau.a + 1/2 * sum(ni), prior.tau.b + 1/2 * sum(sumsq))

    # Update alpha
    eta = rbeta(1, c.alpha + 1, n)
    eps = (prior.alpha.a + n.star - 1) /
        (n*(prior.alpha.b - log(eta)) + prior.alpha.a + n.star - 1)
    if (runif(1) < eps){
        c.alpha = rgamma(1, prior.alpha.a + n.star, prior.alpha.b - log(eta))
    } else {
        c.alpha = rgamma(1, prior.alpha.a + n.star - 1, prior.alpha.b - log(eta))
        }

    # Update sigma (baseline covariance)
    temp.V = prior.sigma.V
    for (j in 1:n.star)
        temp.V = temp.V + (theta.star[j,] - c.mu) %*% t(theta.star[j,] - c.mu)
    inv.c.sigma = rwishart(solve(temp.V), prior.sigma.df + n.star)
    c.sigma = solve(inv.c.sigma)

    # Update mu (baseline mean)
    temp.s = apply(theta.star, 2, sum)
    temp.G = solve(solve(prior.mu.B) + n.star * inv.c.sigma)
    temp.mu = temp.G %*% (solve(prior.mu.B) %*% prior.mu.a + inv.c.sigma %*% temp.s)
#   c.mu = mvrnorm(1, temp.mu, temp.G)
    c.mu = rtmvnorm(1, as.numeric(temp.mu), temp.G, lower = theta.lower, upper = theta.upper,
        algorithm = "gibbs", burn.in.samples = 50, start.value = c.mu)
    c.mu = as.numeric(c.mu)

    # Improve mixing on the unique thetas
    for (j in 1:n.star){
        temp.w = which(clus.w == j)
        cand = mvrnorm(1, theta.star[j,], cand.sig)

        candsumsq = 0
        for (i in temp.w)
            candsumsq = candsumsq + sum((dat[[i]]$xy[,2]-jc(dat[[i]]$fixed, cand))^2)
        currsumsq = 0
        for (i in temp.w)
            currsumsq = currsumsq + sum((dat[[i]]$xy[,2]-jc(dat[[i]]$fixed, c.theta[i,]))^2)

        temp.cand = dmvnorm(cand, c.mu, c.sigma) - 1/(2*c.tau) * candsumsq
        temp.curr = dmvnorm(theta.star[j,], c.mu, c.sigma) - 1/(2*c.tau) * currsumsq
        if (log(runif(1)) < temp.cand - temp.curr)
            c.theta[temp.w,] = matrix(cand, length(temp.w), nparam, byrow = TRUE)
        }

    # Reset the objects
    param.theta[iter,] = c(t(c.theta))
    param.alpha[iter] = c.alpha
    param.mu[iter,] = c.mu
    param.sigma[[iter]] = c.sigma
    param.tau[iter] = c.tau

    if (iter == (nburn + nmcmc))
        cat("\n")
    }






### Remove burnin
param.theta = tail(param.theta, nmcmc)
param.alpha = tail(param.alpha, nmcmc)
param.mu = tail(param.mu, nmcmc)
param.sigma = tail(param.sigma, nmcmc)
param.tau = tail(param.tau, nmcmc)

#accept = tail(accept, nmcmc)
#apply(accept, 2, mean)

### Stack theta matrix
new.theta = matrix(0, nmcmc*n, nparam)
for (i in 1:n)
    new.theta[((i-1)*nmcmc + 1):(i*nmcmc),] = param.theta[,((i-1)*nparam+1):(i*nparam)]

### Some plots of the posteriors
pairs(param.mu, pch = 20)
plot(sapply(param.sigma, function(x) determinant(x)$modulus[1]), type='l')

plot(param.tau, type='l')
plot(param.alpha, type='l')


### May take a while (mean of theta overlain by individual thetas)
pairs(rbind(param.mu, new.theta), pch = 20, col = c(rep("black", nmcmc),
    rep(cols, each = nmcmc)), cex = 0.5)



plot(0, type='n', xlim = c(1, nmcmc), ylim = range(param.theta))
for (i in 1:n)
    matplot(param.theta[,seq((i-1)*nparam+1, i*nparam)], lty = i, type = 'l', add = TRUE)


### Predictions of the observations
for (i in 1:n){
    pred = matrix(0, nmcmc, nrow(dat[[i]]$xy))
    for (j in 1:nmcmc)
        pred[j,] = jc(dat[[i]]$fixed,
            param.theta[j, seq((i-1)*nparam+1, i*nparam)]) +
            rnorm(ni[i], 0, sqrt(param.tau[j]))
    qlines = apply(pred, 2, quantile, c(0.025, 0.975))
#   matplot(dat[[i]]$xy[,1], t(pred), type='l', lty = 1, col = 'steelblue', lwd = 0.5)
    plot(dat[[i]]$xy, lwd = 3, col = cols[i], ylim = range(pred), type = 'l')
    lines(dat[[i]]$xy[,1], qlines[1,], lwd = 3, col = col.mult("gray50", cols[i]))
    lines(dat[[i]]$xy[,1], qlines[2,], lwd = 3, col = col.mult("gray50", cols[i]))
    if (i != n)
        readline()
    }

### Predict a new observation
pred.theta_0 = matrix(0, nmcmc, nparam)
for (i in 1:nmcmc){
    pred.theta_0
    prob = c(param.alpha[i] / (param.alpha[i] + n), rep(1 / (param.alpha[i] + n), n))
    draw = sample(n+1, 1, prob = prob)
    if (draw == 1){
        pred.theta_0[i,] = rtmvnorm(1, param.mu[i,], param.sigma[[i]],
            lower = theta.lower, upper = theta.upper, algorithm = "gibbs",
            burn.in.sample = 50, start.value = param.mu[i,])
    } else {
        pred.theta_0[i,] = param.theta[i, (((draw-1)-1)*nparam + 1):((draw-1)*nparam)]
        }
    }
pairs(pred.theta_0, pch = 20, cex = 0.5)

pred.y_0 = NULL
for (i in 1:n){
    pred.y_0[[i]] = matrix(0, nmcmc, ni[i])
    for (j in 1:nmcmc){
        pred.y_0[[i]][j,] = jc(dat[[i]]$fixed, pred.theta_0[j,]) +
            rnorm(ni[i], 0, sqrt(param.tau[j]))
        }
    }

mpred = lapply(pred.y_0, function(x) apply(x, 2, mean))
qpred = lapply(pred.y_0, function(x) apply(x, 2, quantile, c(0.025, 0.975)))



plot(0, type='n', xlim = c(0, 0.27), ylim = c(-3, 3), xlab = "Plastic Strain",
    ylab = "Stress")
for (i in 1:n){
    lines(dat[[i]]$xy[,1], dat[[i]]$xy[,2], type='l', col = cols[i], lwd = 2)
    lines(dat[[i]]$xy[,1], mpred[[i]], type='l', col = darkcols[i], lwd = 2, lty = 3)
#   lines(dat[[i]]$xy[,1], qpred[[i]][1,], type='l', col = darkcols[i], lwd = 2)
#   lines(dat[[i]]$xy[,1], qpred[[i]][2,], type='l', col = darkcols[i], lwd = 2)
    if (i != n)
        readline()
    }


