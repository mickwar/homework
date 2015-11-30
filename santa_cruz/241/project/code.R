source("~/files/R/mcmc/bayes_functions.R")
library(MASS)

### Load the data (use load("./chen_gray_fake.RData") when real data not available)
source("./read_data.R")
dat = read_data()
n = length(dat)

ni = sapply(dat, function(x) nrow(x$xy))

# Plot the data
cols = rainbow(n)
plot(0, type='n', xlim = c(0, 0.27), ylim = c(0, 0.012), xlab = "Plastic Strain",
    ylab = "Stress")
for (i in 1:n)
    lines(dat[[i]]$xy[,1], dat[[i]]$xy[,2], type='l', col = cols[i], lwd = 2)

### Simplified Johnson-Cook model
jc = function(eps, t, eta){
    A = eta[1]
    B = eta[2]
    n = eta[3]
    C = eta[4]
    m = eta[5]
    strain_rate = t[1]
    temperature = t[2] / 3250 # ~melting (kelvin)
    return ((A + B * eps^n) * (1 + C * log(strain_rate)) * (1 - temperature^m))
    }

### Multivariate normal density (up to a constant)
dmvnorm = function(x, mu, sigma, log = TRUE){
    if (NROW(sigma) == 1){
        p = length(x)
#       out = 0.5*p*log(sigma) - 0.5/sigma * t(x-mu) %*% (x-mu)
        out = 0.5*p*log(sigma) - 0.5/sigma * sum((x-mu)^2)
    } else {
        inv = solve(sigma)
        out = 0.5*determinant(inv)$modulus[1] - 0.5*t(x-mu) %*% inv %*% (x-mu)
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
prior.alpha.a = 100
prior.alpha.b = 5

# G0 (baselin) mean
prior.mu.a = c(1,1,0,1,0)
prior.mu.B = diag(5)

# G0 (basline) covariance
prior.sigma.df = 5
prior.sigma.V = diag(1, 5)



### MCMC
nburn = 10000
nmcmc = 20000

nparam = 5 # Parameters in JC model
param.theta = matrix(0, nburn + nmcmc, nparam * n)                  # JC parameters
param.alpha = double(nburn + nmcmc)                                 # DP precision
param.mu = matrix(0, nburn + nmcmc, nparam)                         # G0 mean
param.sigma = rep(list(matrix(0, nparam, nparam)), nburn + nmcmc)   # G0 variance
param.tau = double(nburn + nmcmc)                                   # Measurement variance

cand.sig = rep(list(0.1*diag(nparam)), n)
accept = matrix(0, nburn + nmcmc, n)
window = 200


param.alpha[1] = 1
param.mu[1,] = c(1,1,0,1,0)
param.sigma[[1]] = diag(1, nparam)
param.tau[1] = 1
param.theta[1,] = rep(param.mu[1,], n)
clus.w = rep(1, n)


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
        n.j.minus = table(clus.w)
        temp.label = as.numeric(names(n.j.minus))
        temp.ind = which(temp.label == clus.w[i])
        n.j.minus[temp.ind] = n.j.minus[temp.ind] - 1

        n.star.minus = length(temp.label) # but not really
        temp.star = double(n.star.minus)
        for (j in 1:n.star.minus)
            temp.star[j] = which(clus.w == temp.label[j])[1]
        theta.star.minus = matrix(c.theta[temp.star,], ncol = nparam)

        A = c.alpha / (n - 1 + c.alpha)
        Bj = n.j.minus / (n - 1 + c.alpha)

        # Draw a candidate value from prior conditional
        temp.samp = sample(n.star.minus + 1, 1, prob = c(A, Bj))
        if (temp.samp == 1){ # Make a draw from G0
#           cand = mvrnorm(1, c.mu, c.sigma)
            cand = mvrnorm(1, c.theta[i,], cand.sig[[i]])
        } else { # Make a draw from the existing groups
            cand = theta.star.minus[temp.samp-1,]
            }

        # Decide whether to accept the candidate
        temp.cand = dmvnorm(dat[[i]]$xy[,2],
            jc(dat[[i]]$xy[,1], c(dat[[i]]$strain_rate, dat[[i]]$temperature), cand),
            c.tau)
        temp.curr = dmvnorm(dat[[i]]$xy[,2],
            jc(dat[[i]]$xy[,1], c(dat[[i]]$strain_rate, dat[[i]]$temperature), c.theta[i,]),
            c.tau)
        if (temp.samp != 1){
            temp.cand = temp.cand +
                dmvnorm(cand, c.mu, c.sigma) - 0
#               dmvnorm(cand, c.theta[i,], cand.sig[[i]])
            temp.curr = temp.curr +
                dmvnorm(c.theta[i,], c.mu, c.sigma) - 0
#               dmvnorm(c.theta[i,], cand, cand.sig[[i]])
        } else {
            temp.cand = temp.cand + 
            }
#           temp.cand = temp.cand - dmvnorm(cand, c.theta[i,], cand.sig[[i]]) - log(A)
        temp.cand - temp.curr
        if (log(runif(1)) < temp.cand - temp.curr){
            c.theta[i,] = cand
            # Acceptance rate? Neal doesn't seem to mention it
            accept[iter,i] = 1

            # Change clus.w
            if (temp.samp == 1){ # a new cluster was formed
                clus.w[i] = which.min(1:n %in% clus.w[-i])
            } else {            # part of existing cluster
                clus.w[i] = temp.samp - 1
                }

            # Re-arrange clus.w for empty clusters
            a = which.min(1:n %in% clus.w)
            b = n - which.max(n:1 %in% clus.w) + 1
            if (b > a)
                clus.w[clus.w == b] = a

            }

        if ((floor(iter/window) == iter/window) && iter <= nburn){
            cand.sig[[i]] = (cand.sig[[i]] +
                autotune(mean(accept[(iter-window+1):iter,i]), k = max(window/50,1.5)) *
                cov(param.theta[(iter-window+1):iter, ((i-1)*nparam+1):(i*nparam)]))/2
            }
        }
    n.j = table(clus.w)
    temp.label = as.numeric(names(n.j))
    n.star = length(temp.label)
    temp.star = double(n.star)
    for (j in 1:n.star)
        temp.star[j] = which(clus.w == temp.label[j])[1]
    theta.star = matrix(c.theta[temp.star,], ncol = nparam)

    # Update tau (measurement variance)
    sumsq = 0
    for (i in 1:n)
        sumsq = sumsq + sum((dat[[i]]$xy[,2]-jc(dat[[i]]$xy[,1],
            c(dat[[i]]$strain_rate, dat[[i]]$temperature), c.theta[i,]))^2)
    c.tau = rinvgamma(1, prior.tau.a + 1/2 * sum(ni), prior.tau.b + 1/2 * sumsq)

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
    c.mu = mvrnorm(1, temp.mu, temp.G)

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

accept = tail(accept, nmcmc)
apply(accept, 2, mean)

### Stack theta matrix
new.theta = matrix(0, nmcmc*n, nparam)
for (i in 1:n)
    new.theta[((i-1)*nmcmc + 1):(i*nmcmc),] = param.theta[,((i-1)*nparam+1):(i*nparam)]

### Some plots of the posteriors
pairs(param.mu, pch = 20)
plot(sapply(param.sigma, function(x) determinant(x)$modulus[1]), type='l')

plot(param.tau, type='l')
plot(param.alpha, type='l')

plot(param.theta[,2], type='l')

### May take a while (mean of theta overlain by individual thetas)
pairs(rbind(param.mu, new.theta), pch = 20, col = c(rep("black", nmcmc),rep(cols, each = nmcmc)))



plot(0, type='n', xlim = c(1, nmcmc), ylim = range(param.theta))
for (i in 1:n)
    matplot(param.theta[,seq((i-1)*nparam+1, i*nparam)], lty = i, type = 'l', add = TRUE)


for (i in 1:n){
    pred = matrix(0, nmcmc, nrow(dat[[i]]$xy))
    for (j in 1:nmcmc)
        pred[j,] = jc(dat[[i]]$xy[,1], c(dat[[i]]$strain_rate, dat[[i]]$temperature),
            param.theta[j, seq((i-1)*nparam+1, i*nparam)]) + 0
#           rnorm(nrow(dat[[i]]$xy), 0, sqrt(param.tau[j]))

    matplot(dat[[i]]$xy[,1], t(pred), type='l', lty = 1, col = 'steelblue', lwd = 0.5)
    lines(dat[[i]]$xy, lwd = 3, col = cols[i])
    readline()
    }
