source("~/files/R/mcmc/bayes_functions.R")
library(MASS)

### Load the data (use load("./chen_gray_fake.RData") when real data not available)
#source("./read_data.R")
#dat = read_data()

 set.seed(1)
 n = 100
 b0 = rnorm(n, 3, 0.3)
 b1 = rnorm(n, 1.5, 0.1)
 truth = 1*(runif(n) < 0.5)
 b2 = ifelse(truth, rnorm(n, 0, 0.0), rnorm(n, -1, 0.1))
 plot(density(b2))

#set.seed(1)
#n = 100
#b0 = rnorm(n, 3, 0.1)
#b1 = rnorm(n, 1.5, 0.3)
#b0 = rep(3, n)
#b1 = rep(1.5, n)
#truth = 1*(runif(n) < 0.5)
#b2 = ifelse(truth, rep(0, n), rep(-1, n))

x = seq(-0.5, 2, length = 10)
y = matrix(0, n, length(x))
for (i in 1:n)
    y[i,] = b0[i] + b1[i]*x + b2[i]*x^2 + rnorm(length(x), 0, 0.05)

cols = rainbow(n)
darkcols = character(n)
for (i in 1:n)
    darkcols[i] = col.mult(cols[i], "gray50")
matplot(x, t(y), type='l', lty = 1, lwd = 0.7, col = cols)

dat = rep(list(NULL), n)
for (i in 1:n)
    dat[[i]] = list("xy"=cbind(x, y[i,]), "fixed"=x)

ni = sapply(dat, function(x) nrow(x$xy))


### Called jc still cause lazy
jc = function(x, theta)
    theta[1]*x^0 + theta[2]*x^1 + theta[3]*x^2

### Multivariate normal density (up to a constant)
dmvnorm = function(x, mu, sigma, log = TRUE){
    if (NROW(sigma) == 1){
        p = length(x)
#       out = 0.5*p*log(sigma) - 0.5/sigma * t(x-mu) %*% (x-mu)
        out = -0.5*p*log(sigma) - 0.5/sigma * sum((x-mu)^2)
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
prior.tau.a = 5
prior.tau.b = 1

# DP precision
prior.alpha.a = 5
prior.alpha.b = 1

# G0 (baselin) mean
prior.mu.a = c(0, 0, 0)
prior.mu.B = diag(10, 3)

# G0 (basline) covariance
prior.sigma.df = 3
prior.sigma.V = diag(1, 3)




### MCMC
nburn = 50000
nmcmc = 1000000

nparam = 3 # Parameters in the polynomial
param.theta = matrix(0, nburn + nmcmc, nparam * n)                  # JC parameters
param.alpha = double(nburn + nmcmc)                                 # DP precision
param.mu = matrix(0, nburn + nmcmc, nparam)                         # G0 mean
param.sigma = rep(list(matrix(0, nparam, nparam)), nburn + nmcmc)   # G0 variance
param.tau = double(nburn + nmcmc)                                   # Measurement variance

# cand.sig = rep(list(0.1*diag(nparam)), n)
# accept = matrix(0, nburn + nmcmc, n)
# window = 200
cand.sig = 0.001 * diag(nparam)


param.alpha[1] = 1
param.mu[1,] = c(0,0,0)
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
            cand = mvrnorm(1, c.mu, c.sigma)
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
    c.mu = mvrnorm(1, temp.mu, temp.G)

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


thin = seq(100, nmcmc, by = 100)
thin.theta = param.theta[thin,]
thin.alpha = param.alpha[thin]
thin.mu = param.mu[thin,]
thin.tau = param.tau[thin]

old = thin.sigma
thin.sigma = rep(list(matrix(0, nparam, nparam)), length(thin))
for (i in 1:length(thin))
    thin.sigma[[i]] = old[[thin[i]]]

save(thin.theta, thin.alpha, thin.sigma, thin.mu, thin.tau, file = "./workspaces/toy_semi.RData")


#accept = tail(accept, nmcmc)
#apply(accept, 2, mean)

### Stack theta matrix
#new.theta = matrix(0, nmcmc*n, nparam)
#for (i in 1:n)
#    new.theta[((i-1)*nmcmc + 1):(i*nmcmc),] = param.theta[,((i-1)*nparam+1):(i*nparam)]

### Some plots of the posteriors
pairs(param.mu, pch = 20)
plot(sapply(param.sigma, function(x) determinant(x)$modulus[1]), type='l')

plot(param.tau, type='l')
plot(param.alpha, type='l')

#plot(param.theta[,3], type='l')

### May take a while (mean of theta overlain by individual thetas)
cols2 = character(n)
cols2[truth == 1] = "green"
cols2[truth == 0] = "red"
m.theta = matrix(apply(param.theta, 2, mean), n, nparam, byrow = TRUE)
pairs(rbind(param.mu, m.theta), pch = 20, col = c(rep("black", nmcmc), cols2), cex = 1.0)



plot(0, type='n', xlim = c(1, nmcmc), ylim = range(param.theta))
for (i in 1:n)
    matplot(param.theta[,seq((i-1)*nparam+1, i*nparam)], lty = i, type = 'l', add = TRUE)


for (i in 1:n){
    pred = matrix(0, nmcmc, nrow(dat[[i]]$xy))
    for (j in 1:nmcmc)
        pred[j,] = jc(dat[[i]]$fixed,
            param.theta[j, seq((i-1)*nparam+1, i*nparam)]) + 
            rnorm(1, 0, sqrt(param.tau[j]))
    qlines = apply(pred, 2, quantile, c(0.025, 0.975))
#   matplot(dat[[i]]$xy[,1], t(pred), type='l', lty = 1, col = 'steelblue', lwd = 0.5)
#   lines(dat[[i]]$xy, lwd = 3, col = cols[i])
    plot(dat[[i]]$xy, lwd = 3, col = cols[i], ylim = range(pred), type='l')
    lines(dat[[i]]$xy[,1], qlines[1,], lwd = 3, col = darkcols[i])
    lines(dat[[i]]$xy[,1], qlines[2,], lwd = 3, col = darkcols[i])
    if (i != n)
        readline()
    }


pred.theta_0 = matrix(0, nmcmc, nparam)
for (i in 1:nmcmc){
    prob = c(param.alpha[i] / (param.alpha[i] + n), rep(1/(param.alpha[i] + n), n))
    draw = sample(n+1, 1, replace = FALSE, prob = prob)
    if (draw == 1){
        pred.theta_0[i,] = mvrnorm(1, param.mu[i,], param.sigma[[i]])
    } else {
        pred.theta_0[i,] = param.theta[i, (((draw-1)-1)*nparam + 1):((draw-1)*nparam)]
        }
    }
pairs(pred.theta_0, pch = 20)

pred.x = seq(min(x), max(x), length = 50)
#pred.x = x
pred.y_0 = t(apply(pred.theta_0, 1,
    function(y) y[1]*pred.x^0 + y[2]*pred.x^1 + y[3]*pred.x^2))
for (i in 1:nmcmc)
    pred.y_0[i,] = pred.y_0[i,] + rnorm(1, 0, sqrt(param.tau[i]))
#qlines = apply(pred.y_0, 2, quantile, c(0.025, 0.975))
qhpd = apply(pred.y_0, 2, function(x) hpd.mult(x, density(x), 0.90, 1000))
for (i in 1:length(pred.x))
    if (length(qhpd[[i]]) == 2)
        qhpd[[i]] = c(qhpd[[i]][1], NA, NA, qhpd[[i]][2])

matplot(pred.x, t(pred.y_0), type='l', lty = 1, lwd = 0.5, col = 'steelblue')
#plot(0, type='n', xlim = range(x), ylim = range(qhpd, na.rm = TRUE))
matplot(x, t(y), type='l', lty = 1, lwd = 1.0, add = TRUE, col = 'gray20')
#for (i in 1:length(pred.x))
#    points(rep(pred.x[i], length(qhpd[[i]])), qhpd[[i]], col = 'blue', pch = 20)
for (i in 1:4)
    lines(pred.x, sapply(qhpd, function(x) x[i]), col = 'blue', lwd = 2)
#lines(pred.x, qlines[2,], col = 'blue')


#h = t(apply(pred.theta_0, 1, function(y) y[1]*x^0 + y[2]*x^1 + y[3]*x^2))
#
#
#param.theta[i, (((draw-1)-1)*nparam + 1):((draw-1)*nparam)]
#
#dmvnorm(
#y[1,]
#dim(y)
#for (i in 1:nmcmc)
#    pred.y_0[i,] = pred.y_0[i,] + rnorm(1, 0, sqrt(param.tau[i]))
