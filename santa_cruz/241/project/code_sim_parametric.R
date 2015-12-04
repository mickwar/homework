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

# G0 (baseline) mean
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
param.mu = matrix(0, nburn + nmcmc, nparam)                         # G0 mean
param.sigma = rep(list(matrix(0, nparam, nparam)), nburn + nmcmc)   # G0 variance
param.tau = double(nburn + nmcmc)                                   # Measurement variance

#cand.sig = diag(1e-6, nparam * n)
cand.sig = rep(list(diag(1e-2, nparam)), n)
accept = matrix(0, nburn + nmcmc, n)
window = 200


param.mu[1,] = prior.mu.a
param.sigma[[1]] = diag(1, nparam)
param.tau[1] = 1
param.theta[1,] = rep(param.mu[1,], n)


for (iter in 2:(nburn + nmcmc)){
    cat("\r", iter, "/", nburn + nmcmc)

    # Current values
    c.theta = matrix(param.theta[iter-1,], n, nparam, byrow = TRUE)
#   c.theta = param.theta[iter-1,]
    c.mu = param.mu[iter-1,]
    c.sigma = param.sigma[[iter-1]]
    c.tau = param.tau[iter-1]

    # Update the thetas
    inv.c.sigma = solve(c.sigma)
    currsumsq = double(n)
    currtrace = double(n)
    for (i in 1:n){
        currsumsq[i] = sum((y[i,]-jc(x,c.theta[i,]))^2)
        currtrace[i] = t((c.theta[i,] - c.mu)) %*% inv.c.sigma %*% (c.theta[i,] - c.mu)
        }

    for (i in 1:n){
        cand = mvrnorm(1, c.theta[i,], cand.sig[[i]])

        candsumsq = sum((y[i,]-jc(x,cand))^2)
        candtrace = t((cand - c.mu)) %*% inv.c.sigma %*% (cand - c.mu)

        temp.cand = -1/(2*c.tau) * candsumsq - 1/2 * candtrace
        temp.curr = -1/(2*c.tau) * currsumsq[i] - 1/2 * currtrace[i]
        if (log(runif(1)) < temp.cand - temp.curr){
            c.theta[i,] = cand
            accept[iter,i] = 1
            currsumsq[i] = candsumsq
            currtrace[i] = candtrace
            }

        if ((floor(iter/window) == iter/window) && iter <= nburn){
            cand.sig[[i]] = (cand.sig[[i]] +
                autotune(mean(accept[(iter-window+1):iter,i]), k = max(window/50,1.5)) *
                cov(param.theta[(iter-window+1):iter,((i-1)*nparam+1):(i*nparam)]))/2
            }
        }


    # Update tau (measurement variance)
    c.tau = rinvgamma(1, prior.tau.a + 1/2 * sum(ni), prior.tau.b + 1/2 * sum(currsumsq))


    # Update sigma (superpopulation covariance
    temp.V = prior.sigma.V
    for (j in 1:n)
        temp.V = temp.V + (c.theta[i,] - c.mu) %*% t(c.theta[i,] - c.mu)
    inv.c.sigma = rwishart(solve(temp.V), prior.sigma.df + n)
    c.sigma = solve(inv.c.sigma)

    # Update mu (baseline mean)
    temp.s = apply(c.theta, 2, sum)
    temp.G = solve(solve(prior.mu.B) + n * inv.c.sigma)
    temp.mu = temp.G %*% (solve(prior.mu.B) %*% prior.mu.a + inv.c.sigma %*% temp.s)
    c.mu = mvrnorm(1, temp.mu, temp.G)

    # Reset the objects
    param.theta[iter,] = c(t(c.theta))
    param.mu[iter,] = c.mu
    param.sigma[[iter]] = c.sigma
    param.tau[iter] = c.tau

    if (iter == (nburn + nmcmc))
        cat("\n")
    }





### Remove burnin
param.theta = tail(param.theta, nmcmc)
param.mu = tail(param.mu, nmcmc)
param.sigma = tail(param.sigma, nmcmc)
param.tau = tail(param.tau, nmcmc)

accept = tail(accept, nmcmc)
apply(accept, 2, mean)


thin = seq(100, nmcmc, by = 100)
thin.theta = param.theta[thin,]
thin.mu = param.mu[thin,]
thin.tau = param.tau[thin]
thin.accept = accept[thin,]

old = thin.sigma
thin.sigma = rep(list(matrix(0, nparam, nparam)), length(thin))
for (i in 1:length(thin))
    thin.sigma[[i]] = old[[thin[i]]]

save(thin.theta, thin.sigma, thin.mu, thin.tau, thin.accept, file = "./workspaces/toy_para.RData")


### Some plots of the posteriors
pairs(param.mu, pch = 20)

plot(sapply(param.sigma, function(x) determinant(x)$modulus[1]), type='l')

plot(param.tau, type='l')



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
for (i in 1:nmcmc)
    pred.theta_0[i,] = mvrnorm(1, param.mu[i,], param.sigma[[i]])
pairs(pred.theta_0, pch = 20, cex = 0.5)

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

#matplot(pred.x, t(pred.y_0), type='l', lty = 1, lwd = 0.5, col = 'steelblue')
plot(0, type='n', xlim = range(x), ylim = range(qhpd, na.rm = TRUE))
matplot(x, t(y), type='l', lty = 1, lwd = 0.5, add = TRUE, col = 'gray20')
#for (i in 1:length(pred.x))
#    points(rep(pred.x[i], length(qhpd[[i]])), qhpd[[i]], col = 'blue', pch = 20)
#for (i in 1:4)
#    lines(pred.x, sapply(qhpd, function(x) x[i]), col = 'blue', lwd = 2)
#lines(pred.x, qlines[2,], col = 'blue')

lines(pred.x, qhpd[1,], col = 'darkblue', lwd = 3)
lines(pred.x, qhpd[2,], col = 'darkblue', lwd = 3)


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
