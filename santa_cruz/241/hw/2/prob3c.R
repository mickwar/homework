library(MASS)
source("~/files/R/mcmc/bayes_functions.R")

### Posterior predictive loss criterion
pplc = function(y, ypred, k = Inf){
    n = length(y)
    vars = apply(ypred, 2, var)
    means = apply(ypred, 2, mean)

    factor = k / (k + 1)
    if (k == Inf)
        factor = 1

    return (sum(vars) + factor*sum((y - means)^2))
    }

dat = read.table("~/files/data/fabric.txt", header = TRUE)
 x = dat$length / 100
#x = as.numeric(scale(dat$length))
y = dat$faults
n = length(y)
ord = order(y)

set.seed(1)
jity = jitter(y)

par(mfrow = c(1,1), mar = c(5.1, 4.1, 4.1, 2.1))
plot(x, y, pch = 20)


### Inverse gamma functions
# Mean = rate / (shape - 1)
# What I'm calling rate for the inverse gamma, Wikipedia calls scale
dinvgamma = function(x, shape, rate, log = FALSE){
    out = shape * log(rate) - lgamma(shape) - (shape + 1) * log(x) - rate / x
    if (log)
        return (out)
    return (exp(out))
    }
rinvgamma = function(n, shape, rate)
    1/rgamma(n, shape = shape, rate = rate)


### MCMC
nburn = 20000
nmcmc = 100000

# Parameter objects
param.theta = matrix(0, nburn + nmcmc, n)
param.alpha = double(nburn + nmcmc)
param.psi = matrix(0, nburn + nmcmc, 2) # [,1] is zeta, [,2] is mu
param.beta = double(nburn + nmcmc)

# other mcmc
sig.beta = 0.01
acc.beta = double(nburn + nmcmc)
sig.psi = diag(0.01, 2)
acc.psi = double(nburn + nmcmc)

window = 200

# Priors
prior.beta.a  = 0   # Normal (mean)
prior.beta.b  = 1   # Normal (sd)
prior.alpha.a = 1   # Gamma (shape)           alpha controls n.star (discreteness of G)
prior.alpha.b = 1/2 # Gamma (rate)
#prior.zeta.a  = 2   # Gamma (shape)
#prior.zeta.b  = 1/4 # Gamma (rate)
#prior.mu.a    = 8   # Gamma (shape)
#prior.mu.b    = 1/2 # Gamma (rate)
 prior.zeta.a  = 1   # Gamma (shape)
 prior.zeta.b  = 1/2 # Gamma (rate)
 prior.mu.a    = 1   # Gamma (shape)
 prior.mu.b    = 1/2 # Gamma (rate)

#p1 = rgamma(50000, prior.zeta.a, prior.zeta.b)
#p2 = rgamma(50000, prior.mu.a, prior.mu.b)
#p3 = rgamma(50000, p1, p1 / p2)
#
#mean(p1); sd(p1)
#mean(p2); sd(p2)
#mean(p3); sd(p3)
#
#plot(density(p1))
#plot(density(p2))
#plot(density(p3))


# Initial values
#param.theta[1,] = y
param.beta[1] = rnorm(1, prior.beta.a, prior.beta.b)
param.alpha[1] = rgamma(1, prior.alpha.a, prior.alpha.b)
param.psi[1,] = c(rgamma(1, prior.zeta.a, prior.zeta.b),
    rgamma(1, prior.mu.a, prior.mu.b))
#param.theta[1,] = rep(0, n)
param.theta[1,] = y


# Iterations
for (iter in 2:(nburn + nmcmc)){
    cat("\r", iter, "/", nburn + nmcmc)

    # Current parameters (used for convenience)
    c.theta = param.theta[iter-1,]
    c.alpha = param.alpha[iter-1]
     c.beta = param.beta[iter-1]
     c.zeta = param.psi[iter-1,1]
       c.mu = param.psi[iter-1,2]

    # Update thetas
    for (i in 1:n){
        n.j.minus = table(c.theta[-i])
        theta.star.minus = as.numeric(names(n.j.minus))
        n.star.minus = length(theta.star.minus)

        temp.p = (c.zeta / c.mu) / (exp(c.beta*x[i]) + c.zeta / c.mu)
        q0 = dnbinom(y[i], size = c.zeta, prob = temp.p)

        qj = dpois(y[i], theta.star.minus*exp(c.beta*x[i]))

        # Probabilities of determining which to draw
        # Calculate A
        A = c.alpha * q0 / (c.alpha * q0 + sum(n.j.minus * qj))
        
        # Calculate B's
        Bj = n.j.minus * qj / (c.alpha * q0 + sum(n.j.minus * qj))

        # Make the update
        draw = sample(n.star.minus + 1, 1, prob = c(A, Bj))
        if (draw == 1){ # Make a draw from h
            c.theta[i] = rgamma(1, y[i] + c.zeta, exp(c.beta*x[i]) + c.zeta/c.mu)
        } else { # Make a draw from the existing groups
            c.theta[i] = theta.star.minus[draw-1]
            }
        }
    n.j = table(c.theta)
    theta.star = as.numeric(names(n.j))
    n.star = length(theta.star)

    # Update beta
    cand = rnorm(1, c.beta, sig.beta)
    temp.curr = dnorm(c.beta, prior.beta.a, prior.beta.b, log = TRUE) +
        sum(dpois(y, c.theta*exp(c.beta*x), log = TRUE))
    temp.cand = dnorm(cand, prior.beta.a, prior.beta.b, log = TRUE) +
        sum(dpois(y, c.theta*exp(cand*x), log = TRUE))
    if (log(runif(1)) < temp.cand - temp.curr){
        c.beta = cand
        acc.beta[iter] = 1
        }

    # Update psi = (zeta, mu)
    cand = mvrnorm(1, c(c.zeta, c.mu), sig.psi)
    if (all(cand > 0)){
        temp.curr = dgamma(c.zeta, prior.zeta.a, prior.zeta.b, log = TRUE) +
            dgamma(c.mu, prior.mu.a, prior.mu.b, log = TRUE) +
            sum(dgamma(theta.star, c.zeta, c.zeta / c.mu, log = TRUE))
        temp.cand = dgamma(cand[1], prior.zeta.a, prior.zeta.b, log = TRUE) +
            dgamma(cand[2], prior.mu.a, prior.mu.b, log = TRUE) +
            sum(dgamma(theta.star, cand[1], cand[1] / cand[2], log = TRUE))
        if (log(runif(1)) < temp.cand - temp.curr){
            c.zeta = cand[1]
            c.mu = cand[2]
            acc.psi[iter] = 1
            }
        }

    # Update alpha
    # Use an auxiliary variable eta to draw a new alpha
    eta = rbeta(1, c.alpha + 1, n)
    eps = (prior.alpha.a + n.star - 1) /
        (n*(prior.alpha.b - log(eta)) + prior.alpha.a + n.star - 1)
    if (runif(1) < eps){
        c.alpha = rgamma(1, prior.alpha.a + n.star, prior.alpha.b - log(eta))
    } else {
        c.alpha = rgamma(1, prior.alpha.a + n.star - 1, prior.alpha.b - log(eta))
        }

    # Extra improvement step (only updating the clusters, not individuals)
    for (j in 1:n.star){
        temp.w = which(c.theta == theta.star[j])
        c.theta[temp.w] = rgamma(1,
            c.zeta + sum(y[temp.w]),
            c.zeta / c.mu + sum(exp(c.beta*x[temp.w])))
        }

    # Put the c.* objects into the regular ones
    param.theta[iter,] = c.theta
    param.alpha[iter] = c.alpha
    param.beta[iter] = c.beta
    param.psi[iter,1] = c.zeta
    param.psi[iter,2] = c.mu

    # Improve candidate sigmas
    if ((floor(iter/window) == iter/window) && (iter <= nburn)){
        sig.beta = autotune(mean(acc.beta[(iter-window+1):iter]),
            target = 0.25, k = window/50) * sig.beta
        sig.psi = autotune(mean(acc.psi[(iter-window+1):iter]),
            target = 0.25, k = window/50) *
            (sig.psi + window * var(param.psi[(iter-window+1):iter,]) / iter)
        }

    if (iter == (nburn + nmcmc))
        cat("\n")
    }

### Burn in
param.theta = tail(param.theta, nmcmc)
param.alpha = tail(param.alpha, nmcmc)
param.beta = tail(param.beta, nmcmc)
param.psi = tail(param.psi, nmcmc)

acc.beta = tail(acc.beta, nmcmc)
acc.psi = tail(acc.psi, nmcmc)

mean(acc.beta)
mean(acc.psi)



### Individual plots for the theta's
#par(mfrow = c(4, 2), mar = c(3.1, 2.1, 2.1, 1.1))
#for (k in 1:n){
#    plot(param.theta[,k], type = 'l')
#    abline(h = y[k], col = 'red', lty = 2)
#    plot(density(param.theta[,k]), main = k)
#    abline(v = y[k], col = 'red', lty = 2)
#    if (k %% 4 == 0)
#        readline()
#    }

### "Trace plot" for clusters
par(mfrow = c(1,1), mar = c(3.1, 4.1, 2.1, 1.1))
plot(0, type='n', xlim = c(1, nrow(param.theta)), ylim = range(param.theta),
    xlab = "MCMC Iteration", ylab = "Unique cluster locations")
for (i in 1:nrow(param.theta)){
#   unq = unique(param.theta[i,])
#   points(rep(i, length(unq)), unq, pch = 20, cex = 0.2)
    points(rep(i, n), param.theta[i,], pch = 20, cex = 0.2)
    }
matplot(param.theta, pch = 20, cex = 0.2)

### Box plots of the theta's
qlines = apply(param.theta, 2, quantile, c(0.025, 0.975))
mline = apply(param.theta, 2, mean)


#pdf("./figs/boxplots_1.pdf", width = 8, height = 8)
par(mfrow = c(1,1), mar = c(3.1, 2.1, 2.1, 1.1))
#boxplot(param.theta, pch = 20, cex = 0.1)
#points(y, col = 'red', pch = 20, cex = 0.5)
boxplot(param.theta[,ord], pch = 20, cex = 0.1)
points(y[ord], col = 'red', pch = 20, cex = 0.5)
lines(qlines[1,ord], col = 'green', lwd = 1.5)
lines(qlines[2,ord], col = 'green', lwd = 1.5)
lines(mline[ord], col = 'darkgreen', lwd = 1.5)
#dev.off()

#pdf("./figs/lines_1.pdf", width = 8, height = 8)
par(mfrow = c(1,1), mar = c(4.1, 2.1, 2.1, 1.1))
plot(1:n, type='n', ylim = c(range(y)),
    xlab = expression(theta ~ "index"), cex.lab = 2,
    main = expression("Posterior cluster locations for each" ~ theta[i]), cex.main = 2)
segments(x0 = 1:n, y0 = qlines[1,ord], y1 = qlines[2,ord], lwd = 1.0)
lines(mline[ord], col = 'darkgreen', lwd = 3.0)
points(y[ord], col = 'red', pch = 20, cex = 0.5)
legend("topleft", box.lty = 0, legend = c("Data (ordered)", "Mean", "95% credible interval"),
    lwd = c(NA, 3, 3), pch = c(20, NA, NA), cex = 1.5, col = c("red", "darkgreen", "black"))
#dev.off()




### Posterior for n*
hpd.alpha = hpd.uni(param.alpha)
clusters = apply(param.theta, 1, function(x) length(unique(x)))

### Get theta groups
library(doMC)
library(fields)
registerDoMC(4)
par.fun = function(j){
    out = double(n)
    temp = table(unlist(apply(param.theta, 1, function(x) which(x == x[j]))))
    out[as.numeric(names(temp))] = temp
    return (out)
    }
group = (foreach(j = 1:n, .combine = rbind) %dopar% par.fun(j)) / nmcmc


#par(mfrow = c(1,1), mar=c(5.1,4.1,4.1,2.1))
#image.plot(group)



#pdf("./figs/nstar_1.pdf", height = 8, width = 8)
par(mfrow = c(1,1), mar = c(4.1, 4.1, 2.1, 1.1))
layout(matrix(c(1,3,2,4),2,2))
plot(table(clusters) / nmcmc, main = expression(n^"*"), cex.main = 2, lwd = 3, ylab="Mass")
a = par("pin")
b = par("plt")
image.plot(group[ord, ord], main = "Cluster groupings heatmap")
par("pin" = a, "plt" = b)
plot(param.alpha, type='l', ylab = expression(alpha), main = "Posterior trace plot")
hpd.plot(density(param.alpha), hpd.alpha, main = expression(alpha), xlab = "", cex.main = 2)
#dev.off()


### Posterior for alpha, phi, mu, tau^2
hpd.beta = hpd.uni(param.beta)
hpd.zeta = hpd.uni(param.psi[,1])
hpd.mu = hpd.uni(param.psi[,2])

#pdf("./figs/posts_1.pdf", height = 8, width = 8)
par(mfrow = c(3,2), mar = c(3.1, 2.1, 2.1, 1.1), oma = c(0,0,2,0))
plot(param.beta, type='l', main = "Traceplot")
hpd.plot(density(param.beta), hpd.beta, main = expression(beta), xlab="", cex.main = 2.0)

plot(param.psi[,1], type='l', main = "Traceplot")
hpd.plot(density(param.psi[,1]), hpd.zeta, main = expression(zeta), xlab="", cex.main = 2)

plot(param.psi[,2], type='l', main = "Traceplot")
hpd.plot(density(param.psi[,2]), hpd.mu, main = expression(mu), xlab="", cex.main = 2)
title(main = "Posterior distributions (with 95% hpd set)", outer = TRUE, cex.main = 2)
#dev.off()
par(oma = c(0,0,0,0))

mean(param.beta)
quantile(param.beta, c(0, 0.025, 0.5, 0.975, 1))

apply(param.psi, 2, mean)
t(apply(param.psi, 2, quantile, c(0, 0.025, 0.5, 0.975, 1)))



#csort = group
##finalorder = double(n)
#for (i in 1:n){
#    cat("\r", i)
#    cord = order(csort[i, i:n], decreasing = TRUE)
#    if (i > 1)
#        cord = c(1:(i-1), cord + i - 1)
#    csort = csort[cord,cord]
##   finalorder[i] = which(cord == i)
#    image.plot(csort)
##   readline()
#    }
#
#image.plot(csort)
#image.plot(group[finalorder, finalorder])


par(mfrow = c(1, 1), mar = c(3.1, 2.1, 2.1, 1.1))
plot(0, type='n', xlim = range(param.theta), ylim = c(0, 3.0))
for (k in 1:n)
    lines(density(param.theta[,k]), col = k, lty = k)

#plot(0, type='n', xlim = range(pback), ylim = c(0, 3.0))
#for (k in 1:n)
#    lines(density(pback[,k]), col = k, lty = k)

### Posterior predictions
### list of theta_stars
n_j = apply(param.theta, 1, table)
theta_star = lapply(n_j, function(x) as.numeric(names(x)))

# A new cluster (to predict further with this, we would need to specify a particular x)
theta_0 = double(nmcmc)
for (i in 1:nmcmc){
    if (floor(i/window) == i/window)
        cat("\r", i, "/", nmcmc)
    prob = c(param.alpha[i] / (param.alpha[i] + n), n_j[[i]] / (param.alpha[i] + n))
    draw = sample(length(prob), 1, replace = FALSE, prob = prob)
    if (draw == 1){
        theta_0[i] = rgamma(1, param.psi[i,1], param.psi[i,1] / param.psi[i,2])
    } else {
        theta_0[i] = theta_star[[i]][draw - 1]
        }
    if (i == nmcmc)
        cat("\n")
    }
y_0 = matrix(0, nmcmc, n)
for (i in 1:n)
    y_0[,i] = rpois(nmcmc, theta_0*exp(param.beta*x[i]))

#hpd.theta_0 = hpd.mult(theta_0, density(theta_0))
#hpd.y_0 = hpd.mult(y_0, density(y_0))

##pdf("./figs/pred_1.pdf", height = 6, width = 12)
#par(mfrow = c(2, 1), mar = c(3.1, 2.1, 2.1, 1.1))
#hpd.plot(density(theta_0), hpd.theta_0, main = expression("Posterior predictive for" ~ theta[0] ~ "(with 95% hpd set)"), cex.main = 1, xlab = expression(theta[0]))
#
## # A new data point (just need theta_0 and phi)
## #par(mfrow = c(1, 1), mar = c(5.1, 4.1, 4.1, 2.1))
# hpd.plot(density(y_0), hpd.y_0, main = expression("Posterior predictive for" ~ y[0] ~ "(with 95% hpd set)"), cex.main = 1, xlab = expression(y[0]))
# lines(density(y), lwd =3)
# legend("topleft", box.lty = 0, legend = "Data", lwd = 3, cex = 1.5)
## #dev.off()


### predictions of observations
# Predictions at each x[i]
#y_0 = matrix(0, nmcmc, n)
#for (i in 1:n)
#    y_0[,i] = rpois(nmcmc, param.theta[,i]*exp(param.beta*x[i]))

q0 = apply(y_0, 2, quantile, c(0.025, 0.975))
m0 = apply(y_0, 2, mean)

### Predictions at each x_i
pdf("./figs/pred_3c.pdf", width = 12, height = 6)
par(mfrow = c(1,2), mar = c(4.1,4.1,2.1, 1.1))
plot(x, y, pch = 1, ylim = range(c(y, q0)), lwd = 1.5, main = "Predictions")
segments(x0 = x, y0 = q0[1,], y1 = q0[2,], col = 'forestgreen')
points(x, m0, col = 'darkgreen', pch = 20)
legend("topleft", box.lty = 0, legend = "Data", pch = 1, cex = 1.5)

### Observed vs. Fitted plot
plot(0, type='n', xlim = range(y), ylim = range(c(y, q0)),
    xlab = "Observed", ylab = "Fitted", main = "Observed vs. Fitted", cex.lab = 1.3)
segments(x0 = jity, y0 = q0[1,], y1 = q0[2,], col = 'forestgreen')
points(jity, m0, col = 'darkgreen', pch = 20)
abline(0, 1, lwd = 3, lty = 2)
dev.off()

pplc(y, y_0, 0)                     # 457.03
pplc(y, y_0, Inf)                   # 590.09
pplc(y, y_0, Inf) - pplc(y, y_0, 0) # 133.06

