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

y = as.numeric(unlist(read.table("./data.txt")))

n = length(y)
ord = order(y)

par(mfrow = c(1,1), mar = c(5.1, 4.1, 4.1, 2.1))
plot(density(y))
curve(0.2*dnorm(x, -5, 1) + 0.5*dnorm(x, 0, 1) + 0.3*dnorm(x, 3.5, 1), add = TRUE, col = 'red')
legend("topleft", legend = c("data", "truth"), col = c("black", "red"), lwd = 1, box.lwd = 0, cex = 1.5)

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


range(rnorm(10000, rnorm(10000, 0, sqrt(3)), rinvgamma(10000, 3, rinvgamma(10000, 5, 10)) / rgamma(10000, 1, 1)))


### MCMC
nburn = 2000
nmcmc = 10000

# Parameter objects
param.theta = matrix(0, nburn + nmcmc, n)
param.alpha = double(nburn + nmcmc)
param.psi = matrix(0, nburn + nmcmc, 2) # [,1] is mu, [,2] is tau^2
param.phi = double(nburn + nmcmc)

# Priors
prior.phi.a   = 3   # Inverse Gamma (shape)
prior.phi.b   = 10  # Inverse Gamma (rate)
prior.alpha.a = 1   # Gamma (shape)           alpha controls n.star (discreteness of G)
prior.alpha.b = 1   # Gamma (rate)
prior.mu.mean = 0   # Normal (mean)
prior.mu.var  = 3   # Normal (variance)
prior.tau2.a  = 3   # Inverse Gamma (shape)
prior.tau2.b  = 10  # Inverse Gamma (rate)


# Initial values
#param.theta[1,] = y
param.phi[1] = rinvgamma(1, prior.phi.a, prior.phi.b)
param.alpha[1] = rgamma(1, prior.alpha.a, prior.alpha.b)
param.psi[1,] = c(rnorm(1, prior.mu.mean, sqrt(prior.mu.var)),
    rinvgamma(1, prior.tau2.a, prior.tau2.b))

# Iterations
for (iter in 2:(nburn + nmcmc)){
    cat("\r", iter, "/", nburn + nmcmc)

    # Current parameters (used for convenience)
    c.theta = param.theta[iter-1,]
    c.alpha = param.alpha[iter-1]
    c.mu = param.psi[iter-1,1]
    c.tau2 = param.psi[iter-1,2]
    c.phi = param.phi[iter-1]

    # Update thetas
    for (i in 1:n){
        n.j.minus = table(c.theta[-i])
        theta.star.minus = as.numeric(names(n.j.minus))
        n.star.minus = length(theta.star.minus)

        # q0 here is a normal density (after annoying integration)
        q0 = dnorm(y[i], c.mu, sqrt(c.phi + c.tau2))

        # as is qj, by construction
        qj = dnorm(y[i], theta.star.minus, sqrt(c.phi))

        # Probabilities of determining which to draw
        # Calculate A
        A = c.alpha * q0 / (c.alpha * q0 + sum(n.j.minus * qj))
        
        # Calculate B's
        Bj = n.j.minus * qj / (c.alpha * q0 + sum(n.j.minus * qj))

        # Make the update
        draw = sample(n.star.minus + 1, 1, prob = c(A, Bj))
        if (draw == 1){ # Make a draw from h
            c.theta[i] = rnorm(1, (c.mu*c.phi + y[i]+c.tau2)/(c.phi + c.tau2),
                sqrt(c.phi*c.tau2/(c.phi + c.tau2)))
        } else { # Make a draw from the existing groups
            c.theta[i] = theta.star.minus[draw-1]
            }
        }
    n.j = table(c.theta)
    theta.star = as.numeric(names(n.j))
    n.star = length(theta.star)

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

    # Update psi = (mu, tau^2)
    S = sum(theta.star)
    c.mu = rnorm(1, (prior.mu.mean*c.tau2 + prior.mu.var * S) /
        (c.tau2 + prior.mu.var * n.star),
        sqrt(c.tau2 * prior.mu.var / (c.tau2 + prior.mu.var * n.star)))
    c.tau2 = rinvgamma(1, prior.tau2.a + n.star/2,
        prior.tau2.b + 1/2 * sum((theta.star - c.mu)^2))

    # Update phi
    c.phi = rinvgamma(1, prior.phi.a + n/2, prior.phi.b + 1/2 * sum((y - c.theta)^2))

    # Extra improvement step (only updating the clusters, not individuals)
    for (j in 1:n.star){
        temp.w = which(c.theta == theta.star[j])
        c.theta[temp.w] = rnorm(1,
            (c.phi*c.mu + c.tau2*sum(y[temp.w])) / (c.phi + n.j[j]*c.tau2),
            sqrt(c.phi*c.tau2 / (c.phi + n.j[j]*c.tau2)))
        }

    # Put the c.* objects into the regular ones
    param.theta[iter,] = c.theta
    param.alpha[iter] = c.alpha
    param.psi[iter,1] = c.mu
    param.psi[iter,2] = c.tau2
    param.phi[iter] = c.phi

    if (iter == (nburn + nmcmc))
        cat("\n")
    }

### Burn in
param.theta = tail(param.theta, nmcmc)
param.alpha = tail(param.alpha, nmcmc)
param.phi = tail(param.phi, nmcmc)
param.psi = tail(param.psi, nmcmc)


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
#pdf("./figs/improve_comp.pdf", height = 9, width = 9)
#par(mfrow = c(1,1), mar = c(3.1, 4.1, 2.1, 1.1))
#plot(0, type='n', xlim = c(1, nrow(pback)), ylim = range(pback),
#    xlab = "MCMC Iteration", ylab = "Unique cluster locations", main = "Without improving Gibbs step")
#for (i in 1:nrow(pback)){
#    unq = unique(pback[i,])
#    points(rep(i, length(unq)), unq, pch = 20, cex = 0.2)
#    }
#plot(0, type='n', xlim = c(1, nrow(param.theta)), ylim = range(param.theta),
#    xlab = "MCMC Iteration", ylab = "Unique cluster locations", main = "With improving Gibbs step")
#for (i in 1:nrow(param.theta)){
#    unq = unique(param.theta[i,])
#    points(rep(i, length(unq)), unq, pch = 20, cex = 0.2)
#    }
#dev.off()

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
plot(1:250, type='n', xlim = c(1, 250), ylim = c(range(y)),
    xlab = expression(theta ~ "index"), cex.lab = 2,
    main = expression("Posterior cluster locations for each" ~ theta[i]), cex.main = 2)
segments(x0 = 1:250, y0 = qlines[1,ord], y1 = qlines[2,ord], lwd = 1.0)
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


#table(clusters) / nmcmc

#plot(clusters, type='l')
#abline(h=3, lty=2)

#cback = apply(pback, 1, function(x) length(unique(x)))
#lines(as.numeric(names(table(cback)))+0.2, table(cback) / nmcmc, main = expression(n^"*"), cex.main = 2, col = 'red', type='h', lwd = 3)
#table(cback) / nmcmc

#plot(density(apply(param.theta[clusters == 3,], 1, unique)))
#x  =apply(param.theta[clusters == 3,], 1, unique)

### Posterior for alpha, phi, mu, tau^2
hpd.phi = hpd.uni(param.phi)
hpd.mu = hpd.uni(param.psi[,1])
hpd.tau2 = hpd.uni(param.psi[,2])

#pdf("./figs/posts_1.pdf", height = 8, width = 8)
par(mfrow = c(3,2), mar = c(3.1, 2.1, 2.1, 1.1), oma = c(0,0,2,0))
plot(param.phi, type='l', main = "Traceplot")
hpd.plot(density(param.phi), hpd.phi, main = expression(phi), xlab="", cex.main = 2.0)

plot(param.psi[,1], type='l', main = "Traceplot")
hpd.plot(density(param.psi[,1]), hpd.mu, main = expression(mu), xlab="", cex.main = 2)

plot(param.psi[,2], type='l', main = "Traceplot")
hpd.plot(density(param.psi[,2]), hpd.tau2, main = expression(tau^2), xlab="", cex.main = 2)
title(main = "Posterior distributions (with 95% hpd set)", outer = TRUE, cex.main = 2)
#dev.off()



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

# A new cluster
theta_0 = double(nmcmc)
for (i in 1:nmcmc){
    cat("\r", i, "/", nmcmc)
    prob = c(param.alpha[i] / (param.alpha[i] + n), n_j[[i]] / (param.alpha[i] + n))
    draw = sample(length(prob), 1, replace = FALSE, prob = prob)
    if (draw == 1){
        theta_0[i] = rnorm(1, param.psi[i, 1], sqrt(param.psi[i, 2]))
    } else {
        theta_0[i] = theta_star[[i]][draw - 1]
        }
    if (i == nmcmc)
        cat("\n")
    }
y_0 = rnorm(nmcmc, theta_0, sqrt(param.phi))

hpd.theta_0 = hpd.mult(theta_0, density(theta_0))
hpd.y_0 = hpd.mult(y_0, density(y_0))

#pdf("./figs/pred_1.pdf", height = 6, width = 12)
par(mfrow = c(1, 2), mar = c(3.1, 2.1, 2.1, 1.1))
hpd.plot(density(theta_0), hpd.theta_0, main = expression("Posterior predictive for" ~ theta[0] ~ "(with 95% hpd set)"), cex.main = 1, xlab = expression(theta[0]))

# A new data point (just need theta_0 and phi)
#par(mfrow = c(1, 1), mar = c(5.1, 4.1, 4.1, 2.1))
hpd.plot(density(y_0), hpd.y_0, main = expression("Posterior predictive for" ~ y[0] ~ "(with 95% hpd set)"), cex.main = 1, xlab = expression(y[0]))
lines(density(y), lwd =3)
legend("topleft", box.lty = 0, legend = "Data", lwd = 3, cex = 1.5)
#dev.off()

#y2 = matrix(0, nmcmc, n)
#for (i in 1:n)
#    y2[,i] = rnorm(nmcmc, param.theta[,i], sqrt(param.phi))

### Replicate for each observation (don't really need to, but it makes pplc work)
y1 = matrix(rep(y_0, n), ncol = n)

pplc(y, y1, 0)                      # 2353.17
pplc(y, y1, Inf)                    # 4707.39
pplc(y, y1, Inf) - pplc(y, y1, 0)   # 2354.22

### Prior predictive
library(MCMCpack)
# Cluster

prior.theta_0 = double(20000)
prior.y_0 = double(length(prior.theta_0))
tvec = seq(-10, 10, length = 1500)
k = length(tvec)
for (i in 1:length(prior.y_0)){
    cat("\r", i, "/", length(prior.y_0))
    temp.phi = rinvgamma(1, prior.phi.a, prior.phi.b)
    temp.alpha = rgamma(1, prior.alpha.a, prior.alpha.b)
    temp.mu = rnorm(1, prior.mu.mean, sqrt(prior.mu.var))
    temp.tau = sqrt(rinvgamma(1, prior.tau2.a, prior.tau2.b))
    dG0 = pnorm(tvec, temp.mu, temp.tau) - pnorm(c(-Inf, tvec[-k]), temp.mu, temp.tau)
    flag = TRUE
    while (flag){
        flag = FALSE
        dG = rdirichlet(1, temp.alpha * dG0)
        if (any(is.na(dG)))
            flag = TRUE
        }
    keep = which(dG != 0)
    if (length(keep) == 1){
        prior.theta_0[i] = tvec[keep]
    } else {
        prior.theta_0[i] = sample(tvec[keep], 1, prob = dG[keep])
        }
    prior.y_0[i] = rnorm(1, prior.theta_0[i], sqrt(temp.phi))
    if (i == length(prior.y_0))
        cat("\n")
    }

pdf("./figs/prior_1.pdf", height = 6, width = 12)
par(mfrow = c(1,2), mar = c(3.1, 2.1, 2.1, 1.1))
plot(density(prior.theta_0[1:i]), main = expression("Prior predictive for" ~ theta[0]), lwd = 2,
    col = 'blue')
#curve(dnorm(x, mean(prior.theta_0[1:i]), sd(prior.theta_0[1:i])), add = TRUE, col = 'red', lwd = 2)
legend("topleft", legend = "KDE of draws", box.lwd = 0, box.lty = 0, col = 'blue', lwd = 2, cex = 1.3)
#legend("topleft", legend = c("KDE of draws", "Approximate normal"),
#    box.lwd = 0, box.lty = 0, col = c('blue', 'red'), lwd = 2, cex = 1.3)

plot(density(prior.y_0[1:i]), main = expression("Prior predictive for" ~ y[0]), lwd = 2,
    col = 'blue')
#curve(dnorm(x, mean(prior.y_0[1:i]), sd(prior.y_0[1:i])), add = TRUE, col = 'red', lwd = 2)
lines(density(y), lwd =3)
legend("topleft", legend = c("KDE of draws", "Data"),
    box.lwd = 0, box.lty = 0, col = c('blue', 'black'), lwd = 2, cex = 1.3)
#legend("topleft", legend = c("KDE of draws", "Approximate normal", "Data"),
#    box.lwd = 0, box.lty = 0, col = c('blue', 'red', 'black'), lwd = 2, cex = 1.3)
dev.off()

