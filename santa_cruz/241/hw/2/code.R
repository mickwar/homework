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


### MCMC
nburn = 000
nmcmc = 1000

# Parameter objects
param.theta = matrix(0, nburn + nmcmc, n)
param.alpha = double(nburn + nmcmc)
param.psi = matrix(0, nburn + nmcmc, 2) # [,1] is mu, [,2] is tau^2
param.phi = double(nburn + nmcmc)

# Joe's Priors
# prior.phi.a   = 1.5  # Inverse Gamma (shape)
# prior.phi.b   = 1 # Inverse Gamma (rate)
# prior.alpha.a = 2  # Gamma (shape)           alpha controls n.star (discreteness of G)
# prior.alpha.b = 3  # Gamma (rate)
# prior.mu.mean = median(y)  # Normal (mean)
# prior.mu.var  = 1  # Normal (variance)
# prior.tau2.a  = 10 # Inverse Gamma (shape)
# prior.tau2.b  = 0.5 # Inverse Gamma (rate)

# Arthur's Priors
# prior.phi.a   = 3  # Inverse Gamma (shape)
# prior.phi.b   = 3 # Inverse Gamma (rate)
# prior.alpha.a = 1  # Gamma (shape)           alpha controls n.star (discreteness of G)
# prior.alpha.b = 1  # Gamma (rate)
# prior.mu.mean = 0  # Normal (mean)
# prior.mu.var  = 3  # Normal (variance)
# prior.tau2.a  = 3  # Inverse Gamma (shape)
# prior.tau2.b  = 3 # Inverse Gamma (rate)

# Priors
prior.phi.a   = 3  # Inverse Gamma (shape)
prior.phi.b   = 10 # Inverse Gamma (rate)
prior.alpha.a = 1  # Gamma (shape)           alpha controls n.star (discreteness of G)
prior.alpha.b = 1  # Gamma (rate)
prior.mu.mean = 0  # Normal (mean)
prior.mu.var  = 3  # Normal (variance)
prior.tau2.a  = 3  # Inverse Gamma (shape)
prior.tau2.b  = 10 # Inverse Gamma (rate)


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
par(mfrow = c(4, 2), mar = c(3.1, 2.1, 2.1, 1.1))
for (k in 1:n){
    plot(param.theta[,k], type = 'l')
    abline(h = y[k], col = 'red', lty = 2)
    plot(density(param.theta[,k]), main = k)
    abline(v = y[k], col = 'red', lty = 2)
    if (k %% 4 == 0)
        readline()
    }

### "Trace plot" for clusters
pdf("./figs/improve_comp.pdf", height = 9, width = 9)
par(mfrow = c(2,1), mar = c(3.1, 4.1, 2.1, 1.1))
plot(0, type='n', xlim = c(1, nrow(pback)), ylim = range(pback),
    xlab = "MCMC Iteration", ylab = "Unique cluster locations", main = "Without improving Gibbs step")
for (i in 1:nrow(pback)){
    unq = unique(pback[i,])
    points(rep(i, length(unq)), unq, pch = 20, cex = 0.2)
    }
plot(0, type='n', xlim = c(1, nrow(param.theta)), ylim = range(param.theta),
    xlab = "MCMC Iteration", ylab = "Unique cluster locations", main = "With improving Gibbs step")
for (i in 1:nrow(param.theta)){
    unq = unique(param.theta[i,])
    points(rep(i, length(unq)), unq, pch = 20, cex = 0.2)
    }
dev.off()

### Box plots of the theta's
qlines = apply(param.theta, 2, quantile, c(0.025, 0.975))
mline = apply(param.theta, 2, mean)

par(mfrow = c(1,1), mar = c(3.1, 2.1, 2.1, 1.1))
#boxplot(param.theta, pch = 20, cex = 0.1)
#points(y, col = 'red', pch = 20, cex = 0.5)
boxplot(param.theta[,ord], pch = 20, cex = 0.1)
points(y[ord], col = 'red', pch = 20, cex = 0.5)
lines(qlines[1,ord], col = 'green', lwd = 1.5)
lines(qlines[2,ord], col = 'green', lwd = 1.5)
lines(mline[ord], col = 'darkgreen', lwd = 1.5)


### Posterior for n*
par(mfrow = c(1,1), mar = c(5.1, 4.1, 4.1, 2.1))
clusters = apply(param.theta, 1, function(x) length(unique(x)))
plot(table(clusters) / nmcmc, main = expression(n^"*"), cex.main = 2, lwd = 3)
#plot(table(clusters) / nmcmc, main = expression(n^"*"), cex.main = 2, xlim = c(1, 35), ylim = c(0,0.15), lwd = 3)
table(clusters) / nmcmc

plot(clusters, type='l')
abline(h=3, lty=2)

cback = apply(pback, 1, function(x) length(unique(x)))
lines(as.numeric(names(table(cback)))+0.2, table(cback) / nmcmc, main = expression(n^"*"), cex.main = 2, col = 'red', type='h', lwd = 3)
table(cback) / nmcmc

#plot(density(apply(param.theta[clusters == 3,], 1, unique)))
#x  =apply(param.theta[clusters == 3,], 1, unique)

### Posterior for alpha, phi, mu, tau^2
par(mfrow = c(4,2), mar = c(3.1, 2.1, 2.1, 1.1))
plot(param.alpha, type='l')
plot(density(param.alpha), main = expression(alpha))

plot(param.phi, type='l')
plot(density(param.phi), main = expression(phi))

plot(param.psi[,1], type='l')
plot(density(param.psi[,1]), main = expression(mu))

plot(param.psi[,2], type='l')
plot(density(param.psi[,2]), main = expression(tau^2))



### Get theta groups
library(doMC)
registerDoMC(4)
par.fun = function(j){
    out = double(n)
    temp = table(unlist(apply(param.theta, 1, function(x) which(x == x[j]))))
    out[as.numeric(names(temp))] = temp
    return (out)
    }
group = (foreach(j = 1:n, .combine = rbind) %dopar% par.fun(j)) / nmcmc


library(fields)
par(mfrow = c(1,1), mar=c(5.1,4.1,4.1,2.1))
image.plot(group)
image.plot(group[ord, ord])

csort = group
#finalorder = double(n)
for (i in 1:n){
    cat("\r", i)
    cord = order(csort[i, i:n], decreasing = TRUE)
    if (i > 1)
        cord = c(1:(i-1), cord + i - 1)
    csort = csort[cord,cord]
#   finalorder[i] = which(cord == i)
    image.plot(csort)
#   readline()
    }

image.plot(csort)
image.plot(group[finalorder, finalorder])


par(mfrow = c(1, 1), mar = c(3.1, 2.1, 2.1, 1.1))
plot(0, type='n', xlim = range(param.theta), ylim = c(0, 3.0))
for (k in 1:n)
    lines(density(param.theta[,k]), col = k, lty = k)

plot(0, type='n', xlim = range(pback), ylim = c(0, 3.0))
for (k in 1:n)
    lines(density(pback[,k]), col = k, lty = k)
