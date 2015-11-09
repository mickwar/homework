y = as.numeric(unlist(read.table("./data.txt")))

#y = c(rnorm(20, -5, 1), rnorm(40, 15, 1))
n = length(y)

plot(density(y))

### Inverse gamma functions
# Mean = rate / (shape - 1)
# What I'm calling rate for the inverse gamma, Wikipedia calls scale
dinvgamma = function(x, shape, rate, log = FALSE){
    if (log)
        return (shape * log(rate) - lgamma(shape) - (shape + 1) * log(x) - rate / x)
    return (rate^shape / gamma(shape) * x^(-shape-1) * exp(-rate / x))
    }
rinvgamma = function(n, shape, rate)
    1/rgamma(n, shape = shape, rate = rate)

### Other functions
#calc.q0 = function(x, mu, tau2, phi){
#    sig = phi*tau2 / (phi + tau2)
##   sqrt(sig/(2*pi*phi*tau2)) * exp(mu^2*phi*(1-phi) - 2*mu*phi*x*tau2 + x^2*tau2*(1-tau2))
#    exp(1/2 * (log(sig) - log(2*pi*phi*tau2)) +
#        mu^2*phi*(1-phi) - 2*mu*phi*x*tau2 + x^2*tau2*(1-tau2))
#    }

#x = y[1]
#mu = param.psi[1,1]
#tau2 = param.psi[1,2]
#phi = param.phi[1]



### MCMC
nburn = 10000
nmcmc = 5000

# Parameter objects
param.theta = matrix(0, nburn + nmcmc, n)
param.alpha = double(nburn + nmcmc)
param.psi = matrix(0, nburn + nmcmc, 2) # [,1] is mu, [,2] is tau^2
param.phi = double(nburn + nmcmc)

# Priors
prior.phi.a   = 3  # Inverse Gamma (shape)
prior.phi.b   = 20 # Inverse Gamma (rate)
prior.alpha.a = 1  # Gamma (shape)           alpha controls n.star (discreteness of G)
prior.alpha.b = 1  # Gamma (rate)
prior.mu.mean = 0  # Normal (mean)
prior.mu.var  = 3  # Normal (variance)
prior.tau2.a  = 3  # Inverse Gamma (shape)
prior.tau2.b  = 20 # Inverse Gamma (rate)

# Initial values
param.theta[1,] = y
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

#       # mean of h
#       mu.star = (c.mu * c.phi + y[i] * c.tau2) / (c.phi + c.tau2)
#       # variance of h
#       sig.star = c.phi * c.tau2 / (c.phi + c.tau2)

        # q0 here is a normal density (after annoying integration)
        q0 = dnorm(y[i], c.mu, sqrt(c.phi + c.tau2))

        # as is qj, by construction
        qj = dnorm(y[i], theta.star.minus, sqrt(c.phi))


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
    theta.star = as.numeric(names(table(c.theta)))
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

    # Put the c. objects into the regular ones
    param.theta[iter,] = c.theta
    param.alpha[iter] = c.alpha
    param.psi[iter,1] = c.mu
    param.psi[iter,2] = c.tau2
    param.phi[iter] = c.phi

    if (iter == (nburn + nmcmc))
        cat("\n")
    }

par(mfrow = c(4, 2), mar = c(3.1, 2.1, 2.1, 1.1))
for (k in 1:n){
    plot(tail(param.theta[,k], nmcmc), type = 'l')
    plot(density(tail(param.theta[,k], nmcmc)))
    if (k %% 4 == 0)
        readline()
    }

par(mfrow = c(1,1), mar = c(5.1, 4.1, 4.1, 2.1))
clusters = apply(tail(param.theta, nmcmc), 1, function(x) length(unique(x)))
plot(table(clusters) / nmcmc)
table(clusters) / nmcmc

par(mfrow = c(1,2), mar = c(3.1, 2.1, 2.1, 1.1))
plot(tail(param.alpha, nmcmc), type='l')
plot(density(tail(param.alpha, nmcmc)))

plot(tail(param.phi, nmcmc), type='l')
plot(density(tail(param.phi, nmcmc)))

plot(tail(param.psi[,1], nmcmc), type='l')
plot(density(tail(param.psi[,1], nmcmc)))

plot(tail(param.psi[,2], nmcmc), type='l')
plot(density(tail(param.psi[,2], nmcmc)))

