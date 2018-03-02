library(MASS)
library(mwBASE)

DATA_DIR = "~/files/tmp/268_data/"

### For generating X
make.sigma = function(k, rho){
    out = matrix(rho, k, k)
    diag(out) = 1
    return (out)
    }
### Normal Inverse Gamma density
dnig = function(beta, sig2, mu, lambda, a, b, log = TRUE){
    k = length(mu)
    out = -(a + k/2 + 1)*log(sig2) - (1/sig2) *
        (b + 1/2*t(beta - mu) %*% lambda %*% (beta - mu))
    if (log)
        return (out)
    return (exp(out))
    }
### Random Inverse Gaussian
rigauss = function(n, mu, lambda){
    nu = rnorm(n)
    y = nu^2
    x = mu + y*mu^2/(2*lambda) - mu/(2*lambda) *
        sqrt(4*mu*lambda*y + mu^2*y^2)
    # account of when mu is large (becomes a levy)
    x = ifelse(x < 0, 0, x)
    ind = which(x == 0)
    k = length(ind)

    z = rnorm(n-k)
    if (k > 0){
        x[ind] = lambda[ind] / qnorm(1 - runif(k)/2)^2
        out = x
        out[-ind] = ifelse(z <= mu[-ind] / (mu[-ind] + x[-ind]),
            x[-ind], mu[-ind]^2 / x[-ind])
    } else {
        out = ifelse(z <= mu / (mu + x), x, mu^2 / x)
        }
    return (out)
    }



### Generate and write the data to m files
n = 10^6    # no. data points   (10^8)
m = 100     # no. of subsamples, subsets    (1000)
            # Make sure n/m is an integer and n/m > k

k = 40     # no. of predictors  (100)
rho = 0.99  # corr. between predictors (0.99)
sigma = 10  # measurement error

beta = c(seq(1, 0.1, by = -0.1), rep(0, k-10))
sig.mat = make.sigma(k, rho)

set.seed(1)
for (i in 1:m){
    X = mvrnorm(n/m, rep(0, k), sig.mat)
    y = rnorm(n/m, mean = X %*% beta, sd = sigma)
    fname = paste0(c(DATA_DIR, "dat_", rep("0", nchar(m) - nchar(i)), i, ".txt"), collapse = "")
    write.table(round(data.frame(y, X), 3), quote = FALSE, row.names = FALSE, file = fname)
    rm(X, y)
    gc()
    }


### Read in data and calculate posterior distributions
### for each data subset
m = length(list.files(DATA_DIR))
fname = paste0(c(DATA_DIR, "dat_", rep("0", nchar(m) - nchar(1)), 1, ".txt"), collapse = "")
z = read.table(fname, header = TRUE)
n = NROW(z)*m
k = NCOL(z) - 1
rm(z)

mu.s = matrix(0, k, m)
lambda.s = rep(list(matrix(0, k, k)), m)
a.s = double(m)
b.s = double(m)
for (i in 1:m){
    fname = paste0(c(DATA_DIR, "dat_", rep("0", nchar(m) - nchar(i)), i, ".txt"), collapse = "")
    z = read.table(fname, header = TRUE)
    y = as.vector(z[,1])
    X = as.matrix(z[,-1])

    xty = t(X) %*% y
    yty = sum(y^2)

    lambda.s[[i]] = t(X) %*% X
    lam.inv = solve(t(X) %*% X)
    mu.s[,i] = lam.inv %*% xty
    a.s[i] = (n/m - k) / 2
    b.s[i] = 1/2 * (yty - t(xty) %*% lam.inv %*% xty)

    rm(X, y, z)
    gc()
    }


### Calculate full posterior (sum of the subsets)
# (Algorithm 1)
mu = mu.s[,1]
lambda = lambda.s[[1]]
a = a.s[1]
b = b.s[1]

for (i in 2:m){
    tmp.mu = mu
    mu = solve(lambda + lambda.s[[i]]) %*% (lambda %*% mu + lambda.s[[i]] %*% mu.s[,i])
    a = a + a.s[i] + k/2
    b = b + b.s[i] + 1/2*t(tmp.mu - mu) %*% lambda %*% (tmp.mu - mu) +
        1/2*t(mu.s[,i] - mu) %*% lambda.s[[i]] %*% (mu.s[,i] - mu)
    lambda = lambda + lambda.s[[i]]
    }
# Non-informative prior




### LASSO
# (Algorithm 2, start with Algorithm 1)

# Penalty terms
pen.vec = rep(10, k)

# Gibbs set up
nburn = 500
nmcmc = 2000
B = nburn + nmcmc
draws.beta = matrix(0, B, k)
draws.sig2 = double(B)
draws.psi = matrix(0, B, k)

draws.sig2[1] = 1
draws.psi[1,] = rep(1, k)
for (i in 2:B){
    lambda.psi = diag(1/draws.psi[i-1,])

    # Get NIG parameters
    mu.bar = solve(lambda + lambda.psi) %*% (lambda %*% mu + lambda.psi %*% rep(0, k))
    a.bar = a + 0
    b.bar = b + 0 + 1/2*t(mu - mu.bar) %*% lambda %*% (mu - mu.bar) +
        1/2*t(rep(0, k) - mu.bar) %*% lambda.psi %*% (rep(0, k) - mu.bar)
    lambda.bar = lambda + lambda.psi

    # Update sigma^2
    draws.sig2[i] = 1/rgamma(1, a.bar, b.bar)

    # Update beta
    draws.beta[i,] = mvrnorm(1, mu.bar, draws.sig2[i] * solve(lambda.bar))

    # Update psi
    draws.psi[i,] = 1/rigauss(k, abs(pen.vec * draws.sig2[i] / draws.beta[i,]), pen.vec^2)
    }


draws.beta = tail(draws.beta, nmcmc)
draws.sig2 = tail(draws.sig2, nmcmc)
draws.psi = tail(draws.psi, nmcmc)


mean.beta = colMeans(draws.beta)
med.beta = apply(draws.beta, 2, median)
qq.beta = apply(draws.beta, 2, quantile, c(0.025, 0.975))

plot(mean.beta, pch = 16, col = 'darkgreen', ylim = range(qq.beta))
points(1:k, med.beta, pch = 16, col = 'dodgerblue')
segments(1:k, qq.beta[1,], 1:k, qq.beta[2,], col = 'forestgreen')

points(1:k, beta, pch = 15)


### SSVS
# (Algorithm 3)
# Notes: It seems the spike part (d) needs to have a really small variance
# (0.001^2) for this data set. This is probably because the observation variance
# is included when calculated the bernoulli probability for gamma.
# In HW1, that variance was not included. Was this a mistake from Qian
# or Raj? I want to say there is a mistake in Qian, since including
# sigma^2 in calculating the probabilities must be dealt with when
# selecting d and D. But maybe not, the variance of beta.hat in linear regression
# does include sigma^2. So...
# How to choose reasonable d and D?
#
# The above comment was written after running the model on data with no 
# correlation between covariates. When there was strong correlation (0.99),
# then every variable was selected by SSVS, the spread of the each beta
# was much larger, those betas closer to 0 had credible intervals crossing 0.
# Basically, this method, as it is, is not sufficient when dealing with
# correlated predictors.

# Variance for the spike (d) and slab (D) components
d = 0.01^2
D = 100

# Beta prior for w (mixture weight)
a.w = rep(1, k)
b.w = rep(1, k)


# Gibbs set up
nburn = 2000
nmcmc = 5000
B = nburn + nmcmc
draws.beta = matrix(0, B, k)
draws.sig2 = double(B)
draws.gamma = matrix(0, B, k)   # Indicators for selected variables
draws.w = matrix(0, B, k)       # Probabilities for indicators

draws.sig2[1] = 1
draws.gamma[1,] = rep(0, k)     # Start with no selected variables
draws.w[1,] = rep(0.5, k)

for (i in 2:B){
    lambda.gamma = diag(1/ifelse(draws.gamma[i-1,] == 1, D, d))

    # Get NIG parameters
    mu.bar = solve(lambda + lambda.gamma) %*% (lambda %*% mu + lambda.gamma %*% rep(0, k))
    a.bar = a + 0
    b.bar = b + 0 + 1/2*t(mu - mu.bar) %*% lambda %*% (mu - mu.bar) +
        1/2*t(rep(0, k) - mu.bar) %*% lambda.gamma %*% (rep(0, k) - mu.bar)
    lambda.bar = lambda + lambda.gamma

    # Update sigma^2
    draws.sig2[i] = 1/rgamma(1, a.bar, b.bar)

    # Update beta
    draws.beta[i,] = mvrnorm(1, mu.bar, draws.sig2[i] * solve(lambda.bar))

    # Update gamma
#   h1 = draws.w[i-1,] * dnorm(draws.beta[i,], 0, sqrt(draws.sig2[i] * D))
#   h2 = (1 - draws.w[i-1,]) * dnorm(draws.beta[i,], 0, sqrt(draws.sig2[i] * d))
    h1 = draws.w[i-1,] * dnorm(draws.beta[i,], 0, sqrt(D))
    h2 = (1 - draws.w[i-1,]) * dnorm(draws.beta[i,], 0, sqrt(d))
    draws.gamma[i,] = rbinom(k, 1, h1 / (h1 + h2))

    # Update w
    draws.w[i,] = rbeta(k, a.w + draws.gamma[i,], b.w + 1 - draws.gamma[i,]) 
    }

draws.beta = tail(draws.beta, nmcmc)
draws.sig2 = tail(draws.sig2, nmcmc)
draws.psi = tail(draws.psi, nmcmc)

colMeans(draws.w)


cols = rep("white", k)
cols[which(colMeans(draws.w) > 0.5)] = "dodgerblue"
boxplot(draws.beta, col = cols)
points(1:k, beta, col = 'firebrick', pch = 15)

plot_hpd(draws.sig2)


### MC^3
# (Algorithm 4)

# Variance for the spike (d) and slab (D) components
d = 0.01^2
D = 100

# Beta prior for w (mixture weight)
a.w = rep(1, k)
b.w = rep(1, k)


# Gibbs set up
nburn = 200
nmcmc = 500
B = nburn + nmcmc
draws.beta = matrix(0, B, k)
draws.sig2 = double(B)
draws.gamma = matrix(0, B, k)   # Indicators for selected variables
# draws.w = matrix(0, B, k)       # Probabilities for indicators

draws.sig2[1] = 1
draws.gamma[1,] = rep(0, k)     # Start with no selected variables
# draws.w[1,] = rep(0.5, k)

for (i in 2:B){
    # 4.1, propose gamma
    gamma.prop = draws.gamma[i-1,]
    tmp.ind = sample(k, 1)
    gamma.prop[tmp.ind] = 1 - gamma.prop[tmp.ind]
    vg = (gamma.prop == 1)
    vo = (gamma.prop == 0)
    kg = sum(vg)
    ko = sum(vo)

    # 4.2, reduce
    beta.g = draws.beta[i-1, gamma.prop == 1]
    beta.o = draws.beta[i-1, gamma.prop == 0]
    sig2 = draws.sig2[i-1]

    mu.g = mu[gamma.prop == 1]
    mu.o = mu[gamma.prop == 0]

    lam.gg = matrix(lambda[vg, vg], kg, kg)
    lam.go = matrix(lambda[vg, vo], kg, ko)
    lam.og = matrix(lambda[vo, vg], ko, kg)
    lam.oo = matrix(lambda[vo, vg], ko, ko)

    # gamma subscripts
    mu.bar = mu.g + solve(lam.gg) %*% lam.go %*% mu.o
    lam.bar = lam.gg
    a.bar = a + ko / 2
    b.bar = b + 1/2 * t(mu.o) %*% (lam.oo - lam.og %*% solve(lam.gg) %*% lam.go) %*% mu.o

    dnig(beta.g, sig2, , lam.bar, a.bar, b.bar)
    dnig(beta.g, sig2, mu.bar, lam.bar, a.bar, b.bar)

    dnig(rep(0, kg), sig2, rep(0, kg), matrix(0, kg, kg), -kg/2, 0) -
        dnig(rep(0, kg), sig2, mu.bar, lam.bar, a.bar, b.bar)


    # Update sigma^2
    draws.sig2[i] = 1/rgamma(1, a.bar, b.bar)

    # Update beta
    draws.beta[i,] = mvrnorm(1, mu.bar, draws.sig2[i] * solve(lambda.bar))

    # Update gamma
#   h1 = draws.w[i-1,] * dnorm(draws.beta[i,], 0, sqrt(draws.sig2[i] * D))
#   h2 = (1 - draws.w[i-1,]) * dnorm(draws.beta[i,], 0, sqrt(draws.sig2[i] * d))
    h1 = draws.w[i-1,] * dnorm(draws.beta[i,], 0, sqrt(D))
    h2 = (1 - draws.w[i-1,]) * dnorm(draws.beta[i,], 0, sqrt(d))
    draws.gamma[i,] = rbinom(k, 1, h1 / (h1 + h2))

    # Update w
    draws.w[i,] = rbeta(k, a.w + draws.gamma[i,], b.w + 1 - draws.gamma[i,]) 
    }
