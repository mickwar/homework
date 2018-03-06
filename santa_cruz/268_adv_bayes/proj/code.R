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
dnig = function(beta, sig2, mu, lambdainv, a, b, log = TRUE){
    k = length(mu)
    out = -(a + k/2 + 1)*log(sig2) - (1/sig2) *
        (b + 1/2*t(beta - mu) %*% lambdainv %*% (beta - mu))
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
### NIG summation
nig_sum = function(mu1, lambda1, a1, b1, mu2, lambda2, a2, b2){
    k = length(mu1)
    mu = solve(lambda1 + lambda2) %*% (lambda1 %*% mu1 + lambda2 %*% mu2)
    a = a1 + a2 + k/2
    b = b1 + b2 + 1/2*t(mu1 - mu) %*% lambda1 %*% (mu1 - mu) +
        1/2*t(mu2 - mu) %*% lambda2 %*% (mu2 - mu)
    lambda = lambda1 + lambda2
    return (list("mu"=mu, "lambda"=lambda, "a"=a, "b"=b))
    }



### Generate and write the data to m files
n = 10^6    # no. data points   (10^8)
m = 100     # no. of subsamples, subsets    (1000)
            # Make sure n/m is an integer and n/m > k

k = 40      # no. of predictors  (100)
rho = 0.99  # corr. between predictors (0.99)
sigma = 10  # measurement error

beta = c(seq(1, 0.1, by = -0.1), rep(0, k-10))
sig.mat = make.sigma(k, rho)

set.seed(1)
for (i in 1:m){
    X = mvrnorm(n/m, rep(0, k), sig.mat)
    y = rnorm(n/m, mean = X %*% beta, sd = sigma)
    fname = paste0(c(DATA_DIR, "dat_", rho, "_", rep("0", nchar(m) - nchar(i)), i, ".txt"), collapse = "")
    write.table(round(data.frame(y, X), 3), quote = FALSE, row.names = FALSE, file = fname)
    rm(X, y)
    gc()
    }


### Read in data and calculate posterior distributions
### for each data subset
#m = length(list.files(DATA_DIR))
m = 100
fname = paste0(c(DATA_DIR, "dat_", rho, "_", rep("0", nchar(m) - nchar(1)), 1, ".txt"), collapse = "")
z = read.table(fname, header = TRUE)
n = NROW(z)*m
k = NCOL(z) - 1
rm(z)


mu.s = matrix(0, k, m)
#lambda.s = rep(list(matrix(0, k, k)), m)
lambda.s = array(0, c(m, k, k))

a.s = double(m)
b.s = double(m)
for (i in 1:m){
    fname = paste0(c(DATA_DIR, "dat_", rho, "_", rep("0", nchar(m) - nchar(i)), i, ".txt"), collapse = "")
    z = read.table(fname, header = TRUE)
    y = as.vector(z[,1])
    X = as.matrix(z[,-1])

    xty = t(X) %*% y
    yty = sum(y^2)

    lambda.s[i,,] = t(X) %*% X
    lam.inv = solve(t(X) %*% X)
    mu.s[,i] = lam.inv %*% xty
    a.s[i] = (n/m - k) / 2
    b.s[i] = 1/2 * (yty - t(xty) %*% lam.inv %*% xty)

    rm(X, y, z, xty, yty)
    gc()
    }


### Calculate full posterior (sum of the subsets)
# (Algorithm 1)
mu = mu.s[,1]
lambda = lambda.s[1,,]
a = a.s[1]
b = b.s[1]

for (i in 2:m){
    tmp = nig_sum(mu, lambda, a, b, mu.s[,i], lambda.s[i,,], a.s[i], b.s[i])
    mu = tmp$mu
    lambda = tmp$lambda
    a = tmp$a
    b = tmp$b
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
    tmp = nig_sum(mu, lambda, a, b, rep(0, k), lambda.psi, 0, 0)
    mu.bar = tmp$mu
    lambda.bar = tmp$lambda
    a.bar = tmp$a
    b.bar = tmp$b
#   mu.bar = solve(lambda + lambda.psi) %*% (lambda %*% mu + lambda.psi %*% rep(0, k))
#   a.bar = a + 0
#   b.bar = b + 0 + 1/2*t(mu - mu.bar) %*% lambda %*% (mu - mu.bar) +
#       1/2*t(rep(0, k) - mu.bar) %*% lambda.psi %*% (rep(0, k) - mu.bar)
#   lambda.bar = lambda + lambda.psi

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
    tmp = nig_sum(mu, lambda, a, b, rep(0, k), lambda.gamma, 0, 0)
    mu.bar = tmp$mu
    lambda.bar = tmp$lambda
    a.bar = tmp$a
    b.bar = tmp$b

#   mu.bar = solve(lambda + lambda.gamma) %*% (lambda %*% mu + lambda.gamma %*% rep(0, k))
#   a.bar = a + 0
#   b.bar = b + 0 + 1/2*t(mu - mu.bar) %*% lambda %*% (mu - mu.bar) +
#       1/2*t(rep(0, k) - mu.bar) %*% lambda.gamma %*% (rep(0, k) - mu.bar)
#   lambda.bar = lambda + lambda.gamma

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
calc.post = function(gamma, mu, lambda, a, b){
    vg = (gamma == 1)
    vo = (gamma == 0)
    kg = sum(vg)
    ko = sum(vo)

    if (kg == 0 || ko == 0)
        return (-Inf)

    # 4.2, reduce
    sig2 = 1

    mu.g = mu[gamma == 1]
    mu.o = mu[gamma == 0]

    lam.gg = matrix(lambda[vg, vg], kg, kg)
    lam.go = matrix(lambda[vg, vo], kg, ko)
    lam.og = matrix(lambda[vo, vg], ko, kg)
    lam.oo = matrix(lambda[vo, vo], ko, ko)

#   # gamma subscripts
     mu.squig = mu.g + solve(lam.gg) %*% lam.go %*% mu.o
    lam.squig = lam.gg
      a.squig = a + ko / 2
      b.squig = b + 1/2 * t(mu.o) %*% (lam.oo - lam.og %*% solve(lam.gg) %*% lam.go) %*% mu.o
#    mu.squig = mu.g
#   lam.squig = lam.gg
#     a.squig = a
#     b.squig = b

    mu.prior = rep(0, kg)
#   lam.prior = 1/D * diag(kg)
    lam.prior = 1/c * XtX[vg, vg]
    a.prior = 0
    b.prior = 0 
    

    tmp = nig_sum(mu.squig, lam.squig, a.squig, b.squig,
        mu.prior, lam.prior, a.prior, b.prior)
    mu.bar = tmp$mu
    lam.bar = tmp$lambda
    a.bar = tmp$a
    b.bar = tmp$b
    
    return (dnig(rep(0, kg), sig2, mu.prior, solve(lam.prior), a.prior, b.prior) -
        dnig(rep(0, kg), sig2, mu.bar, solve(lam.bar), a.bar, b.bar) +
        sum(dbinom(gamma, 1, prob = w.vec, log = TRUE)))
    }

# d = 0.01^2
D = 100
c = D

#prior.mu = rep(0, k)
#prior.lam = 1/D * diag(k)

# Beta prior for w (mixture weight)
# a.w = rep(1, k)
# b.w = rep(1, k)
#w.vec = rep(0.5, k)


# Gibbs set up
nburn = 5000
nmcmc = 20000
B = nburn + nmcmc
# draws.beta = matrix(0, B, k)
# draws.sig2 = double(B)
draws.gamma = matrix(0, B, k)   # Indicators for selected variables
# draws.w = matrix(0, B, k)       # Probabilities for indicators
accept = double(B)

#draws.sig2[1] = 1
#draws.gamma[1,] = rep(0, k)     # Start with no selected variables
#draws.gamma[1,] = rbinom(k, 1, prob = 0.5)     # Start with randomly selected variables
#draws.gamma[1,] = rep(1, k)     # Start with no selected variables
draws.gamma[1,] = c(1, rep(0, k-1))
# draws.w[1,] = rep(0.5, k)

XtX = apply(lambda.s, c(2,3), sum)

# for AIC
w = 1 / (1 + 1/sqrt(1+c) * exp(2*c/(2*(1+c))))

# for BIC
#w = 1 / (1 + 1/sqrt(1+c) * exp(log(n)*c/(2*(1+c))))

w.vec = rep(w, k)



curr.post = calc.post(draws.gamma[1,], mu, lambda, a, b)


for (i in 2:B){
    draws.gamma[i,] = draws.gamma[i-1,]
    if (i %% 100 == 0)
        cat(i, "/", B, "\r")

    # 4.1, propose gamma
    cand.gamma = draws.gamma[i-1,]
    tmp.ind = sample(k, 1)
    cand.gamma[tmp.ind] = 1 - cand.gamma[tmp.ind]

    cand.post = calc.post(cand.gamma, mu, lambda, a, b)

    if (log(runif(1)) <= cand.post - curr.post){
        draws.gamma[i,] = cand.gamma
        curr.post = cand.post
        accept[i] = 1
        }

    if (i == B)
        cat("\n")
    }

draws.gamma = tail(draws.gamma, nmcmc)
accept = tail(accept, nmcmc)

mean(accept)
draws.gamma[1,]
draws.gamma[7000,]

colMeans(draws.gamma)

