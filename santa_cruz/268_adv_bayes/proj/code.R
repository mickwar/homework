library(MASS)
library(mwBASE)

DATA_DIR = "~/files/tmp/268_data/"

### Inverse gamma density
digamma = function(x, a, b, log = TRUE){
    out = a*log(b) -lgamma(a)- (a+1)*log(x) - b/x
    if (log)
        return (out)
    return (exp(out))
    }
dmvnorm = function(x, mu, sigma, log = TRUE){
    out = -length(x)/2*log(2*pi) - 0.5*determinant(sigma)$modulus[1] - 0.5*t(x-mu)%*%solve(sigma)%*%(x-mu)
    if (log)
        return (out)
    return (exp(out))
    }
### For generating X
make.sigma = function(k, rho){
    out = matrix(rho, k, k)
    diag(out) = 1
    return (out)
    }
### Normal Inverse Gamma density
dnig = function(beta, sig2, mu, lambda, a, b, log = TRUE){
    # lambda is precision
    k = length(mu)
    U = chol(lambda)
    Uinv = t(backsolve(U, diag(k)))

    t(beta-mu) %*% lambda %*% (beta-mu)
    t(beta-mu) %*% t(U) %*% U %*% (beta-mu)
    sum((U %*% (beta-mu))^2)
    determinant(U)
    sum(log(diag(U)))

    t(beta-mu) %*% solve(lambda) %*% (beta-mu)
    t(beta-mu) %*% t(Uinv) %*% Uinv %*% (beta-mu)
    sum((Uinv %*% (beta-mu))^2)
    determinant(Uinv)
    sum(log(diag(Uinv)))

    out = -k/2*log(2*pi) + a*log(b) - lgamma(a) -(a+k/2+1)*log(sig2) -
        (1/sig2)*(b+1/2*sum((U %*% (beta - mu))^2)) -sum(log(diag(Uinv)))
    if (log)
        return (out)
    return (exp(out))
    }
### Samples for NIG
rnig = function(n, mu, lambda, a, b){
    k = length(mu)
    U = backsolve(chol(lambda), diag(k))

    sigma2 = 1/rgamma(n, a, b)
#   mu + L %*% rnorm(k, 0, 1)
    tapply(sigma2, 1:n, function(s) t(mu + (sqrt(s)*U) %*% rnorm(k)))
    beta = t(sapply(sigma2, function(s) mu + (sqrt(s)*U) %*% rnorm(k)))
    return (list("beta"=beta, "sigma2"=sigma2))
    }
### Random Inverse Gaussian
rigauss = function(n, mu, lambda){
    nu = rnorm(n)
    y = nu^2
    x = mu + y*mu^2/(2*lambda) - mu/(2*lambda) *
        sqrt(4*mu*lambda*y + mu^2*y^2)
    # account of when mu is large (becomes a levy) (happens when x = 0)
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
    return (list("mu"=mu, "lambda"=lambda, "a"=a, "b"=as.numeric(b)))
    }

dnig(5, 3, 4, 1/0.2, 6, 7)
dnig(5, 3, 4, 0.2, 6, 7)

dnorm(5, 4, sqrt(3 * 0.2), log = TRUE) + digamma(3, 6, 7)

A = matrix(c(5, 4, 4, 7), 2, 2)

dnig(c(4.8, 4.1), 3, c(4, 3), A, 6, 7)
dnig(c(4.8, 4.1), 3, c(4, 3), solve(A), 6, 7)

dmvnorm(c(4.8, 4.1), c(4, 3), 3*solve(A), log = TRUE) + digamma(3, 6, 7)
dmvnorm(c(4.8, 4.1), c(4, 3), 1/3*A, log = TRUE) + digamma(3, 6, 7)



xx = seq(0.001, 100, length = 200)
plot(xx, digamma(xx, 5, 100, log = FALSE), type = 'l')
lines(density(1/rgamma(10000, 5, 100)), col = 'red')



### Generate and write the data to m files
n = 10^6    # no. data points   (10^8)
m = 100     # no. of subsamples, subsets    (1000)
            # Make sure n/m is an integer and n/m > k

k = 40      # no. of predictors  (100)
rho = 0.99  # corr. between predictors (0.99)
sigma = 10  # measurement error

beta = c(seq(1, 0.1, by = -0.1), rep(0, k-10))
sig.mat = make.sigma(k, rho)

# set.seed(1)
# for (i in 1:m){
#     X = mvrnorm(n/m, rep(0, k), sig.mat)
#     y = rnorm(n/m, mean = X %*% beta, sd = sigma)
#     fname = paste0(c(DATA_DIR, "dat_", rho, "_", rep("0", nchar(m) - nchar(i)), i, ".txt"), collapse = "")
#     write.table(round(data.frame(y, X), 3), quote = FALSE, row.names = FALSE, file = fname)
#     rm(X, y)
#     gc()
#     }


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
lambda.s = array(0, c(m, k, k))
a.s = double(m)
b.s = double(m)

# mu.2 = matrix(0, kg, m)
# lambda.2 = array(0, c(m, kg, kg))
# a.2 = double(m)
# b.2 = double(m)

for (i in 1:m){
    fname = paste0(c(DATA_DIR, "dat_", rho, "_", rep("0", nchar(m) - nchar(i)), i, ".txt"), collapse = "")
    z = read.table(fname, header = TRUE)
    y = as.vector(z[,1])
    X = as.matrix(z[,-1])

#   X = X[,vg]
    xty = t(X) %*% y
    yty = sum(y^2)

    lambda.s[i,,] = t(X) %*% X
    lam.inv = solve(t(X) %*% X)
    mu.s[,i] = lam.inv %*% xty
    a.s[i] = (n/m - k) / 2
    b.s[i] = 1/2 * (yty - t(xty) %*% lam.inv %*% xty)

#   lambda.2[i,,] = t(X) %*% X
#   lam.inv = solve(t(X) %*% X)
#   mu.2[,i] = lam.inv %*% xty
#   a.2[i] = (n/m - k) / 2
#   b.2[i] = 1/2 * (yty - t(xty) %*% lam.inv %*% xty)

    rm(X, y, z, xty, yty, lam.inv)
    gc()
    }


### Calculate full posterior (sum of the subsets)
# (Algorithm 1)
mu = mu.s[,1]
lambda = lambda.s[1,,]
a = a.s[1]
b = b.s[1]

# mu.3 = mu.2[,1]
# lambda.3 = lambda.2[1,,]
# a.3 = a.2[1]
# b.3 = b.2[1]

for (i in 2:m){
    tmp = nig_sum(mu, lambda, a, b, mu.s[,i], lambda.s[i,,], a.s[i], b.s[i])
    mu = tmp$mu
    lambda = tmp$lambda
    a = tmp$a
    b = tmp$b
#   tmp = nig_sum(mu.3, lambda.3, a.3, b.3, mu.2[,i], lambda.2[i,,], a.2[i], b.2[i])
#   mu.3 = tmp$mu
#   lambda.3 = tmp$lambda
#   a.3 = tmp$a
#   b.3 = tmp$b
    rm(tmp)
    }
# Non-informative prior
alg1 = rnig(1000, mu, lambda, a, b)




### LASSO
# (Algorithm 2, start with Algorithm 1)

# Penalty terms
pen.vec = rep(10, k)

# Gibbs set up
nburn = 2000
nmcmc = 5000
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


# mean.beta = colMeans(draws.beta)
# med.beta = apply(draws.beta, 2, median)
# qq.beta = apply(draws.beta, 2, quantile, c(0.025, 0.975))
# 
# plot(mean.beta, pch = 16, col = 'darkgreen', ylim = range(qq.beta))
# points(1:k, med.beta, pch = 16, col = 'dodgerblue')
# segments(1:k, qq.beta[1,], 1:k, qq.beta[2,], col = 'forestgreen')
# 
# points(1:k, beta, pch = 15)

alg2 = list("beta"=draws.beta, "sigma2"=draws.sig2, "psi"=draws.psi)



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
draws.gamma = tail(draws.gamma, nmcmc)
draws.w = tail(draws.w, nmcmc)

# colMeans(draws.w)
# 
# 
# cols = rep("white", k)
# cols[which(colMeans(draws.w) > 0.5)] = "dodgerblue"
# boxplot(draws.beta, col = cols)
# points(1:k, beta, col = 'firebrick', pch = 15)
# 
# plot_hpd(draws.sig2)

alg3 = list("beta"=draws.beta, "sigma2"=draws.sig2, "gamma"=draws.gamma, "w"=draws.w)




### MC^3
# (Algorithm 4)
calc.post = function(gamma, mu, lambda, a, b){
    vg = (gamma == 1)
    vo = (gamma == 0)
    kg = sum(vg)
    ko = sum(vo)

    if (kg == 0 || ko == 0)
        return (-Inf)

    sig2 = b/(a-1)

    mu.g = mu[vg]
    mu.o = mu[vo]

    lam.gg = matrix(lambda[vg, vg], kg, kg)
    lam.go = matrix(lambda[vg, vo], kg, ko)
    lam.og = matrix(lambda[vo, vg], ko, kg)
    lam.oo = matrix(lambda[vo, vo], ko, ko)

#   # gamma subscripts
     mu.squig = mu.g + solve(lam.gg) %*% lam.go %*% mu.o
    lam.squig = lam.gg
#     a.squig = a + ko / 2
      a.squig = a + (1-m) * ko / 2
      b.squig = b + 1/2 * t(mu.o) %*% (lam.oo - lam.og %*% solve(lam.gg) %*% lam.go) %*% mu.o
#    mu.squig = mu.g
#   lam.squig = lam.gg
#     a.squig = a
#     b.squig = b


    mu.prior = rep(0, kg)
#   lam.prior = 1/D * diag(kg)
    lam.prior = 1/c * lambda[vg, vg]
    a.prior = 0
    b.prior = 0 
    

    tmp = nig_sum(mu.squig, lam.squig, a.squig, b.squig,
        mu.prior, lam.prior, a.prior, b.prior)
    mu.bar = tmp$mu
    lam.bar = tmp$lambda
    a.bar = tmp$a
    b.bar = tmp$b
    
#   return (dnig(rep(0, kg), sig2, mu.prior, lam.prior, a.prior, b.prior) -
#       dnig(rep(0, kg), sig2, mu.bar, lam.bar, a.bar, b.bar) +
#       sum(dbinom(gamma, 1, prob = w.vec, log = TRUE)))
#   return (dnig(rep(0, kg), sig2, mu.squig, lam.squig, a.squig, b.squig) -
#       dnig(rep(0, kg), sig2, mu.bar, lam.bar, a.bar, b.bar) +
#       sum(dbinom(gamma, 1, prob = w.vec, log = TRUE)))
    return (-log(sig2) - dnig(rep(0, kg), sig2, mu.bar, lam.bar, a.bar, b.bar) +
        sum(dbinom(gamma, 1, prob = w.vec, log = TRUE)))
    }

alg4.draws = function(gammas, mu, lambda, a, b){
    B = NROW(gammas)
    k = NCOL(gammas)
    draws.beta = matrix(0, B, k)
    draws.sig2 = matrix(0, B, k)
    for (b in 1:B){
        vg = (gammas[b,] == 1)
        vo = (gammas[b,] == 0)
        kg = sum(vg)
        ko = sum(vo)

        if (kg == 0 || ko == 0)
            return (-Inf)

        sig2 = b/(a-1)

        mu.g = mu[vg]
        mu.o = mu[vo]

        lam.gg = matrix(lambda[vg, vg], kg, kg)
        lam.go = matrix(lambda[vg, vo], kg, ko)
        lam.og = matrix(lambda[vo, vg], ko, kg)
        lam.oo = matrix(lambda[vo, vo], ko, ko)

    #   # gamma subscripts
         mu.squig = mu.g + solve(lam.gg) %*% lam.go %*% mu.o
        lam.squig = lam.gg
    #     a.squig = a + ko / 2
          a.squig = a + (1-m) * ko / 2
          b.squig = b + 1/2 * t(mu.o) %*% (lam.oo - lam.og %*% solve(lam.gg) %*% lam.go) %*% mu.o
    #    mu.squig = mu.g
    #   lam.squig = lam.gg
    #     a.squig = a
    #     b.squig = b


        mu.prior = rep(0, kg)
    #   lam.prior = 1/D * diag(kg)
        lam.prior = 1/c * lambda[vg, vg]
        a.prior = 0
        b.prior = 0 
        

        tmp = nig_sum(mu.squig, lam.squig, a.squig, b.squig,
            mu.prior, lam.prior, a.prior, b.prior)
        mu.bar = tmp$mu
        lam.bar = tmp$lambda
        a.bar = tmp$a
        b.bar = tmp$b

        tmp = rnig(1, mu.bar, lam.bar, a.bar, b.bar)
        draws.beta[b,vg] = tmp$beta
        draws.sig2[b] = tmp$sigma2
        }
    return (list("beta"=draws.beta, "sigma2"=draws.sig2, "gamma"=gammas))
    }
    

# d = 0.01^2
D = 100
c = D



# Gibbs set up
nburn = 5000
nmcmc = 20000
B = nburn + nmcmc
draws.gamma = matrix(0, B, k)   # Indicators for selected variables
accept = double(B)

draws.gamma[1,] = rbinom(k, 1, prob = 0.5)     # Start with randomly selected variables

#XtX = apply(lambda.s, c(2,3), sum)

# for lazy
#w = 0.5

# for AIC
#w = 1 / (1 + 1/sqrt(1+c) * exp(2*c/(2*(1+c))))

# for BIC
w = 1 / (1 + 1/sqrt(1+c) * exp(log(n)*c/(2*(1+c))))

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
colMeans(draws.gamma)

alg4 = alg4.draws(draws.gamma, mu, lambda, a, b)


### Comine plot
MM = array(0, c(4, k))
QQ = array(0, c(4, 2, k))

MM[1,] = colMeans(alg1$beta)
QQ[1,,] = apply(alg1$beta, 2, hpd_mult, force_uni = TRUE)
MM[2,] = colMeans(alg2$beta)
QQ[2,,] = apply(alg2$beta, 2, hpd_mult, force_uni = TRUE)
MM[3,] = colMeans(alg3$beta)
QQ[3,,] = apply(alg3$beta, 2, hpd_mult, force_uni = TRUE)
MM[4,] = colMeans(alg4$beta)
QQ[4,,] = apply(alg4$beta, 2, hpd_mult, force_uni = TRUE)

plot(mm1, ylim = range(qq1, qq2), col = 'blue', pch = 16)
points(mm2, col = 'green', pch = 16)
points(mm3, col = 'red', pch = 16)
points(mm4, col = 'red', pch = 16)

lines(qq1[1,], lty = 2, col = 'blue')
lines(qq1[2,], lty = 2, col = 'blue')
lines(qq2[1,], lty = 2, col = 'green')
lines(qq2[2,], lty = 2, col = 'green')
lines(qq3[1,], lty = 2, col = 'red')
lines(qq3[2,], lty = 2, col = 'red')
lines(qq4[1,], lty = 2, col = 'red')
lines(qq4[2,], lty = 2, col = 'red')


plot(mm1, ylim = range(qq1,qq2,qq3,qq4), lty = 2, type = 'l', col =)

matplot(t(MM), type = 'l')


