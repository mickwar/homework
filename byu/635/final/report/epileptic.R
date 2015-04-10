library(MASS)
### glmm on epileptic data
dat = read.table("~/files/data/epileptics.txt", header = TRUE)

n = nrow(dat)
p = 4 # number of visits (and measurements)

# set up the X design matrix
X = matrix(0, n*p, 6)
X[,1] = 1                               # intercept
X[,2] = rep(log(dat$Base/4), each=p)    # Base
X[,3] = rep(dat$Trt, each=p)            # Trt
X[,4] = X[,2] * X[,3]                   # Base * Trt interaction
X[,5] = rep(log(dat$Age), each=p)       # Age
X[,6] = rep(c(0,0,0,1), n)              # Indicator for Visit 4

# set up the Z matrix
# Z = cbind(kronecker(diag(n), rep(1, p)),
#     kronecker(diag(n), seq(-3, 3, by = 2)/10),
#     diag(n*p))
# dim(Z)
Z = cbind(1, rep(seq(-3, 3, by = 2)/10, n))
Z1 = cbind(Z, 1)

# set up Y
Y = as.numeric(t(dat[,2:5]))
Y1 = Y + 1 # for the gamma function instead of factorial


### wishart
tr = function(x)
sum(diag(x))

# multivariate gamma function
# p = 1 is univariate gamma
mgamma = function(a, p)
    pi^(p*(p-1)/4) * prod(gamma(a+(1-1:p)/2))
# log multi gamma
lmgamma = function(a, p)
    p*(p-1)/4*log(pi) + sum(lgamma(a+(1-1:p)/2))

# density of Wishart(V, n) at X
dwishart = function(X, V, df, log = FALSE){
    # X is p x p, positive definite (the support)
    # V is p x p, positive definite
    # df > p - 1 is degrees of freedom
    p = nrow(V)
    if (log){
        return ((df-p-1)/2*determinant(X)$modulus[1] - (1/2)*tr(solve(V) %*% X) - (df*p/2)*log(2) -
            (df/2)*determinant(V)$modulus[1] - lmgamma(df/2, p))
    } else {
        return (det(X)^((df-p-1)/2) * exp(-0.5*tr(solve(V) %*% X)) /
            (2^(df*p/2) * det(V)^(df/2) * mgamma(df/2, p)))
        }
    }

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
    return (L %*% A %*% t(A) %*% t(L))
    }

# alternate method of getting random draws
rwish2 = function(V, df){
    require(MASS)
    p = nrow(V)
    X = mvrnorm(df, rep(0, p), V)
    t(X) %*% X
    }
### end wishart

# assuming indepedence among subjects AND measurements
calc.post = function(alpha, G, e, K){
    require(MASS)
    b12 = matrix(rep(mvrnorm(n, rep(0, 2), G), each=2), n*p, 2)
    b3 = rnorm(n*p, 0, sqrt(e))
    b = cbind(b12, b3)
    lmu = X %*% alpha + apply(Z1 * b, 1, sum)
    mu = exp(lmu)
    D = matrix(0, 3, 3)
    D[1:2,1:2] = G
    D[3,3] = e
    V = Z1 %*% D %*% t(Z1) + kronecker(diag(n), K)

    # likelihood
    out = sum(Y * lmu - mu - lgamma(Y1))
    
    # priors
    out = out - 0.5*determinant(V)$modulus[1] - 0.5*t(lmu - X %*% alpha) %*% (solve(V) %*%
        (lmu - X %*% alpha))
    out = out + dwishart(D, rand_cov_prior, 3, log = TRUE)
    out = out + 0.5*determinant(acovinv)$modulus[1] - 0.5*t(alpha) %*% acovinv %*% alpha
    
    return (out)
    }


nburn = 100
nmcmc = 500

alpha = matrix(0, nburn + nmcmc, ncol(X))
G = rep(list(matrix(0, ncol(Z), ncol(Z))), nburn+nmcmc)
e = matrix(0, nburn+nmcmc, 1)
K = rep(list(matrix(0, p, p)), nburn+nmcmc)

#alpha[1,] = rep(0, ncol(X))
#G[[1]] = diag(2)
#e[1,] = 1
#K[[1]] = matrix(0.5, p, p) + 0.5*diag(p)

alpha[1,] = old.alpha
G[[1]] = old.G
e[1,] = old.e
K[[1]] = old.K
#old.alpha = alpha[600,]
#old.G = G[[600]]
#old.e = e[600,]
#old.K = K[[600]]

a.sig = rep(0.000001, ncol(X))
a.accept = matrix(0, nburn + nmcmc, ncol(X))

g.sig = rep(0.000001, 3)
g.accept = matrix(0, nburn + nmcmc, 3)

e.sig = 0.000001
e.accept = matrix(0, nburn + nmcmc, 1)

k.sig = rep(0.000001, sum(1:p))
k.accept = matrix(0, nburn + nmcmc, length(k.sig))

alpha_cov_prior = 1e3 * solve(t(X) %*% X)
acovinv = solve(alpha_cov_prior)
rand_cov_prior = matrix(0.5, 3, 3) + 0.5*diag(3)

post = calc.post(alpha[1,], G[[1]], e[1,], K[[1]])
keep.post = double(nburn+nmcmc)
keep.post[1] = post

### METROPOLOIS
for (i in 2:(nburn+nmcmc)){
    ### update alpha
    alpha[i,] = alpha[i-1,]
    cand.alpha = alpha[i-1,]
    for (j in 1:ncol(X)){
        cand.alpha[j] = rnorm(1, cand.alpha[j], a.sig[j])
        cand.post = calc.post(cand.alpha, G[[i-1]], e[i-1,], K[[i-1]])
        if (log(runif(1)) < cand.post - post){
            alpha[i,j] = cand.alpha[j]
            a.accept[i,j] = 1
            post = cand.post
        } else {
            cand.alpha[j] = alpha[i-1,j]
            }
        }

    ### update G
    G[[i]] = G[[i-1]]
    cand.G = G[[i-1]]
    # first variance
    cand.G[1,1] = rnorm(1, cand.G[1,1], g.sig[1])
    if (det(cand.G) > 0){
        cand.post = calc.post(alpha[i,], cand.G, e[i-1,], K[[i-1]])
        if (log(runif(1)) < cand.post - post){
            G[[i]][1,1] = cand.G[1,1]
            g.accept[i,1] = 1
            post = cand.post
        } else {
            cand.G[1,1] = G[[i]][1,1]
            }
    } else {
        cand.G[1,1] = G[[i]][1,1]
        }
    # second variance
    cand.G[2,2] = rnorm(1, cand.G[2,2], g.sig[2])
    if (det(cand.G) > 0){
        cand.post = calc.post(alpha[i,], cand.G, e[i-1,], K[[i-1]])
        if (log(runif(1)) < cand.post - post){
            G[[i]][2,2] = cand.G[2,2]
            g.accept[i,2] = 1
            post = cand.post
        } else {
            cand.G[2,2] = G[[i]][2,2]
            }
    } else {
        cand.G[2,2] = G[[i]][2,2]
        }
    # covariance
    cand.G[1,2] <- cand.G[2,1] <-  rnorm(1, cand.G[1,2], g.sig[3])
    if (det(cand.G) > 0){
        cand.post = calc.post(alpha[i,], cand.G, e[i-1,], K[[i-1]])
        if (log(runif(1)) < cand.post - post){
            G[[i]][1,2] <- G[[i]][2,1] <- cand.G[1,2]
            g.accept[i,3] = 1
            post = cand.post
        } else {
            cand.G[1,2] <- cand.G[2,1] <- G[[i]][1,2]
            }
    } else {
        cand.G[1,2] <- cand.G[2,1] <- G[[i]][1,2]
        }

    ### update e
    e[i,] = e[i-1,]
    cand.e = rnorm(1, e[i-1,], e.sig)
    if (cand.e > 0){
        cand.post = calc.post(alpha[i,], G[[i]], cand.e, K[[i-1]])
        if (log(runif(1)) < cand.post - post){
            e[i,] = cand.e
            e.accept[i,1] = 1
            post = cand.post
            }
        }

    ### update K
    K[[i]] = K[[i-1]]
    cand.K = K[[i-1]]
    # variances
    for (j in 1:p){
        cand.K[j,j] = rnorm(1, cand.K[j,j], k.sig[2])
        if (det(cand.K) > 0){
            cand.post = calc.post(alpha[i,], G[[i]], e[i,], cand.K)
            if (log(runif(1)) < cand.post - post){
                K[[i]][j,j] = cand.K[j,j]
                k.accept[i,j] = 1
                post = cand.post
            } else {
                cand.K[j,j] = K[[i]][j,j]
                }
        } else {
            cand.K[j,j] = K[[i]][j,j]
            }
        }
    # covariances
    count = 0
    for (j in 1:(p-1)){
        for (k in (j+1):p){
            count = count + 1
            cand.K[j,k] <- cand.K[k,j] <-  rnorm(1, cand.K[j,k], k.sig[p+count])
            if (det(cand.K) > 0){
                cand.post = calc.post(alpha[i,], G[[i]], e[i,], cand.K)
                if (log(runif(1)) < cand.post - post){
                    K[[i]][j,k] <- K[[i]][k,j] <- cand.K[j,k]
                    k.accept[i,p+count] = 1
                    post = cand.post
                } else {
                    cand.K[j,k] <- cand.K[k,j] <- K[[i]][j,k]
                    }
            } else {
                cand.K[j,k] <- cand.K[k,j] <- K[[i]][j,k]
                }
            }
        }

    keep.post[i] = post
    }


plot(keep.post, type='l')

accept = cbind(a.accept, g.accept, e.accept, k.accept)
apply(accept, 2, mean)
