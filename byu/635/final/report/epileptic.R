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
Z = cbind(kronecker(diag(n), rep(1, p)),
    kronecker(diag(n), seq(-3, 3, by = 2)/10),
    diag(n*p))
dim(Z)

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
calc.post = function(alpha, Ginv){
    mu = 
    X 
    }


nburn = 100
nmcmc = 500

alpha = matrix(0, nburn + nmcmc, 6)
Ginv = rep(list(matrix(0, ncol(Z), ncol(Z))), nburn+nmcmc)
alpha = matrix(0, nburn + nmcmc, 6)
