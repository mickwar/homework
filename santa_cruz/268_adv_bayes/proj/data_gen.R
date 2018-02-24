library(MASS)
library(mwBASE)

DATA_DIR = "~/files/tmp/268_data/"

### For generating X
make.sigma = function(k, rho){
    out = matrix(rho, k, k)
    diag(out) = 1
    return (out)
    }


### Generate and write the data to m files
n = 10^6    # no. data points   (10^8)
m = 100     # no. of subsamples, subsets    (1000)
            # Make sure n/m is an integer and n/m > k

k = 40     # no. of predictors  (100)
rho = 0.99  # corr. between predictors
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
    }


### Calculate full posterior (sum of the subsets)
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


lambda[1:6,1:6]

B = 10000
plot_hpd(1/rgamma(B, a, b), axes = FALSE); axis(1); axis(2)
abline(v = sigma^2)
