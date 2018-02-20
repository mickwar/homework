#library(glmnet)
library(MASS)
library(mwBASE)

### Zellner's g-prior
# Model:
# y | X, beta, phi  ~ N(X %*% beta, 1/phi * I)
# beta | phi        ~ N(0, g / phi * solve((t(X) %*% X)))
# p(phi)            = 1/phi
#
# Apparently there is a way to incorporate a variable selection
# with this by conditioning the betas on an indicator variable
# (similar to the spike and slab from HW1).

make.sigma = function(p, rho){
    out = matrix(0, p, p)
    inds = expand.grid(1:p, 1:p)
    fill = rho^apply(inds, 1, function(x) abs(diff(x)))
    out[1:p^2] = fill
    return (out)
    }

nburn = 1000
nmcmc = 5000



for (n in c(50, 200)){
    p = 20
    rho = 0.6
    sigma = make.sigma(p, rho)

    g = max(n, p^2)
    beta = c(rep(3, 5), rep(0, p-5))   

    # simulate data
    set.seed(1)
    X = mvrnorm(n, rep(0, p), sigma)
    y = rnorm(n, mean = X %*% beta, sd = 1)

    ixtx = solve(t(X) %*% X)
    beta.hat = ixtx %*% (t(X) %*% y)

    params.beta = matrix(0, nburn + nmcmc, p)
    params.phi = double(nburn + nmcmc)
    params.phi[1] = 1

    q = g / (1 + g)
    post.beta.mu = g / (1 + g) * beta.hat
    post.beta.sig = g / (1 + g) * beta.hat

    for (i in 2:(nburn + nmcmc)){
        # Update beta
        params.beta[i,] = mvrnorm(1, q * beta.hat,
            q / params.phi[i-1] * ixtx)

        xbeta = X %*% params.beta[i,]
        # Update phi
        params.phi[i] = rgamma(1, (n+p)/2, 0.5*sum((y - xbeta)^2) + 1/(2*g)*sum(xbeta^2))
        }

    params.beta = tail(params.beta, nmcmc)
    params.phi = tail(params.phi, nmcmc)

    }

plot_hpd(params.phi, col1 = 'dodgerblue', main = "Precision",
    axes = FALSE, xlab = expression(phi), ylab = "Density")
axis(1); axis(2)
plot_hpd(1/params.phi, col1 = 'dodgerblue', main = "Variance",
    axes = FALSE, xlab = expression(1/phi), ylab = "Density")
axis(1); axis(2)

boxplot(params.beta, col = c(rep(4, 5), rep(2, 15)), axes = FALSE,
    main = "Coefficients", xlab = "Beta Index")
axis(1); axis(2)
points(beta, col = 'green', pch = 15)



### Random Forest and BART
library(randomForest)

p = 5
n = 1000

set.seed(1)
X = matrix(rnorm(n*p), n, p)
y = 10 * sin(pi * X[,1] * X[,2]) +
    20 * ifelse(X[,2] > 0.5, X[,2] - 0.5, 0)^2 +
    10 * X[,4] + rnorm(n, 0, sd = sqrt(0.5))

B = 1000
pred.rf = matrix(0, B, n)

for (ntree in c(10, 500)){
    for (b in 1:B){
        rf = randomForest(x = X, y = y, ntree = ntree, mtry = sqrt(p))
        pred.rf[b,] = predict(rf)
        }

boxplot(pred.rf)

points(y, col = 'red')

plot(y, colMeans(pred.rf, na.rm = TRUE))
abline(0, 1)

plot(importance(rf))
varImpPlot(rf, n.var = min(p, 10), pch = 15)

