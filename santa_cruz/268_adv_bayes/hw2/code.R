library(xtable)
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

out = matrix(0, 2, 4)
counter = 0

for (n in c(50, 200)){
    counter = counter + 1
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

    out[counter, 1:2] = quantile(params.beta[,1], c(0.025, 0.975))
    out[counter, 3:4] = quantile(params.beta[,10], c(0.025, 0.975))
    }


xtable(out)



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
library(BART)

p = c(200, 100)
n = c(200, 500)
ntree = c(10, 500)

#out = matrix(0, 16, 7)
out = data.frame(matrix(0, 16, 8))
counter = 0

for (i0 in 1:2){
    for (i1 in 1:2){
        for (i2 in 1:2){
            for (i3 in 1:2){

                counter = counter + 1

                out[counter, 2] = n[i1]
                out[counter, 3] = p[i1]
                out[counter, 4] = ntree[i2]
                out[counter, 5] = ""

                set.seed(1)
                X = matrix(rnorm(n[i1]*p[i1]), n[i1], p[i1])
                y = 10 * sin(pi * X[,1] * X[,2]) +
                    20 * ifelse(X[,2] > 0.5, X[,2] - 0.5, 0)^2 +
                    10 * X[,4] + rnorm(n[i1], 0, sd = sqrt(0.5))


                if (i3 == 2){
                    out[counter, 5] = "Noise"
                    # With extra noise on 1, 2, and 9
                    X[,1] = X[,1] + rnorm(n, 0, sqrt(0.1))
                    X[,2] = X[,2] + rnorm(n, 0, sqrt(0.1))
                    X[,9] = X[,9] + rnorm(n, 0, sqrt(0.1))
                    }

                if (i0 == 1){
                    out[counter, 1] = "RF"
                    # Random Forest
                    rf = randomForest(x = X, y = y, ntree = ntree[i2], mtry = sqrt(p[i1]))
                    pred.rf = predict(rf, X, predict.all = TRUE)

                    ints = apply(pred.rf$individual, 1, quantile, c(0.025, 0.975))
                    means = rowMeans(pred.rf$individual)
                    coverage = mean(y > ints[1,] & y < ints[2,])
                    lengths = mean(apply(ints, 2, diff))
                    mspe = mean((y - means)^2)
                } else {
                    out[counter, 1] = "BART"
                    # BART
                    ba = wbart(X, y, ntree = ntree[i2], nskip = 1000, ndpost = 2000)
                    pred.ba = predict(ba, X)

                    ints = apply(pred.ba, 2, quantile, c(0.025, 0.975))
                    means = colMeans(pred.ba)
                    coverage = mean(y > ints[1,] & y < ints[2,])
                    lengths = mean(apply(ints, 2, diff))
                    mspe = mean((y - means)^2)
                    }

                plot(y, means)
                segments(y, ints[1,], y, ints[2,])
                abline(0, 1)

                out[counter, 6] = coverage
                out[counter, 7] = lengths
                out[counter, 8] = mspe
                }
            }
        }
    }

colnames(out) = c("Model", "n", "p", "ntree", "Noise", "Coverage", "Length", "MSPE")
out

xtable(out)
