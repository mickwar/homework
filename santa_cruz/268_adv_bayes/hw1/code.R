library(glmnet)
library(MASS)
library(mwBASE)

make.sigma = function(p, rho){
    out = matrix(0, p, p)
    inds = expand.grid(1:p, 1:p)
    fill = rho^apply(inds, 1, function(x) abs(diff(x)))
    out[1:p^2] = fill
    return (out)
    }

spikeslab = function(dat, nmcmc, nburn){
    y = dat$y
    X = cbind(1, dat$X)
    n = length(y)
    p = dat$p
    gj = dat$g.vec      # g_j
    tauj = dat$tau.vec  # tau_j^2

    XtX = t(X) %*% X
    Xty = t(X) %*% y

    # what are good values to choose? similarly for gj and tauj?
    # hyperpriors
    a.w = rep(1, p+1)
    b.w = rep(1, p+1)
    a.sig = 1
    b.sig = 1
#   a.w = dat$a.w
#   b.w = dat$b.w
#   a.sig = dat$a.sig
#   b.sig = dat$b.sig

    # +1 is for intercept
    param.beta = matrix(0, nburn + nmcmc, p+1)
    param.gamma = matrix(0, nburn + nmcmc, p+1)
    param.w = matrix(0, nburn + nmcmc, p+1)
    param.sig2 = double(nburn + nmcmc)

    # starting values
    param.w[1,] = 0.5
    param.sig2[1] = 1
    
    # sample (all are known closed-form posteriors)
    for (i in 2:(nburn + nmcmc)){
        if (floor(i / 100) == i/100)
            cat(i, "/", nburn + nmcmc, "\r")
        # update beta
        Dinv = diag(1/((param.gamma[i-1,]*(gj-1)+1) * tauj))
        sigma = solve(XtX / param.sig2[i-1] + Dinv)
        mu = sigma %*% Xty / param.sig2[i-1]
        param.beta[i,] = mvrnorm(1, mu, sigma)

        # update gamma
        h1 = param.w[i-1,] * dnorm(param.beta[i,], 0, gj*sqrt(tauj))
        h2 = (1-param.w[i-1,]) * dnorm(param.beta[i,], 0, sqrt(tauj))
        param.gamma[i,] = rbinom(p+1, 1, h1 / (h1 + h2))

        # update w
        param.w[i,] = rbeta(p+1, a.w + param.gamma[i,], b.w + 1 - param.gamma[i,]) 

        # update sigma^2
        param.sig2[i] = 1/rgamma(1, a.sig + n/2, b.sig + 0.5*sum((y - X %*% param.beta[i,])^2))

        }
    cat("\n")

    # removed burn-in
    param.beta = tail(param.beta, nmcmc)
    param.gamma = tail(param.gamma, nmcmc)
    param.w = tail(param.w, nmcmc)
    param.sig2 = tail(param.sig2, nmcmc)

    # output
    out = list("beta" = param.beta, "gamma" = param.gamma,
        "w" = param.w, "sig2" = param.sig2)
    return (out)

    }



vec.n = c(500, 400)
vec.p = c(100, 300)
vec.rho = c(0, 0.6)

i1 = 2
i2 = 1
i3 = 2
i4 = 1

for (i1 in 1:2){    # beta
    for (i2 in 1:2){    # rho
        for (i3 in 1:2){    # p
            for (i4 in 1:2){    # n

                # set factors
                n = vec.n[i4]
                p = vec.p[i3]

                sigma = make.sigma(p, vec.rho[i2])
                if (i1 == 1){
                    beta = c(rep(3, 5), rep(0, p-5))   
                } else {
                    beta = c(rep(5, 5), rep(-2, 5),
                        rep(0.5, 5), rep(0, p-15))
                    }

                cat(paste0("beta_type_", i1), vec.rho[i2], vec.p[i3], vec.n[i4], "\n")

                # including the intercept
                beta.star = c(0, beta)

                # simulate data
                set.seed(1)
                X = mvrnorm(n, rep(0, p), sigma)
                y = rnorm(n, mean = X %*% beta, sd = 1)

                # lasso
                alpha = 1
                cv.lasso = cv.glmnet(X, y, alpha = alpha)
                mod.lasso = glmnet(X, y, lambda = cv.lasso$lambda.min, alpha = alpha)

                # compare betas
                plot(c(0, beta), col = "red", pch = 19, xlab = "predictor index",
                    ylab = "coefficients", cex.axis = 1.5, cex.lab = 2)
                lines(coef(mod.lasso), type = "h")

                # prediction at observed data
#               y.pred.lasso1 = predict(cv.lasso, newx=X)
#               mean((y.pred.lasso1 - y)^2)/var(y)

                y.pred.lasso = predict(mod.lasso, newx = X)
#               mean((y.pred.lasso - y)^2)/var(y)  # some kind of MSE

#               plot(y.pred.lasso1, y.pred.lasso)
#               abline(0, 1)


                # ridge
                alpha = 0
                cv.ridge = cv.glmnet(X, y, alpha = alpha)
                mod.ridge = glmnet(X, y, lambda = cv.ridge$lambda.min, alpha = alpha)

                # compare betas
                plot(c(0, beta), col = "red", pch = 19, xlab = "predictor index",
                    ylab = "coefficients", cex.axis = 1.5, cex.lab = 2)
                lines(coef(mod.ridge), type = "h")
                lines(0:p + 1.25, coef(mod.lasso), type = "h", col = 'green')

                # prediction at observed data
                y.pred.ridge = predict(mod.ridge, newx = X)


                plot(y.pred.lasso, y.pred.ridge)
                abline(0, 1)

                # difference on the predictions
                mean((y.pred.lasso - y.pred.ridge)^2)
                abline(0, 1)

                # mse's on the beta
                mean((coef(mod.lasso) - beta.star)^2)
                mean((coef(mod.ridge) - beta.star)^2)


                ### spike and slab
                dat = list("y" = y, "X" = X, "p" = p,
                    "g.vec" = rep(100, p+1), "tau.vec" = rep(0.01, p+1))

                mcmc = spikeslab(dat, nmcmc = 2000, nburn = 500)

                cbind(apply(mcmc$beta, 2, mean), c(0, beta))
                apply(mcmc$gamma, 2, mean)
                apply(mcmc$w, 2, mean)
                keep = which(apply(mcmc$w, 2, mean) > 0.5)
                mean(mcmc$sig2)

                cols = (apply(mcmc$w, 2, mean) > 0.5)*4+1   # color in the selected variables
                boxplot(mcmc$beta, col = cols, xlim = c(1, 40))
                points(beta.star, pch = 16, col = 'firebrick')

                # keep these (first element is intercept)
                which(apply(mcmc$w, 2, mean) > 0.5)

                apply(mcmc$beta, 2, mean)
                apply(mcmc$beta, 2, sd)

                # prediction
                X.star = cbind(1, X)
                keep = 1:101

                mcmc$gamma
                y.pred = mcmc$beta[,keep] %*% t(X.star[,keep])
                mean((apply(y.pred, 2, mean) - y)^2)
                plot(y, apply(y.pred, 2, mean))
                points(y, apply(y.pred, 2, mean), col = 'red')
                abline(0, 1)

                mean((y.pred - matrix(rep(y, each = 2000), 2000))^2)
                # metric?
               


                }
            }
        }
    }
