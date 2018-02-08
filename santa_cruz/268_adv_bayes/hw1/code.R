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
    require(MASS)
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

get_ppl = function(dat, params, plot = TRUE){
    y = dat$y
    X = cbind(1, dat$X)
    n = nrow(X)
    nmcmc = nrow(params$beta)

    # Posterior predictive for each x
    pred.y = params$beta %*% t(X) + matrix(rnorm(n * nmcmc, 0, sqrt(params$sig2)), nmcmc, n)

    # hpds for the predictions
    mm = apply(pred.y, 2, mean)
    hpds.y = apply(pred.y, 2, hpd_mult, force_uni = TRUE)
    if (plot){
        plot(y, mm, type = 'n', ylim = range(hpds.y))
        segments(x0 = y, y0 = hpds.y[1,], y1 = hpds.y[2,], col = 'dodgerblue')
        points(y, mm, pch = '-', col = 'firebrick', cex = 1.5)
        abline(0, 1)
        }

    # hpds for betas
    hpds.beta = apply(params$beta, 2, hpd_mult, force_uni = TRUE)

    # posterior predictive loss
    # goodness of fit and penalty
    gof = sum((y - mm)^2)
    pen = sum(apply(pred.y, 2, var))

    return (list("pred.y" = pred.y, "gof" = gof, "pen" = pen,
        "hpd.y" = hpds.y, "hpd.beta" = hpds.beta))
    }



vec.n = c(500, 400)
vec.p = c(100, 300)
vec.rho = c(0, 0.6)

mod.lasso = rep(list(NULL), 16)
mod.ridge = rep(list(NULL), 16)
mcmc = rep(list(NULL), 16)
preds = rep(list(NULL), 16)
counter = 0

for (i1 in 1:2){    # beta
    for (i2 in 1:2){    # rho
        for (i3 in 1:2){    # p
            for (i4 in 1:2){    # n

                counter = counter + 1

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
                mod.lasso[[counter]] = glmnet(X, y, lambda = cv.lasso$lambda.min, alpha = alpha)

                # compare betas
                plot(c(0, beta), col = "red", pch = 19, xlab = "predictor index",
                    ylab = "coefficients", cex.axis = 1.5, cex.lab = 2)
                lines(coef(mod.lasso[[counter]]), type = "h")

                # prediction at observed data
#               y.pred.lasso1 = predict(cv.lasso, newx=X)
#               mean((y.pred.lasso1 - y)^2)/var(y)

                y.pred.lasso = predict(mod.lasso[[counter]], newx = X)
#               mean((y.pred.lasso - y)^2)/var(y)  # some kind of MSE

#               plot(y.pred.lasso1, y.pred.lasso)
#               abline(0, 1)


                # ridge
                alpha = 0
                cv.ridge = cv.glmnet(X, y, alpha = alpha)
                mod.ridge[[counter]] = glmnet(X, y, lambda = cv.ridge$lambda.min, alpha = alpha)

                # compare betas
                plot(c(0, beta), col = "red", pch = 19, xlab = "predictor index",
                    ylab = "coefficients", cex.axis = 1.5, cex.lab = 2)
                lines(coef(mod.ridge[[counter]]), type = "h")
                lines(0:p + 1.25, coef(mod.lasso[[counter]]), type = "h", col = 'green')

                # prediction at observed data
                y.pred.ridge = predict(mod.ridge[[counter]], newx = X)


                plot(y.pred.lasso, y.pred.ridge)
                abline(0, 1)

                # difference on the predictions
                mean((y.pred.lasso - y.pred.ridge)^2)
                abline(0, 1)

                # mse's on the beta
                mean((coef(mod.lasso[[counter]]) - beta.star)^2)
                mean((coef(mod.ridge[[counter]]) - beta.star)^2)


                ### spike and slab
                dat = list("y" = y, "X" = X, "p" = p,
                    "g.vec" = rep(100, p+1), "tau.vec" = rep(0.01, p+1))

                mcmc[[counter]] = spikeslab(dat, nmcmc = 2000, nburn = 500)
                preds[[counter]] = get_ppl(dat, mcmc[[counter]])

                Lj = apply(preds[[counter]]$hpds, 2, diff)
                ind = 2:ifelse(i1 == 1, 6, 16)

                mean(Lj[ind])   # M non-zero
                mean(Lj[-ind])  # M zero

                # cbind(apply(mcmc$beta, 2, mean), c(0, beta))
                # apply(mcmc$gamma, 2, mean)
                # apply(mcmc$w, 2, mean)
                keep = which(apply(mcmc[[counter]]$w, 2, mean) > 0.5)
                # mean(mcmc$sig2)

                cols = (apply(mcmc[[counter]]$w, 2, mean) > 0.5)*4+1   # color in the selected variables
                boxplot(mcmc[[counter]]$beta, col = cols, xlim = c(1, 40))
                points(beta.star, pch = 16, col = 'firebrick')

                xx = 1:30
                plot(0, type='n', xlim = range(xx-1), axes = FALSE,
                    ylim = range(range(coef(mod.lasso[[counter]])), range(coef(mod.ridge[[counter]])),
                        range(mcmc[[counter]]$beta)), xlab = "Beta Index", ylab = "Value")
                axis(1); axis(2)
                points(xx-1, coef(mod.lasso[[counter]])[xx], pch = 16, col = 'green')
                points(xx-1+0.15, coef(mod.ridge[[counter]])[xx], pch = 16, col = 'blue')
                points(xx-1+0.30, apply(mcmc[[counter]]$beta[,xx], 2, mean), pch = 16, col = 'red')

                }
            }
        }
    }



par(mfrow = c(4, 4), mar = c(0,0,0,0), oma = c(4,10,10,4))
xlim = c(0, 25)
ylim = c(-2.5, 5.5)
for (counter in 1:16){

    xx = 1:30
    plot(0, type='n', xlim = xlim, ylim = ylim, axes = FALSE,
        xlab = "Beta Index", ylab = "Value")
    points(xx-1, coef(mod.lasso[[counter]])[xx], pch = 16, col = 'green')
    points(xx-1+0.15, coef(mod.ridge[[counter]])[xx], pch = 16, col = 'blue')
    points(xx-1+0.30, apply(mcmc[[counter]]$beta[,xx], 2, mean), pch = 16, col = 'red')

    # Draw axes
    if ((counter <= 4) && (counter %% 2 == 1))
        axis(3, lwd = -1, lwd.ticks = 1)
    if ((counter >= 13) && (counter %% 2 == 0))
        axis(1, lwd = -1, lwd.ticks = 1)
    if (counter %% 4 == 0)
        axis(4)
    box()

    # Put in thicker lines in the middle
    if (counter >= 5 && counter <= 8)
        axis(1, labels = FALSE, lwd.ticks = 0, lwd = 5,
            at = xlim + c(-1, 1)*diff(xlim)/5)
    if ((counter %% 4) == 2)
        axis(4, labels = FALSE, lwd.ticks = 0, lwd = 5,
            at = ylim + c(-1, 1)*diff(ylim)/5)
    }
cex = 1.25
mtext(c("p=100", "p=300"), side = 3, outer = TRUE, at = c(0.25, 0.75),
    line = 7, cex = cex)
mtext(c("n=500", "n=400"), side = 3, outer = TRUE, at = c(0.125, 0.375, 0.625, 0.875),
    line = 4, cex = cex)
mtext(c("BetaType1", "BetaType2"), side = 2, outer = TRUE, at = rev(c(0.25, 0.75)),
    line = 9.5, las = 1, cex = cex, adj = 0)
mtext(c("Rho=0", "Rho=0.6"), side = 2, outer = TRUE, at = rev(c(0.125, 0.375, 0.625, 0.875)),
    line = 6.3, las = 1, cex = cex, adj = 0)
mtext(expression(beta[j]), outer = TRUE, side = 3, at = -0.120, line = 4, cex = 2.5, adj = 0)



sas = matrix(0, 6, 16)
for (counter in 1:16){
    Lj = apply(preds[[counter]]$hpd.beta, 2, diff)
    ind = rep(FALSE, ifelse(((counter+1) %/% 2) %% 2, 100, 300)+1)
    ind[2:ifelse(counter <= 8, 6, 16)] = TRUE

    keep = (apply(mcmc[[counter]]$w, 2, mean) > 0.5)

    sas[1, counter] = mean(Lj[ind])   # M non-zero
    sas[2, counter] = mean(Lj[-ind])  # M zero
    sas[3, counter] = mean(keep == ind)   # Proportion of correctly selected/non-selected variables
    sas[4, counter] = preds[[counter]]$gof
    sas[5, counter] = preds[[counter]]$pen
    sas[6, counter] = preds[[counter]]$gof + preds[[counter]]$pen
    }

rownames(sas) = c("M_non-zero", "M_zero", "Select", "GOF", "PEN", "PPL")

# Variables correctly selected nearly every time
sas = sas[-3,]

sas
