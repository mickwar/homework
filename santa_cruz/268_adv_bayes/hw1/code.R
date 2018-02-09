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

spikeslab = function(dat, nmcmc, nburn, mh = FALSE){
    require(MASS)
    require(mwBASE)

    y = dat$y
    X = as.matrix(dat$X)
    n = NROW(X)
    p = NCOL(X)
    gj = dat$g.vec      # g_j
    tauj = dat$tau.vec  # tau_j^2

    XtX = t(X) %*% X
    Xty = t(X) %*% y

    # what are good values to choose? similarly for gj and tauj?
    # hyperpriors
    a.w = rep(1, p)
    b.w = rep(1, p)
    a.sig = 1
    b.sig = 1
#   a.w = dat$a.w
#   b.w = dat$b.w
#   a.sig = dat$a.sig
#   b.sig = dat$b.sig

    # +1 is for intercept
    param.beta = matrix(0, nburn + nmcmc, p)
    param.gamma = matrix(0, nburn + nmcmc, p)
    param.w = matrix(0, nburn + nmcmc, p)
    param.sig2 = double(nburn + nmcmc)

    # starting values
    param.w[1,] = 0.5
    param.sig2[1] = 1

    # m-h set up
    if (mh == TRUE){
        cand_sig = diag(1/p, p)
        cand_chol = chol(cand_sig)
        accept = double(nburn + nmcmc)
        calc.beta = function(beta, gamma, sig2, y, X){
            # likelihood
            out = -p/2*log(sig2)-sum((y - X %*% beta)^2) / (2 * sig2)

            # prior
            D.vec = (gamma*(gj-1)+1) * tauj
            out = out -sum(0.5*log(D.vec)) -0.5*sum(D.vec * beta^2)
            return (out)
            }
        window = 200
        }
    
    # sample
    for (i in 2:(nburn + nmcmc)){
        if (floor(i / 100) == i/100)
            cat(i, "/", nburn + nmcmc, "\r")

        # M-H updates are faster, but greater issues with convergence,
        # requiring more samples

        # update beta
        if (mh == TRUE){
            param.beta[i,] = param.beta[i-1,]
            cand = t(param.beta[i-1,] + rnorm(p) %*% cand_chol)
            if (log(runif(1)) <= calc.beta(cand, param.gamma[i-1,], param.sig2[i-1], y, X) -
                calc.beta(param.beta[i-1,], param.gamma[i-1,], param.sig2[i-1], y, X)){
                param.beta[i,] = cand
                accept[i] = 1
                }
            if ((floor(i/window) == i/window) && (i <= nburn)){
                cand_sig = autotune(mean(accept[(i-window+1):i]), target = 0.25, k = window / 50) *
                    (cand_sig + var(param.beta[(i-window+1):i,])/i)
                cand_chol = chol(cand_sig)
                }
        } else {
            Dinv = diag(1/((param.gamma[i-1,]*(gj-1)+1) * tauj))
            sigma = solve(XtX / param.sig2[i-1] + Dinv)
            mu = sigma %*% Xty / param.sig2[i-1]
            param.beta[i,] = mvrnorm(1, mu, sigma)
            }

        # update gamma
        h1 = param.w[i-1,] * dnorm(param.beta[i,], 0, gj*sqrt(tauj))
        h2 = (1-param.w[i-1,]) * dnorm(param.beta[i,], 0, sqrt(tauj))
        param.gamma[i,] = rbinom(p, 1, h1 / (h1 + h2))

        # update w
        param.w[i,] = rbeta(p, a.w + param.gamma[i,], b.w + 1 - param.gamma[i,]) 

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
    X = dat$X
    n = NROW(X)
    p = NCOL(X)
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

                par(mfrow = c(2,2))

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

#               # including the intercept (no, that's dumb)
#               beta.star = c(0, beta)

                # simulate data
                set.seed(1)
                X = mvrnorm(n, rep(0, p), sigma)
                y = rnorm(n, mean = X %*% beta, sd = 1)

                # lasso
                alpha = 1
                cv.lasso = cv.glmnet(X, y, alpha = alpha, intercept = FALSE)
                mod.lasso[[counter]] = glmnet(X, y, lambda = cv.lasso$lambda.min, alpha = alpha,
                    intercept = FALSE)

                # compare betas
#               plot(beta, col = "red", pch = 19, xlab = "predictor index",
#                   ylab = "coefficients", cex.axis = 1.5, cex.lab = 2)
#               lines(coef(mod.lasso[[counter]])[-1], type = "h")

                # prediction at observed data
                y.pred.lasso = predict(mod.lasso[[counter]], newx = X)
#               plot(y, y.pred.lasso)
#               mean((y.pred.lasso - y)^2)/var(y)  # some kind of MSE



                # ridge
                alpha = 0
                cv.ridge = cv.glmnet(X, y, alpha = alpha, intercept = FALSE)
                mod.ridge[[counter]] = glmnet(X, y, lambda = cv.ridge$lambda.min, alpha = alpha,
                    intercept = FALSE)

                # compare betas
#               plot(beta, col = "red", pch = 19, xlab = "predictor index",
#                   ylab = "coefficients", cex.axis = 1.5, cex.lab = 2)
#               lines(coef(mod.ridge[[counter]]), type = "h")
#               lines(coef(mod.lasso[[counter]]), type = "h", col = 'green')

                # prediction at observed data
                y.pred.ridge = predict(mod.ridge[[counter]], newx = X)


#               plot(y.pred.lasso, y.pred.ridge)
#               abline(0, 1)

                # difference on the predictions
#               mean((y.pred.lasso - y.pred.ridge)^2)
#               abline(0, 1)

                # mse's on the beta
                mean((coef(mod.lasso[[counter]])[-1] - beta)^2)
                mean((coef(mod.ridge[[counter]])[-1] - beta)^2)

#               mean((coef(mod.ridge[[16]])[-1] - beta)^2) / var(y)

                mean((y.pred.lasso - y)^2)
                mean((y.pred.lasso - y)^2) / var(y)

                mean((y.pred.ridge - y)^2)
                mean((y.pred.ridge - y)^2) / var(y)


                ### spike and slab
                dat = list("y" = y, "X" = X, "p" = p,
                    "g.vec" = rep(100, p), "tau.vec" = rep(0.01, p))

                mcmc[[counter]] = spikeslab(dat, nmcmc = 2000, nburn = 500)
                preds[[counter]] = get_ppl(dat, mcmc[[counter]])

                Lj = apply(preds[[counter]]$hpd.beta, 2, diff)
                ind = 1:ifelse(i1 == 1, 5, 15)

                mean(Lj[ind])   # M non-zero
                mean(Lj[-ind])  # M zero

                # cbind(apply(mcmc$beta, 2, mean), c(0, beta))
                # apply(mcmc$gamma, 2, mean)
                # apply(mcmc$w, 2, mean)
                keep = which(apply(mcmc[[counter]]$w, 2, mean) > 0.5)
                # mean(mcmc$sig2)

                cols = (apply(mcmc[[counter]]$w, 2, mean) > 0.5)*4+1   # color in the selected variables
                boxplot(mcmc[[counter]]$beta, col = cols, xlim = c(1, 40))
                points(beta, pch = 16, col = 'firebrick')

#               xx = 1:30
#               plot(0, type='n', xlim = range(xx-1), axes = FALSE,
#                   ylim = range(range(coef(mod.lasso[[counter]])), range(coef(mod.ridge[[counter]])),
#                       range(mcmc[[counter]]$beta)), xlab = "Beta Index", ylab = "Value")
#               axis(1); axis(2)
#               points(xx-1, coef(mod.lasso[[counter]])[xx], pch = 16, col = 'green')
#               points(xx-1+0.15, coef(mod.ridge[[counter]])[xx], pch = 16, col = 'blue')
#               points(xx-1+0.30, apply(mcmc[[counter]]$beta[,xx], 2, mean), pch = 16, col = 'red')

                plot(mcmc[[counter]]$beta[,1], type = 'l')
                plot(mcmc[[counter]]$sig2, type = 'l')

                }
            }
        }
    }



pdf("figs/betas.pdf", width = 8, height = 8)
par(mfrow = c(4, 4), mar = c(0,0,0,0), oma = c(4,10,10,4))
xlim = c(1, 20)
ylim = c(-2.5, 5.5)
for (counter in 1:16){

    xx = 2:30
    plot(0, type='n', xlim = xlim, ylim = ylim, axes = FALSE,
        xlab = "Beta Index", ylab = "Value")
    points(xx-1, coef(mod.lasso[[counter]])[xx], pch = 16, col = 'green')
    points(xx-1+0.15, coef(mod.ridge[[counter]])[xx], pch = 16, col = 'blue')
    points(xx-1+0.30, apply(mcmc[[counter]]$beta[,xx-1], 2, mean), pch = 16, col = 'red')

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
dev.off()



sas = matrix(0, 16, 12)
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
                # simulate data
                set.seed(1)
                X = mvrnorm(n, rep(0, p), sigma)
                y = rnorm(n, mean = X %*% beta, sd = 1)



                pred.lasso = predict(mod.lasso[[counter]], newx = X)
                pred.ridge = predict(mod.ridge[[counter]], newx = X)
                pred.spike = colMeans(preds[[counter]]$pred.y)

                Lj = apply(preds[[counter]]$hpd.beta, 2, diff)
                ind = rep(FALSE, ifelse(((counter+1) %/% 2) %% 2, 100, 300))
                ind[1:ifelse(counter <= 8, 5, 15)] = TRUE
                keep = (apply(mcmc[[counter]]$w, 2, mean) > 0.5)

                j = 1
                sas[counter, j] = mean((pred.lasso - y)^2)    # Mean squared prediction error

                j = j + 1
                sas[counter, j] = mean((pred.lasso - y)^2) / var(y)   # Standardize MSPE

                j = j + 1
                sas[counter, j] = mean((pred.ridge - y)^2)

                j = j + 1
                sas[counter, j] = mean((pred.ridge - y)^2) / var(y)

                j = j + 1
                sas[counter, j] = mean((pred.spike - y)^2)

                j = j + 1
                sas[counter, j] = mean((pred.spike - y)^2) / var(y)

                j = j + 1
                sas[counter, j] = mean(Lj[ind])   # M non-zero

                j = j + 1
                sas[counter, j] = mean(Lj[-ind])  # M zero

                j = j + 1
                sas[counter, j] = mean(keep == ind)   # Proportion of correctly selected/non-selected variables

                j = j + 1
                sas[counter, j] = preds[[counter]]$gof

                j = j + 1
                sas[counter, j] = preds[[counter]]$pen

                j = j + 1
                sas[counter, j] = preds[[counter]]$gof + preds[[counter]]$pen
                }
            }
        }
    }

colnames(sas) = c("Lasso_MSE", "Lasso_MSEV", "Ridge_MSE", "Ridge_MSEV",
    "Spike_MSE", "Spike_MSEV", "M_non-zero", "M_zero", "Select", "GOF", "PEN", "PPL")

# Variables correctly selected nearly every time
#sas = sas[-3,]
#sas = sas[, -which(colnames(sas) == "Select")]

xtable(sas)

                counter = counter + 1
for (counter in 1:16){
    ind = rep(FALSE, ifelse(((counter+1) %/% 2) %% 2, 100, 300))
    ind[2:ifelse(counter <= 8, 6, 16)] = TRUE
    keep = (apply(mcmc[[counter]]$w, 2, mean) > 0.5)
    }

### Last bullet point
# set factors
n = 500
p = 100
sigma = make.sigma(p, 0.6)
beta = c(rep(3, 5), rep(0, p-5))   

# simulate data
set.seed(1)
X = mvrnorm(n, rep(0, p), sigma)
y = rnorm(n, mean = X %*% beta, sd = 1)

X.pred = mvrnorm(50, rep(0, p), sigma)
y.pred = rnorm(50, mean = X.pred %*% beta, sd = 1)

# This is combo that corresponds to the scenario of interest
counter = 5
preds[[counter]]$pred.y

params = mcmc[[counter]]
nmcmc = nrow(params$beta)
post.pred = params$beta %*% t(X.pred) + matrix(rnorm(50 * nmcmc, 0, sqrt(params$sig2)), nmcmc, 50)

# MSPE
mean((y.pred - colMeans(post.pred))^2)
