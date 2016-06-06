source("~/files/R/mcmc/bayes_functions.R")
dat = read.table("~/files/data/fye_2016_data.txt", header = TRUE)

### Add 1 (to handle any zeroes if any) and take log transform
dat[,1] = log(dat[,1]+1)
names(dat)[1]="log(ozone+1)"

###
y = dat[,1]
rad = dat[,2]
temp = dat[,3]
wind = dat[,4]
n = length(y)

design.matrices = list(
    cbind("intercept"=1, temp, wind, "temp*wind"=temp*wind),
    cbind("intercept"=1, rad, temp, wind, "temp*wind"=temp*wind))
qnum = length(design.matrices)
cols.A = c('firebrick', 'dodgerblue')
cols.B = c('darkred', 'darkblue')

nburn = 200
nmcmc = 10000

for (model.n in 1:qnum){
    X = design.matrices[[model.n]]

    p = NCOL(X)
    prior.m = rep(0, p)
    prior.g = n
    prior.a = 0
    prior.b = 0

    A = t(X) %*% X
    max(eigen(A)$values) / min(eigen(A)$values)
    chol.A = t(chol(solve(A)))
    post.mean = 1/(prior.g + 1) * (prior.g * solve(t(X) %*% X) %*% t(X) %*% y + prior.m)

    param.beta = matrix(0, nburn + nmcmc, p)
    param.sig2 = double(nburn + nmcmc)
    param.sig2[1] = 1

    for (i in 2:(nburn + nmcmc)){
        param.beta[i,] = post.mean +
            sqrt(prior.g * param.sig2[i-1] / (prior.g + 1)) * chol.A %*% rnorm(p, 0, 1)
        param.sig2[i] = 1/rgamma(1,
            prior.a + n/2 + p/2, 
            prior.b + 0.5*sum((y - X %*% param.beta[i,])^2) +
                0.5/prior.g * t(param.beta[i,] - prior.m) %*% A %*% (param.beta[i,] - prior.m))
        }

    # Truncate
    param.beta = tail(param.beta, nmcmc)
    param.sig2 = tail(param.sig2, nmcmc)

    ### Posterior predictions
    pred.y = matrix(0, n, nmcmc)
    for (i in 1:nmcmc)
        pred.y[,i] = rnorm(n, X %*% param.beta[i,], 1*sqrt(param.sig2[i]))

    hpd.beta = apply(param.beta, 2, hpd.uni)
    hpd.sig2 = hpd.uni(param.sig2)

#   hpd.plot(density(param.sig2), hpd.sig2, main = expression("Model Variance"),
#       xlab = expression(sigma^2))
    pdf(paste0("figs/post_", model.n, ".pdf"), width = 3*ceiling(sqrt(p)), height = 3*floor(sqrt(p)))
    par(mfrow = c(floor(sqrt(p)), ceiling(sqrt(p))), mar = c(4.1, 2.1, 2.1, 1.1),
        oma = c(0, 0, 2, 0))
    for (i in 1:p){
        hpd.plot(density(param.beta[,i]), hpd.beta[,i], main = colnames(X)[i],
            xlab = bquote(beta[.(i-1)]), col1 = cols.A[model.n], border = NA,
            axes = FALSE)
        axis(1); axis(2)
        abline(v = 0, col = 'forestgreen', lwd = 2)
        }
    title(main = paste0("M", model.n, " marginal posteriors"), outer = TRUE)
    par(mfrow = c(1,1), mar = c(5.1, 4.1, 4.1, 2.1), oma = c(0,0,0,0))
    dev.off()

    mm = apply(pred.y, 1, mean)
    qq = apply(pred.y, 1, quantile, c(0.025, 0.975))
    pdf(paste0("figs/pred_obsfit_", model.n, ".pdf"), width = 6, height = 6)
    plot(y, mm, pch = 20, ylim = range(qq), col = cols.B[model.n], xlab = "Observed",
        ylab = "Predicted", main = paste0("M", model.n, " model predictions"),
        cex.lab = 1.3, axes = FALSE)
    abline(0, 1)
    axis(1); axis(2)
    segments(x0 = y, y0 = qq[1,], x1 = y, y1 = qq[2,], col = cols.A[model.n])
    points(y, mm, pch = 20, col = cols.B[model.n])
    dev.off()
    }

