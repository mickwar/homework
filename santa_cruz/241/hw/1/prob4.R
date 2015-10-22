library(MCMCpack)
source("~/files/R/mcmc/bayes_functions.R")

cdf = function(x, t){
    out = length(t)
    for (i in 1:length(t))
        out[i] = mean(x <= t[i])
    return (out)
    }
gen.y1 = function(n)
    rpois(n, 5)
gen.y2 = function(n){
    out = double(n)
    mix = sample(1:2, n, replace = TRUE, prob = c(0.7, 0.3))
    out = ifelse(mix == 1, rpois(sum(mix == 1), 3), out)
    out = ifelse(mix == 2, rpois(sum(mix == 2), 11), out)
    return (out)
    }

n = 300

set.seed(1)
#y = gen.y1(n)
 y = gen.y2(n)

# initial candidate sigma
sig = diag(2)

nj = as.vector(table(y))
ystar = as.numeric(names(table(y))) # sort(unique(y))
nstar = length(ystar)

calc.like = function(alpha, lambda){
    logzpar = function(z, m)
        lgamma(z + m) - lgamma(z)

    nstar * log(alpha) - logzpar(alpha, n) +
        sum(dpois(ystar, lambda, log = TRUE)) + 
        sum(logzpar(alpha * dpois(ystar, lambda, log = FALSE) + 1, nj - 1))
    }

F0.dist.post = function(t, alpha, lambda)
    alpha/(alpha+n)*ppois(t, lambda) + 1/(alpha+n)*sapply(t, function(t) sum(y <= t))

tvec = seq(0, 30, by = 1)

nburn = 10000
nmcmc = 25000
window = 500

fparam = matrix(0, nburn + nmcmc, length(tvec))
al.param = matrix(0, nburn + nmcmc, 2)
al.accept = double(nburn + nmcmc)

a.alpha = 1
b.alpha = 1/10
a.lambda = 5/4
b.lambda = 1/4

### Initial values
al.param[1, 1] = rgamma(1, a.alpha, b.alpha)
al.param[1, 2] = rgamma(1, a.lambda, b.lambda)
fparam[1,] = F0.dist.post(c(tvec[-1], Inf), al.param[1, 1], al.param[1, 2]) - 
    F0.dist.post(tvec, al.param[1, 1], al.param[1, 2])
fparam[1,] = rdirichlet(1, (al.param[1,1] + n)*fparam[1,])

keep.det = double(nburn / window)
keep.acc = double(nburn / window)

### mcmc
par(mfrow = c(1, 1))
for (i in 2:(nburn + nmcmc)){
    cat("\r", i, "/", nburn+nmcmc)
    # update alpha, lambda
    al.param[i,] = al.param[i-1,]
    cand = mvrnorm(1, al.param[i-1,], sig)
    if (all(cand > 0)){
        post = calc.like(al.param[i-1, 1], al.param[i-1, 2]) +
            dgamma(al.param[i-1, 1], a.alpha, b.alpha, log = TRUE) +
            dgamma(al.param[i-1, 2], a.lambda, b.lambda, log = TRUE)
        cand.post = calc.like(cand[1], cand[2]) +
            dgamma(cand[1], a.alpha, b.alpha, log = TRUE) +
            dgamma(cand[2], a.lambda, b.lambda, log = TRUE)
        if (log(runif(1)) <= cand.post - post){
            al.param[i,] = cand
            al.accept[i] = 1
            }
        }
    if ((floor(i/window) == i/window) && (i <= nburn)){
        keep.acc[i/window] = mean(al.accept[(i-window+1):i])
#       sig = (sig + autotune(mean(al.accept[(i-window+1):i]), k = max(1.5, window / 50)) *
#           cov(al.param[(i-window+1):i,])) / 2
        sig = autotune(mean(al.accept[(i-window+1):i]), target = 0.234, k = window/50) *
            (sig + window * var(al.param[(i-window+1):i,]) / i)
        keep.det[i/window] = determinant(sig)$modulus[1]
        plot(keep.det, type='l', pch = 20, ylim = range(keep.det[1:(i/window)]))
        text(1:(nburn/window), keep.det, round(keep.acc, 2))
        }

    # update F (fparam)
    fparam[i,] = F0.dist.post(c(tvec[-1]-0.5, Inf), al.param[i, 1], al.param[i, 2]) - 
        F0.dist.post(tvec-0.5, al.param[i, 1], al.param[i, 2])
    fparam[i,] = rdirichlet(1, (al.param[i, 1] + n)*fparam[i,])
    if (i == (nburn + nmcmc))
        cat("\n")
    }

al.param = tail(al.param, nmcmc)
fparam = tail(fparam, nmcmc)
al.accept = tail(al.accept, nmcmc)

### joint posterior of alpha and lambda, colored by iteration number
plot(al.param, col = rgb(seq(0, 1, length = nmcmc), 0, 0), pch = 20)

qlines = apply(fparam, 2, quantile, c(0, 0.025, 0.5, 0.975, 1))
qlines[3,] = apply(fparam, 2, mean)
cdfparam = apply(fparam, 1, cumsum)

mean(al.accept)

apply(al.param, 2, mean)
var(al.param)

### posterior predictive
temp.1 = sample(y, nmcmc, replace = TRUE)
temp.2 = rpois(nmcmc, al.param[,2])
pred = ifelse(runif(nmcmc) <= al.param[,1] / (al.param[,1] + n), temp.2, temp.1)

### 95% hpd
hpds = apply(al.param, 2, hpd.uni)

### Alpha and Lambda
#pdf("figs/prob4a1.pdf", width = 9, height = 9)
#pdf("figs/prob4b1.pdf", width = 9, height = 9)
par(mfrow = c(2,2), mar = c(4.1, 4.1, 3.1, 1.1))
plot(al.param[,1], type = 'l', main = paste0("Acceptance rate: ", round(mean(al.accept), 3)),
    xlab = "Iteration", ylab = expression(alpha))
dens = density(al.param[,1])
mp = mean(al.param[,1])
hpd.plot(dens, hpds[,1], xlab = expression(alpha),
    main=expression("Marginal posterior for" ~ alpha))
lines(rep(bound(mp, dens), 2), c(0, bound(mp, dens, FALSE)), col = 'red', lwd = 2)
legend("topright", box.lwd = NA, legend = c(NA, paste0("Mean = ", round(mp,2)),
    paste0("Lower HPD = ", round(hpds[1,1], 2)), paste0("Upper HPD = ", round(hpds[2,1], 2))),
    col = c(NA, 'red', col.mult('dodgerblue', 'gray50'), col.mult('dodgerblue', 'gray50')),
    lwd = c(NA, 2, 2, 2))

dens = density(al.param[,2])
mp = mean(al.param[,2])
plot(al.param[,2], type = 'l', main = paste0("Acceptance rate: ", round(mean(al.accept), 3)),
    xlab = "Iteration", ylab = expression(lambda))
hpd.plot(dens, hpds[,2], xlab = expression(lambda),
    main=expression("Marginal posterior for" ~ lambda))
lines(rep(bound(mp, dens), 2), c(0, bound(mp, dens, FALSE)), col = 'red', lwd = 2)
legend("topright", box.lwd = NA, legend = c(NA, paste0("Mean = ", round(mp,2)),
    paste0("Lower HPD = ", round(hpds[1,2], 2)), paste0("Upper HPD = ", round(hpds[2,2], 2))),
    col = c(NA, 'red', col.mult('dodgerblue', 'gray50'), col.mult('dodgerblue', 'gray50')),
    lwd = c(NA, 2, 2, 2))
#dev.off()


### pmf and cdf
#pdf("figs/prob4a2.pdf", height = 9, width = 9)
#pdf("figs/prob4b2.pdf", height = 9, width = 9)
par(mfrow = c(1,1))
plot(tvec, qlines[4,], pch = "_", col = 'darkgreen', cex = 2,
    main = "Posterior p.m.f. from DP", xlab = "y", ylab = "mass")
segments(x0 = tvec, y0 = qlines[2,], y1 = qlines[4,], col = 'darkgreen', lwd = 1)
points(tvec, qlines[3,], pch = 1, col = 'green', lwd = 3)
lines(as.numeric(names(table(pred)))+0.2, table(pred)/nmcmc, type='h', col = 'blue', lwd = 2)
points(tvec, qlines[2,], pch = "_", col = 'darkgreen', cex = 2)
#lines(as.numeric(names(table(y)))+0.3, table(y)/n, col = 'blue', type='h', lwd = 2)
points(as.numeric(names(table(y)))-0.2, table(y)/n, col = 'black', pch = 20, lwd = 2)
legend("topright", border = NA, box.lty = 0,
    legend = c("", "Data (frequency)", "Posterior mean", "95% posterior DP pointwise intervals",
        "Posterior predictive distribution"),
    col = c(NA, "black", "green", "darkgreen", "blue"), pch = c(NA, 20, 1, NA, NA),
    lty = c(NA, NA, NA, 1, 1), lwd = c(NA, NA, 3, 1, 2), cex = 1.3)
#dev.off()

# matplot(tvec, cdfparam, type = 's', lty = 1, col = rgb(0.0, 0.7, 0.0),
#     main = "Posterior c.d.f. from DP", xlab = "y", ylab = "Cumulative mass")
# plot(ecdf(y), add = TRUE, verticals = TRUE, col.01line = NA)
# legend("bottomright", border = NA, box.lty = 0, legend = c("Data", "Posterior draws", ""),
#     col = c("black", rgb(0, 0.7, 0), NA), pch = c(20, NA, NA), lty = c(1, 1, NA),
#     lwd = c(1, 1, NA))

