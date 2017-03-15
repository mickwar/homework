library(mwBASE)
dat = read.table("./california_pr_obs19501999.txt", header = TRUE)

# Get December, January, and February (but not leap day)
dat.ind = grep("(-12-)|(-01-)|(-02-(?!29))", dat[,1], perl = TRUE)

n = length(dat[,2])
ndays = 90  # Number of winter days
nyears = length(dat.ind) / ndays

plot(dat.ind, dat[dat.ind, 2], type='l', xlim = c(0, 800))

plot(dat[,2], type='l', xlim = c(0, 1400))

win = 30
qq = sapply(1:(n - win), function(i) quantile(dat[i:(i+win), 2], 0.95))
xx = (1:(n-win)) + floor(win/2)

plot(dat[,2], type='l', xlim = c(0, 1800))
lines(xx, qq, col = 'red')



X = cbind(1, cos(2*pi/365*xx), sin(2*pi/365*xx))
bhat = solve(t(X) %*% X) %*% t(X) %*% log(qq)
pred.x = 1:n
pred.y = exp(cbind(1, cos(2*pi/365*pred.x), sin(2*pi/365*pred.x)) %*% bhat)
plot(dat[,2], type='l', xlim = c(1500, 3000))
plot(dat[,2], type='l', xlim = c(0, 1500))
lines(xx, qq, col = 'red')
lines(pred.x, pred.y, col = 'darkgreen', lwd = 3)

plot(dat[,2] - pred.y, xlim = c(0, 1500))
new.ind = which(dat[,2] - pred.y > 0)
plot(density(dat[new.ind, 2]))

length(new.ind)

plot(xx, log(qq), xlim = c(0, 1500), type='l')
lines(pred.x, log(pred.y), col = 'darkgreen', lwd = 3)


winter = dat[dat.ind, 2]

#winter = exp(rnorm(4500))
#winter = runif(4500)

threshold = as.vector(quantile(winter, 0.95))
ind = which(winter > threshold)

exceedance = winter[ind] - threshold

y = list("exceedance" = exceedance, "threshold" = threshold,
    "n" = length(winter), "n_u" = length(exceedance),
    "nyears" = nyears, "ndays" = ndays)

calc.post = function(dat, params){
    y = dat$exceed
    u = dat$threshold
    n = dat$n
    n_u = dat$n_u
    ndays = dat$ndays
    nyears = dat$nyears
    mu = params[1]
    sig = params[2]
    ksi = params[3]
    
    if (sig <= 0) return (-Inf)
    if (1 + ksi*(u - mu)/sig <= 0) return (-Inf)
    if (min(1 + ksi*(y - mu + u)/sig) <= 0) return (-Inf)

    if (ksi == 0) ksi = 1e-10

    # likelihood
    out = -nyears*(1 + ksi*(u - mu) / sig)^(-1/ksi) - n_u*log(sig) +
        (-1/ksi - 1)*sum(log(1 + ksi*(y - mu + u)/sig))

    # priors
    out = out + dnorm(ksi, 0, 1, log = TRUE)
    out = out + dnorm(mu, 0, 1000, log = TRUE)
    out = out - log(sig)

    return (out)
    }

out = mcmc_sampler(y, calc.post, 3, nburn = 50000, nmcmc = 250000)

mean(out$accept)
plot(out$params[,3])

mu = out$params[,1]
sig = out$params[,2]
ksi = out$params[,3]
zeta = rbeta(length(mu), y$n_u + 1, y$n - y$n_u + 1)

squig = sig + ksi*(y$threshold - mu)

plot(density(mu))
plot(density(squig))
plot(density(ksi))


param.star = cbind(ksi, squig, zeta)
return.period = sort(unique(c(20, 30, 50, exp(seq(log(0.1), log(100), length = 100)))))

ZZ = apply(param.star, 1,
    function(x)
        y$threshold[1] + x[2]/x[1] * (((return.period*ndays)*x[3])^x[1] - 1)
    )
Zm = apply(ZZ, 1, mean)
Zq = apply(ZZ, 1, quantile, c(0.025, 0.975))

color = "dodgerblue"

par(mfrow = c(1, 1), mar = c(5.1, 4.1, 4.1, 2.1))
xaxis = NULL
for (i in floor(log(0.1, 10)):(ceiling(log(100, 10))-1))
    xaxis = c(xaxis, seq(1, 9, by = 1)*10^i)
xaxis = c(xaxis, 10^ceiling(log(100, 10)))
xlab = xaxis
xlab[-seq(1, length(xaxis), by = 9)] = NA
plot(log(return.period), Zm, ylim = range(Zq),
    type='n', axes = FALSE,
    xlab = paste0("Return Period (# of winters)"),
    ylab = paste0("Return Level (m/days)"),
#   main = "Mean Return Level Plot", col = col_mult(output$base.color, 'gray50'),
        cex.lab = 1.2, cex.main = 1.5)
axis(1, at = log(xaxis), labels = xlab); axis(2)

polygon(c(log(return.period), rev(log(return.period))),
    c(Zq[1,], rev(Zq[2,])), col = col_fade(color, 0.5),
    border = NA)
lines(log(return.period), Zq[1,],
    col = col_fade(color, 0.5), lwd = 2, lend = 1)
lines(log(return.period), Zq[2,],
    col = col_fade(color, 0.5), lwd = 2, lend = 1)

points(log(tail(1/((y$n:1)/(y$n+1)) / ndays, y$n_u)),
    y$threshold + sort(y$exceed), pch = 16,
    col = col_fade("black", 0.7), type='o')




pred = (runif(length(squig))^(-ksi) - 1) * squig / ksi
hist(y$exceed, col = 'gray', breaks = 20, freq = FALSE)
lines(density(pred, n = 10000), col = 'green', lwd = 3)

