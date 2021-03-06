library(mwBASE)
dat = read.table("./california_pr_obs19501999.txt", header = TRUE)

# Get December, January, and February (but not leap day)
dat.ind = grep("(-12-)|(-01-)|(-02-(?!29))", dat[,1], perl = TRUE)

ndays = 90  # Number of winter days
nyears = length(dat.ind) / ndays

plot(dat.ind, dat[dat.ind, 2], type='l', xlim = c(0, 800))

winter = dat[dat.ind, 2]

threshold = as.vector(quantile(winter, 0.90))
ind = which(winter > threshold)

exceedance = winter[ind] - threshold

y = list("exceedance" = exceedance, "threshold" = threshold,
    "n" = length(winter), "n_u" = length(exceedance),
    "nyears" = nyears, "ndays" = ndays)

calc.post = function(dat, params){
    y = dat$exceed
    n = dat$n
    n_u = dat$n_u
    ndays = dat$ndays
    nyears = dat$nyears
    sig = params[1]
    ksi = params[2]
    
    if (sig <= 0) return (-Inf)
    if (min(1 + ksi*y/sig) <= 0) return (-Inf)

    if (ksi == 0) ksi = 1e-10

    # likelihood
    out = -n_u*log(sig) + (-1/ksi - 1)*sum(log(1 + ksi*y/sig))

    # priors
    out = out + dnorm(ksi, 0, 1, log = TRUE)
    out = out - log(sig)

    return (out)
    }

out = mcmc_sampler(y, calc.post, 2, nburn = 50000, nmcmc = 50000)

mean(out$accept)
plot(out$params[,1])

sig = out$params[,1]
ksi = out$params[,2]
zeta = rbeta(length(sig), y$n_u + 1, y$n - y$n_u + 1)

plot(density(sig))
plot(density(ksi))


param.star = cbind(ksi, sig, zeta)
return.period = sort(unique(c(20, 30, 50, exp(seq(log(0.1), log(100), length = 100)))))

ZZ = apply(param.star, 1,
    function(x)
        y$threshold[1] + x[2]/x[1] * (((return.period*ndays)*x[3])^x[1] - 1)
    )
Zm = apply(ZZ, 1, mean)
Zq = apply(ZZ, 1, quantile, c(0.025, 0.975))

color = "black"

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
    col = col_fade("green", 0.5), type='o')

