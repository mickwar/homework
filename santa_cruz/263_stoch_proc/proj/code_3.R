library(mwBASE)
tmp = read.table("./california_pr_obs19501999.txt", header = TRUE)

tt = tmp[,1]
y = tmp[,2]
n = length(y)
nyears = 50

plot(y, type='l', xlim = c(0, 1400))

win = 30
qq = sapply(1:(n - win), function(i) quantile(y[i:(i+win)], 0.95))
xx = (1:(n-win)) + floor(win/2)

plot(y, type='l', xlim = c(0, 1800))
lines(xx, qq, col = 'red')

X = cbind(1, cos(2*pi/365.25*xx), sin(2*pi/365.25*xx))
bhat = solve(t(X) %*% X) %*% t(X) %*% log(qq)
pred.x = 1:n
pred.y = exp(cbind(1, cos(2*pi/365.25*pred.x), sin(2*pi/365.25*pred.x)) %*% bhat)
plot(y, type='l', xlim = c(1500, 3000))
lines(xx, qq, col = 'red')
lines(pred.x, pred.y, col = 'darkgreen', lwd = 3)

plot(y - pred.y, xlim = c(0, 1500))
ind = which(y - pred.y > 0)
plot(density(y[ind]))

dat = list("y" = y, "threshold" = pred.y, "ind" = ind, "nyears" = nyears,
    "exceed" = (y - pred.y)[ind])



calc.post = function(dat, params){
    y = dat$y[dat$ind]
    u_vec = dat$threshold[dat$ind]
    n = length(dat$y)
    tt = (1:n)[dat$ind]
    nyears = dat$nyears

    # mu
    mu_int = params[1]
    mu_a = params[2]
    mu_b = params[3]

    # sigma
    sig_int = params[4]
    sig_a = params[5]
    sig_b = params[6]

    # ksi
    ksi_int = params[7]
    ksi_a = params[8]
    ksi_b = params[9]

    mu_vec = mu_int + mu_a*cos(2*pi/365.25 * tt) + mu_b*sin(2*pi/365.25 * tt)
    sig_vec = exp(sig_int + sig_a*cos(2*pi/365.25 * tt) + sig_b*sin(2*pi/365.25 * tt))
    ksi_vec = ksi_int + ksi_a*cos(2*pi/365.25 * tt) + ksi_b*sin(2*pi/365.25 * tt)
    
    ksi_vec = ifelse(ksi_vec == 0, 1e-10, ksi_vec)

    A = (1 + ksi_vec*(u_vec - mu_vec)/sig_vec)
    B = (1 + ksi_vec*(y - mu_vec)/sig_vec)

    if (min(A) <= 0) return (-Inf)
    if (min(B) <= 0) return (-Inf)


    # likelihood
    out = sum(-nyears * A^(-1/ksi_vec)-log(sig_vec)-(1/ksi_vec+1)*log(B))

    # priors
    out = out + dnorm(mu_int, 0, 100, log = TRUE)
    out = out + dnorm(mu_a, 0, 100, log = TRUE)
    out = out + dnorm(mu_b, 0, 100, log = TRUE)

    out = out + dnorm(sig_int, 0, 100, log = TRUE)
    out = out + dnorm(sig_a, 0, 100, log = TRUE)
    out = out + dnorm(sig_b, 0, 100, log = TRUE)

    out = out + dnorm(ksi_int, 0, 100, log = TRUE)
    out = out + dnorm(ksi_a, 0, 100, log = TRUE)
    out = out + dnorm(ksi_b, 0, 100, log = TRUE)

    return (out)
    }

out = mcmc_sampler(dat, calc.post, 9, nburn = 50000, nmcmc = 50000,
    chain_init = double(9))

mean(out$accept)
plot(out$params[,6], type='l')

mu_vec = apply(out$params[,1:3], 1,
    function(p) p[1] + p[2]*cos(2*pi/365*1:365) + p[3]*sin(2*pi/365*1:365))
sig_vec = apply(out$params[,4:6], 1,
    function(p) exp(p[1] + p[2]*cos(2*pi/365*1:365) + p[3]*sin(2*pi/365*1:365)))
ksi_vec = apply(out$params[,7:9], 1,
    function(p) p[1] + p[2]*cos(2*pi/365*1:365) + p[3]*sin(2*pi/365*1:365))

plot(1:365, apply(mu_vec, 1, mean))
plot(1:365, apply(sig_vec, 1, mean))
plot(1:365, apply(ksi_vec, 1, mean))

plot(pred.y[1:365])


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

