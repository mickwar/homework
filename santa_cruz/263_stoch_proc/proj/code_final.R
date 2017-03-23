library(mwBASE)
tmp = read.table("./california_pr_obs19501999.txt", header = TRUE)

rm.leap = grep("-02-29", tmp[,1])
#rm.leap = c(grep("-02-29", tmp[,1]), grep("1982|1983|1998", tmp[,1]))
# Remove those pesky leap years
tt = tmp[-rm.leap, 1]
y = tmp[-rm.leap, 2]

n = length(y)
nyears = 50

#tt = 1:18262
#y = rnorm(length(tt))


# tmp = substr(tt, 5, 10)
# tmp = unique(tmp)[-366]
# qq = double(365)
# for (i in 1:365)
#     qq[i] = quantile(y[grep(tmp[i], tt)], 0.95)
# xx = 1:365
# plot(xx, qq)
    
pdf("./figs/data.pdf", width = 8, height = 8)
plot(head(y, 2000), type='l', ylab = "Total precipitation", axes = FALSE, xlab = "Year",
    main = "California precipitation", cex.main = 2, cex.lab = 1.4)
axis(2)
axis(1, at = 365*(0:7), labels = 1950 + 0:7)
dev.off()

win = 30
qq = sapply(1:(n - win), function(i) quantile(y[i:(i+win)], 0.95))
xx = (1:(n-win)) + floor(win/2)
make.X = function(x, p = 1){
    out = matrix(1, length(x), 1)
    for (i in 1:p){
        out = cbind(out, cos(2*pi/(365/i)*x), sin(2*pi/(365/i)*x))
        }
    return (out)
    }
p = 2
X = make.X(xx, p)
bhat = solve(t(X) %*% X) %*% t(X) %*% log(qq)
pred.x = 1:n
pred.y = exp(make.X(pred.x, p) %*% bhat)


pdf("./figs/threshold.pdf", width = 8, height = 8)
plot(head(y, 2000), type='l', ylab = "Total precipitation", axes = FALSE, xlab = "Year",
    main = "California precipitation", cex.main = 2, cex.lab = 1.4)
axis(2)
axis(1, at = 365*(0:7), labels = 1950 + 0:7)
lines(pred.x, pred.y, lwd = 3, col = 'dodgerblue')
legend("topright", lwd = 3, col = 'dodgerblue', legend = "Threshold", bty = 'n', cex = 2)
dev.off()


## Get exceedances index
ind = which(y - pred.y > 0)
dat = list("y" = y, "threshold" = pred.y, "ind" = ind, "nyears" = nyears,
    "exceed" = (y - pred.y)[ind])

# proportion of exceedances
length(ind) / n

# should be approximately uniform
hist(ind / n, col = 'gray', breaks = 20)

# number of exceedances by day/month
plot(table(substr(tt[ind], 5, 10)) / length(ind), pch = 16, type='l')

# number of exceedances by month
plot(table(substr(tt[ind], 5, 8)) / length(ind), bty = 'n')

# number of exceedances by year
plot(table(substr(tt[ind], 1, 4)) / length(ind), bty = 'n')

pdf("./figs/exceedance_loc.pdf", width = 12, height = 8)
par(mfrow = c(1,2))
plot(table(substr(tt[ind], 5, 8)) / length(ind), bty = 'n', ylab = "", xlab = "Month")
plot(table(substr(tt[ind], 1, 4)) / length(ind), bty = 'n', ylab = "", xlab = "Year")
title("Proportion of exceedance by time", outer = TRUE, line = -2, cex.main = 2)
par(mfrow = c(1,1))
dev.off()


calc.like = function(dat, params){
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
    ksi_vec = ksi_int + ksi_a*cos(2*pi/365 * tt) + ksi_b*sin(2*pi/365 * tt)
#   ksi_vec = rep(ksi_int, length(tt))

    mu_vec = mu_int + mu_a*cos(2*pi/365 * tt) + mu_b*sin(2*pi/365 * tt)
    sig_vec = exp(sig_int + sig_a*cos(2*pi/365 * tt) + sig_b*sin(2*pi/365 * tt))
    
    ksi_vec = ifelse(ksi_vec == 0, 1e-10, ksi_vec)

    A = (1 + ksi_vec*(u_vec - mu_vec)/sig_vec)
    B = (1 + ksi_vec*(y - mu_vec)/sig_vec)

    if (min(A) <= 0) return (-Inf)
    if (min(B) <= 0) return (-Inf)


    # likelihood
    out = sum(-nyears * A^(-1/ksi_vec)-log(sig_vec)-(1/ksi_vec+1)*log(B))

    return (out)
    }

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
    ksi_vec = ksi_int + ksi_a*cos(2*pi/365 * tt) + ksi_b*sin(2*pi/365 * tt)
#   ksi_vec = rep(ksi_int, length(tt))

    mu_vec = mu_int + mu_a*cos(2*pi/365 * tt) + mu_b*sin(2*pi/365 * tt)
    sig_vec = exp(sig_int + sig_a*cos(2*pi/365 * tt) + sig_b*sin(2*pi/365 * tt))
    
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
par(mfrow = c(3,3), mar = c(4.1, 2.1, 2.1,1.1))
for (i in 1:9)
    plot(out$params[,i], type='l', bty = 'n')
par(mfrow = c(1, 1), mar = c(5.1, 4.1, 4.1, 2.1))


### Plots of the time-vary parameters
mu_vec = apply(out$params[,1:3], 1,
    function(p) p[1] + p[2]*cos(2*pi/365*1:365) + p[3]*sin(2*pi/365*1:365))
sig_vec = apply(out$params[,4:6], 1,
    function(p) exp(p[1] + p[2]*cos(2*pi/365*1:365) + p[3]*sin(2*pi/365*1:365)))
ksi_vec = apply(out$params[,7:9], 1,
    function(p) p[1] + p[2]*cos(2*pi/365*1:365) + p[3]*sin(2*pi/365*1:365))
#ksi_vec = matrix(rep(out$params[,7], each = 365), 365)

tvar_mu = rbind(apply(mu_vec, 1, mean), apply(mu_vec, 1, quantile, c(0.025, 0.975)))
tvar_sig = rbind(apply(sig_vec, 1, mean), apply(sig_vec, 1, quantile, c(0.025, 0.975)))
tvar_ksi = rbind(apply(ksi_vec, 1, mean), apply(ksi_vec, 1, quantile, c(0.025, 0.975)))

pdf("./figs/post_mu.pdf", width = 8, height = 8)
plot(1:365, tvar_mu[1,], lwd = 3, ylim = range(tvar_mu), cex.main = 2,
    main = expression(mu), axes = FALSE, xlab = "", ylab = "")
lines(1:365, tvar_mu[2,])
lines(1:365, tvar_mu[3,])
abline(v = seq(1, 365, length = 13), lty = 2, col = 'gray')
axis(2)
axis(1, 15 + seq(1, 335, length = 12), month.abb, lwd = 0)
dev.off()

pdf("./figs/post_sig.pdf", width = 8, height = 8)
plot(1:365, tvar_sig[1,], lwd = 3, ylim = range(tvar_sig), cex.main = 2,
    main = expression(sigma), axes = FALSE, xlab = "", ylab = "")
lines(1:365, tvar_sig[2,])
lines(1:365, tvar_sig[3,])
abline(v = seq(1, 365, length = 13), lty = 2, col = 'gray')
axis(2)
axis(1, 15 + seq(1, 335, length = 12), month.abb, lwd = 0)
dev.off()

pdf("./figs/post_ksi.pdf", width = 8, height = 8)
plot(1:365, tvar_ksi[1,], lwd = 3, ylim = range(tvar_ksi), cex.main = 2,
    main = expression(xi), axes = FALSE, xlab = "", ylab = "")
lines(1:365, tvar_ksi[2,])
lines(1:365, tvar_ksi[3,])
abline(v = seq(1, 365, length = 13), lty = 2, col = 'gray')
axis(2)
axis(1, 15 + seq(1, 335, length = 12), month.abb, lwd = 0)
dev.off()


### Diagnostic plots
k = length(ind)
stand_y = matrix(0, k, nrow(out$params))
for (i in 1:k){
    tmp_mu = out$params[,1] + out$params[,2]*cos(2*pi/365*ind[i]) +
        out$params[,3]*sin(2*pi/365*ind[i])
    tmp_sig = exp(out$params[,4] + out$params[,5]*cos(2*pi/365*ind[i]) +
        out$params[,6]*sin(2*pi/365*ind[i]))
    tmp_ksi = out$params[,7] + out$params[,8]*cos(2*pi/365*ind[i]) +
        out$params[,9]*sin(2*pi/365*ind[i])
#   tmp_ksi = out$params[,7]
    tmp_squig = tmp_sig + tmp_ksi*(pred.y[ind[i]] - tmp_mu)

    stand_y[i,] = 1/tmp_ksi*log(1 + tmp_ksi*(y[ind[i]] - pred.y[ind[i]])/tmp_squig)
    }

y_mm1 = apply(1-exp(-stand_y), 1, mean)
y_qq1 = apply(1-exp(-stand_y), 1, quantile, c(0.025, 0.975))
y_mm2 = apply(stand_y, 1, mean)
y_qq2 = apply(stand_y, 1, quantile, c(0.025, 0.975))

pdf("./figs/diag.pdf", width = 12, height = 8)
par(mfrow = c(1, 2))
# Probability plot
ord = order(y_mm1)
plot((1:k)/(k+1), y_mm1[ord], bty = 'n', main = "Probability plot",
    xlab = "Empirical", ylab = "Model", pch = 16, cex.main = 2, cex.lab = 1.4)
lines((1:k)/(k+1), y_qq1[1, ord])
lines((1:k)/(k+1), y_qq1[2, ord])
abline(0, 1)

# Quantile plot
ord = order(y_mm2)
plot(-log(1-(1:k)/(k+1)), y_mm2[ord], ylim = range(y_qq2), main = "QQ plot",
    bty='n', xlab = "Empirical", ylab = "Model", pch = 16, cex.main = 2, cex.lab = 1.4)
lines(-log(1-(1:k)/(k+1)), y_qq2[1, ord])
lines(-log(1-(1:k)/(k+1)), y_qq2[2, ord])
abline(0, 1)
par(mfrow = c(1, 1))
dev.off()

# DIC
d.theta = apply(out$params, 1, function(p) calc.like(dat, p))
mean(d.theta) + 0.5*var(d.theta)
mean(d.theta)
0.5*var(d.theta)


### Return levels
u.year = exp(make.X(1:365, p) %*% bhat)
#u.year = pred.y[1:365]
zeta = matrix(0, 365, nrow(out$params))
for (i in 1:365)
    zeta[i,] = (1 + ksi_vec[i,]*(u.year[i] - mu_vec[i,])/sig_vec[i,])^(-1/ksi_vec[i,])



get_zm = function(m){
    z = matrix(0, nrow(out$params), 365)
    for (i in 1:365){
#       tmp_squig = sig_vec[i,] + ksi_vec[i,] * (u.year[i] - mu_vec[i,])
#       z[,i] =  u.year[i] + tmp_squig/ksi_vec[i,]*( (m^ksi_vec[i,]) - 1)
#       z[,i] =  u.year[i] + tmp_squig/ksi_vec[i,]*( (zeta[i,]*m)^(ksi_vec[i,]) - 1)
        z[,i] =  mu_vec[i,] + sig_vec[i,]/ksi_vec[i,]*( (m)^(ksi_vec[i,]) - 1)
        }
    return (z)
    }

z.10 = get_zm(10)
z.20 = get_zm(20)
z.50 = get_zm(50)
daily.10 = aggregate(y, list(substr(tt, 5, 10)), function(x) tail(sort(x), 10))

pdf("./figs/return10.pdf", width = 8, height = 8)
z.mm = apply(z.10, 2, mean)
z.qq = apply(z.10, 2, quantile, c(0.025, 0.975))
plot(z.mm, pch = 16, ylim = c(0, max(z.10, y)), axes = FALSE, type='l', lwd = 3,
    xlab = "", ylab = "Total precipitation", main = "10-year return level")
#segments(1:365, z.qq[1,], 1:365, z.qq[2,])
for (i in 1:365)
    points(rep(i, 10), unlist(daily.10[i,][2]), col = rainbow(10), pch = 16)
lines(daily.10$x[,10])
abline(v = seq(1, 365, length = 13), lty = 2, col = 'gray')
axis(2)
axis(1, 15 + seq(1, 335, length = 12), month.abb, lwd = 0)
polygon(c(1:365, 365:1), c(z.qq[1,], rev(z.qq[2,])), col = col_fade("black", 1.0), border = FALSE)
dev.off()

pdf("./figs/return20.pdf", width = 8, height = 8)
z.mm = apply(z.20, 2, mean)
z.qq = apply(z.20, 2, quantile, c(0.025, 0.975))
plot(z.mm, pch = 16, ylim = c(0, max(z.20, y)), axes = FALSE, type='l', lwd = 3,
    xlab = "", ylab = "Total precipitation", main = "20-year return level")
#segments(1:365, z.qq[1,], 1:365, z.qq[2,])
for (i in 1:365)
    points(rep(i, 10), unlist(daily.10[i,][2]), col = rainbow(10), pch = 16)
lines(daily.10$x[,10])
abline(v = seq(1, 365, length = 13), lty = 2, col = 'gray')
axis(2)
axis(1, 15 + seq(1, 335, length = 12), month.abb, lwd = 0)
polygon(c(1:365, 365:1), c(z.qq[1,], rev(z.qq[2,])), col = col_fade("black", 1.0), border = FALSE)
dev.off()

pdf("./figs/return50.pdf", width = 8, height = 8)
z.mm = apply(z.50, 2, mean)
z.qq = apply(z.50, 2, quantile, c(0.025, 0.975))
plot(z.mm, pch = 16, ylim = c(0, max(z.50, y)), axes = FALSE, type='l', lwd = 3,
    xlab = "", ylab = "Total precipitation", main = "50-year return level")
#segments(1:365, z.qq[1,], 1:365, z.qq[2,])
for (i in 1:365)
    points(rep(i, 10), unlist(daily.10[i,][2]), col = rainbow(10), pch = 16)
lines(daily.10$x[,10])
abline(v = seq(1, 365, length = 13), lty = 2, col = 'gray')
axis(2)
axis(1, 15 + seq(1, 335, length = 12), month.abb, lwd = 0)
polygon(c(1:365, 365:1), c(z.qq[1,], rev(z.qq[2,])), col = col_fade("black", 1.0), border = FALSE)
dev.off()

## tmp = sort(z.10[,c(1:59, 335:365)])
## tmp = tail(tmp, nrow(z.10))
## t.mm = mean(tmp)
## t.qq = quantile(tmp, c(0.025, 0.975))
## lines(1:59, rep(t.mm, 59), col = 'blue', lwd = 3)
## lines(335:365, rep(t.mm, 31), col = 'blue', lwd = 3)
## lines(1:59, rep(t.qq[1], 59), col = 'blue')
## lines(335:365, rep(t.qq[1], 31), col = 'blue')
## lines(1:59, rep(t.qq[2], 59), col = 'blue')
## lines(335:365, rep(t.qq[2], 31), col = 'blue')
## 
## tmp = sort(z.10[,152:243])
## tmp = tail(tmp, nrow(z.10))
## t.mm = mean(tmp)
## t.qq = quantile(tmp, c(0.025, 0.975))
## lines(152:243, rep(t.mm, 92), col = 'red', lwd = 3)
## lines(152:243, rep(t.qq[1], 92), col = 'red')
## lines(152:243, rep(t.qq[2], 92), col = 'red')
## 
## tmp = sort(z.10[,60:151])
## tmp = tail(tmp, nrow(z.10))
## t.mm = mean(tmp)
## t.qq = quantile(tmp, c(0.025, 0.975))
## lines(60:151, rep(t.mm, 92), col = 'green', lwd = 3)
## lines(60:151, rep(t.qq[1], 92), col = 'green')
## lines(60:151, rep(t.qq[2], 92), col = 'green')
## 
## tmp = sort(z.10[,244:334])
## tmp = tail(tmp, nrow(z.10))
## t.mm = mean(tmp)
## t.qq = quantile(tmp, c(0.025, 0.975))
## lines(244:334, rep(t.mm, 91), col = 'orange', lwd = 3)
## lines(244:334, rep(t.qq[1], 91), col = 'orange')
## lines(244:334, rep(t.qq[2], 91), col = 'orange')


## # 5-, 10-, 20-, 50-year at the 15th of every month.
## mid = round(seq(15, 350, length = 12))
## zeta = matrix(0, nrow(out$params), 365)
## 
## for (i in 1:365)
##     zeta[,i] = (1 + ksi_vec[i,] * 
##         (pred.y[i] - mu_vec[i,]) / sig_vec[i,]) ^ (-1/ksi_vec[i,])
## 
## range(colMeans(zeta))
## 
## 
## 
## 
## plot(pred.y[1:365])
## 
## 
## mu = out$params[,1]
## sig = out$params[,2]
## ksi = out$params[,3]
## zeta = rbeta(length(mu), y$n_u + 1, y$n - y$n_u + 1)
## 
## squig = sig + ksi*(y$threshold - mu)
## 
## plot(density(mu))
## plot(density(squig))
## plot(density(ksi))
## 
## 
## param.star = cbind(ksi, squig, zeta)
## return.period = sort(unique(c(20, 30, 50, exp(seq(log(0.1), log(100), length = 100)))))
## 
## ZZ = apply(param.star, 1,
##     function(x)
##         y$threshold[1] + x[2]/x[1] * (((return.period*ndays)*x[3])^x[1] - 1)
##     )
## Zm = apply(ZZ, 1, mean)
## Zq = apply(ZZ, 1, quantile, c(0.025, 0.975))
## 
## color = "dodgerblue"
## 
## par(mfrow = c(1, 1), mar = c(5.1, 4.1, 4.1, 2.1))
## xaxis = NULL
## for (i in floor(log(0.1, 10)):(ceiling(log(100, 10))-1))
##     xaxis = c(xaxis, seq(1, 9, by = 1)*10^i)
## xaxis = c(xaxis, 10^ceiling(log(100, 10)))
## xlab = xaxis
## xlab[-seq(1, length(xaxis), by = 9)] = NA
## plot(log(return.period), Zm, ylim = range(Zq),
##     type='n', axes = FALSE,
##     xlab = paste0("Return Period (# of winters)"),
##     ylab = paste0("Return Level (m/days)"),
## #   main = "Mean Return Level Plot", col = col_mult(output$base.color, 'gray50'),
##         cex.lab = 1.2, cex.main = 1.5)
## axis(1, at = log(xaxis), labels = xlab); axis(2)
## 
## polygon(c(log(return.period), rev(log(return.period))),
##     c(Zq[1,], rev(Zq[2,])), col = col_fade(color, 0.5),
##     border = NA)
## lines(log(return.period), Zq[1,],
##     col = col_fade(color, 0.5), lwd = 2, lend = 1)
## lines(log(return.period), Zq[2,],
##     col = col_fade(color, 0.5), lwd = 2, lend = 1)
## 
## points(log(tail(1/((y$n:1)/(y$n+1)) / ndays, y$n_u)),
##     y$threshold + sort(y$exceed), pch = 16,
##     col = col_fade("black", 0.7), type='o')
## 
## 
## 
## 
## pred = (runif(length(squig))^(-ksi) - 1) * squig / ksi
## hist(y$exceed, col = 'gray', breaks = 20, freq = FALSE)
## lines(density(pred, n = 10000), col = 'green', lwd = 3)
## 
