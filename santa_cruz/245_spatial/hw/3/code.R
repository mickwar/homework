library(maps)
library(geoR)
dat = read.table("./northeast_data.txt", header = TRUE)

y = log(dat$Ozone)
lon = dat$Longitude
lat = dat$Latitude
ele = dat$Altitude  # Meters above sea level
loc = cbind(lon, lat, ele)



### Plot data points on map
col.vals = (y - min(y)) / diff(range(y))
cols = rgb(col.vals, 0, 1-col.vals)
pdf("./figs/data.pdf", width = 9, height = 9)
map("state", xlim = range(loc[,1])+c(-1,2), ylim = range(loc[,2]+c(-1, 1)),
    lwd = 2)
points(loc, col = cols, pch = 16, cex = 1.5)
points(x = rep(max(loc[,1]) + 0.67, 100),
    y = seq(min(loc[,2]), max(loc[,2]), length = 100),
    pch = 15, col = rgb(seq(0, 1, length = 100), 0, seq(1, 0, length = 100)),
    cex = 2.5)
lines(rep(max(loc[,1]) + 0.97, 2), range(loc[,2]))
for (i in 1:6){
    lines(max(loc[,1]) + c(0.97, 1.07), rep(min(loc[,2])+(i-1)*diff(range(loc[,2])/5), 2))
#   text(max(loc[,1])+1.15, min(loc[,2])+(i-1)*diff(range(loc[,2])/5),
#       round(quantile(y, (i-1)/5), 3), pos=4)
    text(max(loc[,1])+1.12, min(loc[,2])+(i-1)*diff(range(loc[,2])/5),
        round(seq(min(y), max(y), length = 6)[i], 3), pos=4)
    }
title(main = "Northeast Stations and log Ozone Measurements", cex.main = 1.5)
dev.off()

### Explore trend and anisotropies
# plot(lon, y, bty = 'n', pch = 16, col = cols)
# plot(lat, y, bty = 'n', pch = 16, col = cols)
# plot(ele, y, bty = 'n', pch = 16, col = cols)

mod = step(lm(y ~ . + .^2 + I(loc^2), data = data.frame(loc)), scope = list("lower" = lm(y ~ 1),
    "upper" = lm(y ~ . + .^2 + I(loc^2), data = data.frame(loc))),
    direction = "both")
summary(mod)

yhat = predict(mod)

pdf("./figs/regression.pdf", height = 15, width = 9)
par(mfrow = c(3, 1))
plot(lon, y, type = 'n', bty = 'n', main = "Longitude", cex.main = 1.5,
    xlab = "Longitude", ylab = "log Ozone")
segments(x0 = lon, y0 = y, y1 = yhat, col = "gray80")
points(lon, y, pch = 16, col = "gray70")
points(lon, yhat, pch = 16, col = "black")

plot(lat, y, type = 'n', bty = 'n', main = "Latitude", cex.main = 1.5,
    xlab = "Latitude", ylab = "log Ozone")
segments(x0 = lat, y0 = y, y1 = yhat, col = "gray80")
points(lat, y, pch = 16, col = "gray70")
points(lat, yhat, pch = 16, col = "black")

plot(ele, y, type = 'n', bty = 'n', main = "Altitude", cex.main = 1.5,
    xlab = "Altitude", ylab = "log Ozone")
segments(x0 = ele, y0 = y, y1 = yhat, col = "gray80")
points(ele, y, pch = 16, col = "gray70")
points(ele, yhat, pch = 16, col = "black")
par(mfrow = c(1, 1))
dev.off()

# Residuals
res = y - yhat


pdf("./figs/variogram.pdf", height = 6, width = 12)
par(mfrow = c(1, 2))
tmp = variog(coords = loc[,-3], data = res, trend = y ~ 1 + loc[,3])
plot(tmp$u, tmp$v, pch = 16, ylim = c(0, max(tmp$v)), bty = 'n',
    main = "Binned Semivariogram for log Ozone", cex.main = 1.5,
    xlab = "Distance", ylab = "Semivariance")
tmp = variog4(coords = loc[,-3], data = res, trend = y ~ 1 + loc[,3])
tmp = variog4(coords = loc[,-3], data = res, uvec = 13)
tmp.xlim = range(sapply(tmp[-5], function(x) range(x$u)))
tmp.ylim = c(0, max(sapply(tmp[-5], function(x) max(x$v, na.rm = TRUE))))
plot(0, type='n', bty='n', xlim = tmp.xlim, ylim = tmp.ylim,
    main = "Directional Semivariogram for log Ozone", cex.main = 1.5,
    xlab = "Distance", ylab = "Semivariance")
for (i in 1:4)
    lines(tmp[[2]]$u, tmp[[i]]$v, col = i, lty = i)
legend(tmp.xlim[1], tmp.ylim[2], legend = expression(0^o, 45^o, 90^o, 135^o),
    lty = 1:4, col = 1:4)
par(mfrow = c(1, 1))
dev.off()

### LSE fit to Matern correlation
obs = variog(coords = loc[,-3], data = res, trend = y ~ 1 + loc[,3])

# C(u) = sig^2 * rho(u) + kappa * (u == 0)              covariance
# gamma(u) = sig^2 * (1 - rho(u)) + kappa * (u > 0)     semi-variogram
semiv = function(u, phi, nu, sig2, kappa)
    sig2 * (1-matern(u, phi, nu)) + kappa
make.f = function(nu){
    f.opt = function(p){
        phi = p[1]
        sig2 = p[2]
        kappa = p[3]
        sum( (semiv(obs$u, phi, nu, sig2, kappa) - obs$v)^2)
        }
    return (f.opt)
    }

pdf("./figs/matern.pdf", width = 12, height = 12)
par(mfrow = c(2, 2))
for (nu in c(0.5, 1, 1.5, 2.5)){
    pars = optim(rep(1, 3), make.f(nu), lower = c(0, 0, 0), method = "L-BFGS-B")$par
    plot(obs$u, obs$v, ylim = c(0, max(obs$v)), pch = 16, bty = 'n')
    lines(obs$u, semiv(obs$u, pars[1], nu, pars[2], pars[3]), col = 'red')
    text(rep(1.25, 4), seq(0.0080, 0.0065, length = 4), c(round(pars, 3), nu),
        pos = 4, cex = 1.5)
    text(rep(0.5, 4), seq(0.0080, 0.0065, length = 4), expression(sigma^2, phi, kappa, nu),
        pos = 4, cex = 1.5)
    }
par(mfrow = c(1, 1))
dev.off()

