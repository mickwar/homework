library(maps)
library(geoR)
library(MASS)
autotune = function(accept, target = 0.25, k = 2.5)
    (1+(cosh(accept-target)-1)*(k-1)/(cosh(target-
        ceiling(accept-target))-1))^sign(accept-target)
dat = read.table("./northeast_data.txt", header = TRUE)
dat = dat[-31,]

y = log(dat$Ozone)
n = length(y)
lon = dat$Longitude
lat = dat$Latitude
ele = log(dat$Altitude+1)  # Meters above sea level
loc = cbind(lon, lat, ele)

### Plot data points on map
col.vals = (y - min(y)) / diff(range(y))
cols = rgb(col.vals, 0, 1-col.vals)
# pdf("./figs/data.pdf", width = 9, height = 9)
map("state", xlim = range(loc[,1])+c(-1,2), ylim = range(loc[,2])+c(-1, 1),
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
# dev.off()

### Explore trend and anisotropies
 plot(lon, y, bty = 'n', pch = 16, col = cols)
 plot(lat, y, bty = 'n', pch = 16, col = cols)
 plot(ele, y, bty = 'n', pch = 16, col = cols)

# mod = step(lm(y ~ . + .^2 + I(loc^2), data = data.frame(loc)), scope = list("lower" = lm(y ~ 1),
#     "upper" = lm(y ~ . + .^2 + I(loc^2), data = data.frame(loc))),
#     direction = "both", k = log(n))
#mod2 = lm(y ~ lon + lat + ele + I(lon^2) + I(lat^2) + I(ele^2) + lat*ele, data = data.frame(loc))

### Based on previous analysis
# mod = lm(y ~ 1 + lat + ele + I(lat^2) + lat*ele, data = data.frame(loc))
# mod = lm(y ~ 1 + lat + I(lat^2) + I(ele^2), data = data.frame(loc))
mod = lm(y ~ 1 + lon + lat + lon*lat, data = data.frame(loc))
summary(mod)

yhat = predict(mod)


# pdf("./figs/regression.pdf", height = 15, width = 9)
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
# dev.off()

# Residuals
# res = y - yhat


# pdf("./figs/variogram.pdf", height = 6, width = 12)
# par(mfrow = c(1, 2))
# tmp = variog(coords = loc[,-3], data = res, trend = y ~ 1 + loc[,3])
# plot(tmp$u, tmp$v, pch = 16, ylim = c(0, max(tmp$v)), bty = 'n',
#     main = "Binned Semivariogram for log Ozone", cex.main = 1.5,
#     xlab = "Distance", ylab = "Semivariance")
# tmp = variog4(coords = loc[,-3], data = res, trend = y ~ 1 + loc[,3])
# tmp = variog4(coords = loc[,-3], data = res, uvec = 13)
# tmp.xlim = range(sapply(tmp[-5], function(x) range(x$u)))
# tmp.ylim = c(0, max(sapply(tmp[-5], function(x) max(x$v, na.rm = TRUE))))
# plot(0, type='n', bty='n', xlim = tmp.xlim, ylim = tmp.ylim,
#     main = "Directional Semivariogram for log Ozone", cex.main = 1.5,
#     xlab = "Distance", ylab = "Semivariance")
# for (i in 1:4)
#     lines(tmp[[2]]$u, tmp[[i]]$v, col = i, lty = i)
# legend(tmp.xlim[1], tmp.ylim[2], legend = expression(0^o, 45^o, 90^o, 135^o),
#     lty = 1:4, col = 1:4)
# par(mfrow = c(1, 1))
# # dev.off()
# 
# # ### LSE fit to Matern correlation
# obs = variog(coords = loc[,-3], data = res, trend = y ~ 1 + loc[,3])
# 
# # C(u) = sig^2 * rho(u) + kappa * (u == 0)              covariance
# # gamma(u) = sig^2 * (1 - rho(u)) + kappa * (u > 0)     semi-variogram
# semiv = function(u, phi, nu, sig2, kappa)
#     sig2 * (1-matern(u, phi, nu)) + kappa
# make.f = function(nu){
#     f.opt = function(p){
#         phi = p[1]
#         sig2 = p[2]
#         kappa = p[3]
#         sum( (semiv(obs$u, phi, nu, sig2, kappa) - obs$v)^2)
#         }
#     return (f.opt)
#     }
# 
# # pdf("./figs/matern.pdf", width = 12, height = 12)
# par(mfrow = c(2, 2))
# for (nu in c(0.5, 1, 1.5, 2.5)){
#     pars = optim(rep(1, 3), make.f(nu), lower = c(0, 0, 0), method = "L-BFGS-B")$par
#     plot(obs$u, obs$v, ylim = c(0, max(obs$v)), pch = 16, bty = 'n')
#     lines(obs$u, semiv(obs$u, pars[1], nu, pars[2], pars[3]), col = 'red')
#     text(rep(11, 4), seq(0.0050, 0.0035, length = 4), c(round(pars, 3), nu),
#         pos = 4, cex = 1.5)
#     text(rep(9.5, 4), seq(0.0050, 0.0035, length = 4), expression(sigma^2, phi, kappa, nu),
#         pos = 4, cex = 1.5)
#     }
# par(mfrow = c(1, 1))
# # dev.off()



### MCMC
# D = cbind(1, loc[,2], loc[,3], loc[,2]^2, loc[,2] * loc[,3])
#D = cbind(1, lon, lat, lon*lat, lon^2, lat^2, ele)
D = cbind(1, lon, lat, lon*lat)
X = y
k = NCOL(D)
dists = as.matrix(dist(loc[,1:2], upper = TRUE, diag= TRUE))
make.K = function(psi, gamma2, nu)
    1/gamma2*matern(dists, psi, nu) + diag(n)

# For tau^2 (observational error)
prior.a = 3
prior.b = 1

# For psi (range)
prior.c = 10
prior.d = 0.1

# For gamma^2 (= tau^2 / sigma^2, obs variance / spatial variance)
prior.e = 3
prior.f = 1

nburn = 20000
nmcmc = 30000
window = 500
param.beta = matrix(0, nburn + nmcmc, k)
param.tau2 = double(nburn + nmcmc)
param.psi = double(nburn + nmcmc)
param.gamma2 = double(nburn + nmcmc)

param.tau2[1] = 1
param.psi[1] = 1
param.gamma2[1] = 1
param.nu[1] = nu.states[1]
cand.sig = 1e-3 * diag(2)
accept = double(nburn + nmcmc)
vec.like = double(nburn + nmcmc)


# system.time({
# for (j in 1:100){
#     K = make.K(param.psi[1], param.gamma2[1], nu)
#     U = chol(K)
#     U.inv = backsolve(U, diag(n))
#     DU = t(D) %*% U.inv
#     CDU = tcrossprod(DU)
#     beta.hat = solve(CDU) %*% (DU %*% (t(U.inv) %*% X))
#     S2 = sum((t(X - D %*% beta.hat) %*% U.inv)^2)
#     Det1 = -sum(log(diag(U)))
#     Det2 = -0.5*determinant(CDU)$modulus[1]
#     post.curr = Det1 + Det2 + (-(n-k)/2-prior.a) + log(S2 + 2*prior.b) +
#         dgamma(param.psi[1], prior.c, prior.d, log = TRUE) +
#         dgamma(param.gamma2[1], prior.e, prior.f, log = TRUE)
#     }
#     })

nu = 1.5    # Based on DIC
K = make.K(param.psi[1], param.gamma2[1], param.nu[1])
K.inv = solve(K)
beta.hat = solve(t(D) %*% K.inv %*% D) %*% t(D) %*% K.inv %*% X
S2 = t(X - D %*% beta.hat) %*% K.inv %*% (X - D %*% beta.hat)
Det1 = -0.5*determinant(K)$modulus[1] 
Det2 = -0.5*determinant(t(D) %*% K.inv %*% D)$modulus[1]
post.curr = Det1 + Det2 + (-(n-k)/2-prior.a) + log(S2 + 2*prior.b) +
    dgamma(param.psi[1], prior.c, prior.d, log = TRUE) +
    dgamma(param.gamma2[1], prior.e, prior.f, log = TRUE) +
    prior.nu(param.nu[1])
    
vec.like[1] = Det1 - n/2*log(param.tau2[1]) -
    1/(2*param.tau2[1])*t(X - D %*% param.beta[1,]) %*% K.inv %*% (X - D %*% param.beta[1,])



for (i in 2:(nburn + nmcmc)){
    cat(i, "/", nburn + nmcmc, "\r")
    param.beta[i,] = param.beta[i-1,]
    param.tau2[i] = param.tau2[i-1]
    param.psi[i] = param.psi[i-1]
    param.gamma2[i] = param.gamma2[i-1]
    param.nu[i] = param.nu[i-1]
    vec.like[i] = vec.like[i-1]

    cand = mvrnorm(1, c(param.psi[i-1], param.gamma2[i-1]), cand.sig)
    if (all(cand > 0)){
        cand.K = make.K(cand[1], cand[2], nu)
        cand.K.inv = chol2inv(chol(cand.K))
        cand.DKD = t(D) %*% cand.K.inv %*% D
        
        cand.beta.hat = solve(cand.DKD) %*% t(D) %*% cand.K.inv %*% X
        cand.S2 = t(X - D %*% cand.beta.hat) %*% cand.K.inv %*% (X - D %*% cand.beta.hat)
        cand.Det1 = -0.5*determinant(cand.K)$modulus[1]
        cand.Det2 = -0.5*determinant(cand.DKD)$modulus[1]

        post.cand = cand.Det1 + cand.Det2 + (-(n-k)/2-prior.a) + log(cand.S2 + 2*prior.b) + 
            dgamma(cand[1], prior.c, prior.d, log = TRUE) +
            dgamma(cand[2], prior.e, prior.f, log = TRUE)
#           prior.nu(cand[3])


        if (log(runif(1)) <= post.cand - post.curr){
            param.psi[i] = cand[1]
            param.gamma2[i] = cand[2]
            param.nu[i] = cand[3]
            accept[i] = 1
            K = cand.K
            K.inv = cand.K.inv
            beta.hat = cand.beta.hat
            S2 = cand.S2
            Det1 = cand.Det1
            Det2 = cand.Det2
            DKD = cand.DKD
            post.curr = post.cand

            param.tau2[i] = 1/rgamma(1, prior.a + (n-k)/2, prior.b + S2/2)
            param.beta[i,] = mvrnorm(1, beta.hat, param.tau2[i] * solve(DKD))

            vec.like[i] = Det1 - n/2*log(param.tau2[i]) -
                1/(2*param.tau2[i])*t(X - D %*% param.beta[i,]) %*% K.inv %*% (X - D %*% param.beta[i,])
            }
        }
    if ((floor(i/window) == i/window) && (i <= nburn))
        cand.sig = autotune(mean(accept[(i-window+1):i]), target = 0.234, k = window / 50) *
            (cand.sig + window * var(cbind(param.psi, param.gamma2)[(i-window+1):i,]) / i)

    }

param.beta = tail(param.beta, nmcmc)
param.tau2 = tail(param.tau2, nmcmc)
param.psi = tail(param.psi, nmcmc)
param.gamma2 = tail(param.gamma2, nmcmc)
param.nu = tail(param.nu, nmcmc)
vec.like = tail(vec.like, nmcmc)
accept = tail(accept, nmcmc)
sig2 = param.tau2 / param.gamma2


mean(accept)
par(mfrow = c(1,1))

plot(vec.like, type = 'l')
plot(param.psi, type = 'l')
plot(param.tau2, type = 'l')
plot(param.gamma2, type = 'l')
plot(sig2, type = 'l')

plot(param.nu, type = 'l')
abline(h = nu.states, lty = 2)

plot(nu.states, table(param.nu) / nmcmc, ylim = c(0,1))
points(nu.states, exp(prior.nu(nu.states)), col = 'blue')

plot(param.beta[,c(1,2)], pch = 16)

colMeans(param.beta)
coef(mod)

plot(density(param.psi))
plot(density(param.tau2))
plot(density(param.gamma2))
curve(dgamma(x, prior.e, prior.f), add = TRUE, col = 'blue')
plot(density(sig2))

mean(param.tau2)
var(param.tau2)

mean(sig2)
var(sig2)
summary(mod)

mean(param.psi)


### Predictions
#post.betahat = colMeans(param.beta)
pred.x = D %*% colMeans(param.beta)
post.K = make.K(mean(param.psi), mean(param.gamma2), nu)
mse = t(X - pred.x) %*% post.K %*% (X - pred.x)

plot(X, pred.x, pch = 16, bty = 'n')
abline(0, 1)

### DIC
D.theta = -2*vec.like
var(D.theta)/2 + mean(D.theta)
var(D.theta)/2
mean(D.theta)
var(D.theta)/2 + mean(D.theta)
