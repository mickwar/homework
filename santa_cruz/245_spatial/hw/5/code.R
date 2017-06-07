library(rgl)
library(maps)
library(geoR)
library(spBayes)
obs = read.table("./obs.txt", header = TRUE)
knots = read.table("./knots.txt", header = TRUE)
krig = read.table("./krig.txt", header = TRUE)



# Remove outliers and leverage points(make them NA)
obs$ozone[31] = NA
obs$ozone[72] = NA
obs$ozone[77] = NA
obs$ozone[103] = NA

# Scaled
sc.X = cbind(scale(obs$kmEast), scale(obs$kmNorth),
    scale(obs$alt))
sc.knots = cbind(scale(knots$kmEast), scale(knots$kmNorth))






ind = which(!is.na(obs$ozone))
y = obs$ozone[ind]
X = sc.X[ind,]
X = cbind(X, X[,1]*X[,2], X[,1]*X[,3], X[,2]*X[,3], X[,1]*X[,2]*X[,3],
    X[,1]^2, X[,2]^2, X[,3]^2)
X = data.frame(X)
mod = step(lm(y ~ 1, data = X), scope = list("lower" = lm(y ~ 1, data = X),
    "upper" = lm(y ~ ., data = X)), direction = "both")
summary(mod)
# Ozone: kmNorth, kmNorth^2, (kmEast * kmNorth * alt)

ind = which(!is.na(obs$no2))
y = obs$no2[ind]
X = sc.X[ind,]
X = cbind(X, X[,1]*X[,2], X[,1]*X[,3], X[,2]*X[,3], X[,1]*X[,2]*X[,3],
    X[,1]^2, X[,2]^2, X[,3]^2)
X = data.frame(X)
mod = step(lm(y ~ 1, data = X), scope = list("lower" = lm(y ~ 1, data = X),
    "upper" = lm(y ~ ., data = X)), direction = "both")
summary(mod)
# NO2: kmNorth, alt, kmEast* kmNorth * alt

ind = which(!is.na(obs$pm25))
y = obs$pm25[ind]
X = sc.X[ind,]
X = cbind(X, X[,1]*X[,2], X[,1]*X[,3], X[,2]*X[,3], X[,1]*X[,2]*X[,3],
    X[,1]^2, X[,2]^2, X[,3]^2)
X = data.frame(X)
mod = step(lm(y ~ 1, data = X), scope = list("lower" = lm(y ~ 1, data = X),
    "upper" = lm(y ~ ., data = X)), direction = "both")
summary(mod)
# PM2.5: kmEast, kmNorth, kmEast, (kmNorth^2)



### Ozone predictive process
ind = which(!is.na(obs$ozone))
y = obs$ozone[ind]
X = sc.X[ind,]
X = cbind(X[,2], X[,2]^2, X[,1]*X[,2]*X[,3])

#mod = spLM(obs$ozone[ind] ~ 1 + X, coords = as.matrix(obs[ind, c(4, 5)]),
coord = sc.X[ind,1:2]
mod = spLM(obs$ozone[ind] ~ 1+X, coords = coord, knots = sc.knots,
    n.samples = 10000,
    cov.model = "matern",
    priors = list("beta.flat", "sigma.sq.ig" = c(1, 0.5), "tau.sq.ig" = c(1, 0.01),
        "phi.Unif" = c(1e-6, 100), "nu.Unif" = c(0.5, 3.5)),
    starting = list("beta" = coef(lm(y ~ X)), "sigma.sq" = 6, "tau.sq" = 0.03,
        "phi" = 7, "nu" = 1.3),
    tuning = list("sigma.sq" = 0.1, "tau.sq" = 0.1, "phi" = 0.1, "nu" = 0.1))


str(mod)
apply(mod$acceptance, 1, mean)

colMeans(mod$p.theta.samples)
var(mod$p.theta.samples)

apply(mod$p.theta.samples, 2, function(x) length(unique(x)) / length(x))
var(mod$p.theta.samples)

ss = seq(0, 1, length = nrow(mod$p.theta.samples))
pairs(as.matrix(mod$p.theta.samples), col = rgb(ss, 0, 0))

#pred.X = cbind(obs$kmEast, obs$kmNorth, obs$alt)
#pred.X = cbind(1, X[,2], X[,1]*X[,3], X[,2]*X[,3], X[,2]^2)
pred.X = sc.X
pred.X = cbind(1, pred.X[,2], pred.X[,2]^2, pred.X[,1]*pred.X[,2]*pred.X[,3])

pred.samples = spPredict(mod, pred.coords = cbind(scale(obs$kmEast), scale(obs$kmNorth)), pred.covars = pred.X)
pred.ozone = apply(pred.samples$p.y.predictive.samples, 1, mean)

pred.samples = spPredict(mod, pred.coords = obs[,4:5], pred.covars = matrix(1, nrow(obs), 1))
pred.ozone = apply(pred.samples$p.y.predictive.samples, 1, mean)
plot(density(pred.ozone[5,]))


lonlat = obs[, c(1, 2)]
col.obs = (y - min(y)) / diff(range(y))
col.pred = (pred.ozone - min(pred.ozone)) / diff(range(pred.ozone))

col.pred = rep(3, length(pred.ozone))
col.pred[-ind] = 2

plot3d(lonlat[,1], lonlat[,2], pred.ozone, col = col.pred, type = 'n')
segments3d(rep(lonlat[,1], each = 2), rep(lonlat[,2], each = 2), c(t(cbind(40, pred.ozone))),
    col = rep(col.pred, each = 2))
segments3d(rep(lonlat[ind,1]+0.03, each = 2), rep(lonlat[ind,2]+0.03, each = 2), c(t(cbind(40, y))),
    col = 'blue')
points3d(knots[,1], knots[,2], rep(40, nrow(knots)), col = 'orange')

str(mod)
dim(mod$p.theta.samples)
plot(mod$p.theta.samples[,1])
colMeans(tail(mod$p.theta.samples, 5000))


betas = spRecover(mod, get.w = FALSE)

coef(lm(y ~ X))
apply(as.matrix(betas$p.beta.recover.samples), 2, mean)
apply(as.matrix(betas$p.beta.recover.samples), 2, sd)
