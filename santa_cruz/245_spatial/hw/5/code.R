library(rgl)
library(maps)
library(geoR)
library(spBayes)
library(fields)
library(akima)
fix.poly = function(z){
    z = cbind(z$x, z$y)

    v = t(cbind(c(1, which(is.na(z[,1]))+1),
        c(which(is.na(z[,1]))-1, nrow(z))))

    island = NULL
    out = z[v[1,1]:v[2,1],]
    prev.ind = 2
    tmp = as.matrix(v[,-1])
    while (NCOL(tmp) > 0){
        ind = which(sapply(c(tmp), function(x) all(z[v[prev.ind],] == z[x,])))
        if (length(ind) == 0){
            island[[length(island)+1]] = unique(out)
            ind = c(2, 1)
            out = z[tmp[ind[1], ind[2]]:tmp[3-ind[1], ind[2]],]
            prev.ind = which(v == tmp[3-ind[1], ind[2]])
            tmp = as.matrix(tmp[,-ind[2]])
        } else {
            ind = which(tmp == tmp[ind], arr.ind = TRUE)
            out = rbind(out, z[tmp[ind[1], ind[2]]:tmp[3-ind[1], ind[2]],])
            prev.ind = which(v == tmp[3-ind[1], ind[2]])
            tmp = as.matrix(tmp[,-ind[2]])
            }
        }
    if (length(island) > 0){
        island[[length(island)+1]] = unique(out)
        out = island
    } else {
        out = unique(out)
        }
        
    return (out)
    }
states = c("Maine", "Vermont", "New Hampshire", "Connecticut",
    "Massachusetts", "New Jersey", "Delaware", "Rhode Island",
    "Maryland", "District Of Columbia")
plot.nice = function(xy, z, zscale, nlevels = 20, ...){
    df = data.frame("x" = xy[,1], "y" = xy[,2], "z" = z)
    fld = with(df, interp(x = x, y = y, z = z))
    xlim = c(min(xy[,1]), max(xy[,1])+diff(range(xy[,1]))*0.2)
    ylim = range(xy[,2])
    offx = diff(range(xy[,1]))*0.2
    if (missing(zscale))
        zscale = range(fld$z, na.rm = TRUE)
    zlevels = seq(zscale[1], zscale[2], length = nlevels)
    plot(0, type='n', xlim = c(min(xy[,1]), max(xy[,1])+offx),
        ylim = range(xy[,2]),
        xlab = "Longitude", ylab = "Latitude", bty = 'n')
    .filled.contour(x = fld$x, y = fld$y, z = fld$z,
        levels = zlevels,
        col = tim.colors(nlevels))
    rect(rep(xlim[2] - offx*0.50 , nlevels), head(seq(ylim[1], ylim[2], length = nlevels+1), nlevels),
        rep(xlim[2] - offx*0.75, nlevels), tail(seq(ylim[1], ylim[2], length = nlevels+1), nlevels),
        col = tim.colors(nlevels))
    text(x = rep(xlim[2] - offx*0.25, 7), y = seq(ylim[1], ylim[2], length = 5),
        labels = round(seq(zscale[1], zscale[2], length = 5), 2))
    title(...)
    map("state", add = TRUE)
    }


obs = read.table("./obs.txt", header = TRUE)
knots = read.table("./knots.txt", header = TRUE)
krig = read.table("./krig.txt", header = TRUE)



# Remove outliers and leverage points (make them NA)
obs$ozone[31] = NA
obs$ozone[72] = NA
obs$ozone[77] = NA
obs$ozone[103] = NA
obs$no2[104] = NA

# Remove rows that now have all NA in ozone, no2, pm25
rm = apply(obs[,6:8], 1, function(x) all(is.na(x)))
obs = obs[!rm,]

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
# NO2: kmNorth, alt, kmEast* kmNorth * alt, kmEast^2

ind = which(!is.na(obs$pm25))
y = obs$pm25[ind]
X = sc.X[ind,]
X = cbind(X, X[,1]*X[,2], X[,1]*X[,3], X[,2]*X[,3], X[,1]*X[,2]*X[,3],
    X[,1]^2, X[,2]^2, X[,3]^2)
X = data.frame(X)
mod = step(lm(y ~ 1, data = X), scope = list("lower" = lm(y ~ 1, data = X),
    "upper" = lm(y ~ ., data = X)), direction = "both")
summary(mod)
# PM2.5: kmEast, kmNorth, alt, (kmNorth^2)


### Ozone predictive process
ind = which(!is.na(obs$ozone))
y = obs$ozone[ind]
X = sc.X[ind,]
X = cbind(X[,2], X[,2]^2, X[,1]*X[,2]*X[,3])

#mod = spLM(obs$ozone[ind] ~ 1 + X, coords = as.matrix(obs[ind, c(4, 5)]),
coord = sc.X[ind,1:2]
ozone.pp = spLM(obs$ozone[ind] ~ 1+X, coords = coord, knots = sc.knots,
    n.samples = 10000,
    cov.model = "matern",
    priors = list("beta.flat", "sigma.sq.ig" = c(1, 0.5), "tau.sq.ig" = c(1, 0.01),
        "phi.Unif" = c(1e-6, 100), "nu.Unif" = c(0.5, 3.5)),
    starting = list("beta" = coef(lm(y ~ X)), "sigma.sq" = 6, "tau.sq" = 0.07,
        "phi" = 7.5, "nu" = 1.4),
    tuning = list("sigma.sq" = 0.1, "tau.sq" = 0.1, "phi" = 0.1, "nu" = 0.1))
ozone.mod = spLM(obs$ozone[ind] ~ 1+X, coords = coord, knots = sc.knots,
    n.samples = 10000,
    cov.model = "matern", modified.pp = TRUE,
    priors = list("beta.flat", "sigma.sq.ig" = c(1, 0.5), "tau.sq.ig" = c(1, 0.01),
        "phi.Unif" = c(1e-6, 100), "nu.Unif" = c(0.5, 3.5)),
    starting = list("beta" = coef(lm(y ~ X)), "sigma.sq" = 6, "tau.sq" = 0.07,
        "phi" = 7.5, "nu" = 1.4),
    tuning = list("sigma.sq" = 0.1, "tau.sq" = 0.1, "phi" = 0.1, "nu" = 0.1))
pred.X = sc.X
pred.X = cbind(1, pred.X[,2], pred.X[,2]^2, pred.X[,1]*pred.X[,2]*pred.X[,3])

ozone.pred.pp.samples = spPredict(ozone.pp, pred.coords = sc.X[,1:2], pred.covars = pred.X)
ozone.pred.mod.samples = spPredict(ozone.mod, pred.coords = sc.X[,1:2], pred.covars = pred.X)
pred.ozone.pp = apply(ozone.pred.pp.samples$p.y.predictive.samples, 1, mean)
pred.ozone.mod = apply(ozone.pred.mod.samples$p.y.predictive.samples, 1, mean)


### NO2 predictive process
# NO2: kmNorth, alt, kmEast* kmNorth * alt, kmEast^2
ind = which(!is.na(obs$no2))
y = obs$no2[ind]
X = sc.X[ind,]
X = cbind(X[,2], X[,3], X[,1]*X[,2]*X[,3], X[,1]^2)

#mod = spLM(obs$ozone[ind] ~ 1 + X, coords = as.matrix(obs[ind, c(4, 5)]),
coord = sc.X[ind,1:2]
no2.pp = spLM(y ~ 1+X, coords = coord, knots = sc.knots,
    n.samples = 5000, cov.model = "matern", modified.pp = FALSE,
    priors = list("beta.flat", "sigma.sq.ig" = c(1, 1), "tau.sq.ig" = c(1, 1),
        "phi.Unif" = c(1e-6, 100), "nu.Unif" = c(0.5, 7.5)),
    starting = list("beta" = coef(lm(y ~ X)), "sigma.sq" = 6, "tau.sq" = 0.07,
        "phi" = 7, "nu" = 1.4),
    tuning = list("sigma.sq" = 0.1, "tau.sq" = 0.5, "phi" = 0.1, "nu" = 0.1))
no2.mod = spLM(y ~ 1+X, coords = coord, knots = sc.knots,
    n.samples = 10000, cov.model = "matern", modified.pp = TRUE,
    priors = list("beta.flat", "sigma.sq.ig" = c(1, 1), "tau.sq.ig" = c(1, 1),
        "phi.Unif" = c(1e-6, 100), "nu.Unif" = c(0.5, 7.5)),
    starting = list("beta" = coef(lm(y ~ X)), "sigma.sq" = 6, "tau.sq" = 0.07,
        "phi" = 7, "nu" = 1.4),
    tuning = list("sigma.sq" = 0.1, "tau.sq" = 0.5, "phi" = 0.1, "nu" = 0.1))
pred.X = sc.X
pred.X = cbind(1, pred.X[,2], pred.X[,3], pred.X[,1]*pred.X[,2]*pred.X[,3], pred.X[,1]^2)

no2.pred.pp.samples = spPredict(no2.pp, pred.coords = sc.X[,1:2], pred.covars = pred.X)
no2.pred.mod.samples = spPredict(no2.mod, pred.coords = sc.X[,1:2], pred.covars = pred.X)
pred.no2.pp = apply(no2.pred.pp.samples$p.y.predictive.samples, 1, mean)
pred.no2.mod = apply(no2.pred.mod.samples$p.y.predictive.samples, 1, mean)

plot.nice(lonlat, pred.no2.pp, zscale = range(pred.no2.pp, pred.no2.mod),
    nlevels = 20)
plot.nice(lonlat, pred.no2.mod, zscale = range(pred.no2.pp, pred.no2.mod),
    nlevels = 20)

par(mfrow = c(2, 1))
plot.nice(lonlat, pred.no2.pp, zscale = range(pred.no2.pp, y),
    nlevels = 20)
plot.nice(lonlat[ind,], y, zscale = range(pred.no2.pp, y),
    nlevels = 20)
par(mfrow = c(1, 1))








mod = no2.mod
ss = seq(0, 1, length = nrow(mod$p.theta.samples))
pairs(as.matrix(mod$p.theta.samples), col = rgb(ss, 0, 0))
















plot(pred.ozone.pp, pred.ozone.mod)
abline(0, 1)

df.ozone.pp = data.frame("lon" = lonlat[,1], "lat" = lonlat[,2], "z" = pred.ozone.pp)
df.ozone.mod = data.frame("lon" = lonlat[,1], "lat" = lonlat[,2], "z" = pred.ozone.mod)
fld.ozone.pp = with(df.ozone.pp, interp(x = lon, y = lat, z = z))
fld.ozone.mod = with(df.ozone.mod, interp(x = lon, y = lat, z = z))



par(mfrow = c(2,1))







#pred.X = cbind(obs$kmEast, obs$kmNorth, obs$alt)
#pred.X = cbind(1, X[,2], X[,1]*X[,3], X[,2]*X[,3], X[,2]^2)
pred.X = sc.X
pred.X = cbind(1, pred.X[,2], pred.X[,2]^2, pred.X[,1]*pred.X[,2]*pred.X[,3])



lonlat = obs[, c(1, 2)]
col.obs = (y - min(y)) / diff(range(y))
col.pred = (pred.ozone - min(pred.ozone)) / diff(range(pred.ozone))

col.pred = rep(3, length(pred.ozone))
col.pred[-ind] = 2
pp = apply(pred.samples[[1]], 1, quantile, c(0.025, 0.975))

plot(y, pred.ozone[ind], ylim = range(pp[,ind]))
abline(0, 1)
segments(y, pp[1,ind], y, pp[2,ind])

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


### Form new observation matrix 
new.obs = obs
ind = which(is.na(new.obs$ozone))
new.obs$ozone[ind] = pred.ozone[ind]
