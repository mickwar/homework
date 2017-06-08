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
obs$no2[23] = NA
obs$no2[49] = NA
obs$no2[104] = NA

# Remove rows that now have all NA in ozone, no2, pm25
rm = apply(obs[,6:8], 1, function(x) all(is.na(x)))
obs = obs[!rm,]

# Scaled
sc.X = cbind(scale(obs$kmEast), scale(obs$kmNorth),
    scale(obs$alt))
sc.knots = cbind(scale(knots$kmEast), scale(knots$kmNorth))

#
lonlat = as.matrix(obs[,1:2])






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
# NO2: kmNorth, alt, kmEast* kmNorth

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
    n.samples = 10000, cov.model = "matern", modified.pp = FALSE,
    priors = list("beta.flat", "sigma.sq.ig" = c(1, 0.5), "tau.sq.ig" = c(1, 0.01),
        "phi.Unif" = c(1e-6, 100), "nu.Unif" = c(0.5, 3.5)),
    starting = list("beta" = coef(lm(y ~ X)), "sigma.sq" = 6, "tau.sq" = 0.07,
        "phi" = 7.5, "nu" = 1.4),
    tuning = list("sigma.sq" = 0.1, "tau.sq" = 0.1, "phi" = 0.1, "nu" = 0.1))
ozone.mod = spLM(obs$ozone[ind] ~ 1+X, coords = coord, knots = sc.knots,
    n.samples = 10000, cov.model = "matern", modified.pp = TRUE,
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
# NO2: kmNorth, alt, kmEast* kmNorth
ind = which(!is.na(obs$no2))
y = obs$no2[ind]
X = sc.X[ind,]
# X = cbind(X[,2], X[,3], X[,1]*X[,2]*X[,3], X[,1]^2)
X = cbind(X[,2], X[,3], X[,1]*X[,2])

#mod = spLM(obs$ozone[ind] ~ 1 + X, coords = as.matrix(obs[ind, c(4, 5)]),
coord = sc.X[ind,1:2]
no2.pp = spLM(y ~ 1+X, coords = coord, knots = sc.knots,
    n.samples = 10000, cov.model = "matern", modified.pp = FALSE,
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
# pred.X = cbind(1, pred.X[,2], pred.X[,3], pred.X[,1]*pred.X[,2]*pred.X[,3], pred.X[,1]^2)
pred.X = cbind(1, pred.X[,2], pred.X[,3], pred.X[,1]*pred.X[,2])

no2.pred.pp.samples = spPredict(no2.pp, pred.coords = sc.X[,1:2], pred.covars = pred.X)
no2.pred.mod.samples = spPredict(no2.mod, pred.coords = sc.X[,1:2], pred.covars = pred.X)
pred.no2.pp = apply(no2.pred.pp.samples$p.y.predictive.samples, 1, mean)
pred.no2.mod = apply(no2.pred.mod.samples$p.y.predictive.samples, 1, mean)

### PM25 predictive process
# PM25: kmEast, kmNorth, alt, (kmNorth^2)
ind = which(!is.na(obs$pm25))
y = obs$pm25[ind]
X = sc.X[ind,]
X = cbind(X[,1], X[,2], X[,3], X[,2]^2)

#mod = spLM(obs$ozone[ind] ~ 1 + X, coords = as.matrix(obs[ind, c(4, 5)]),
coord = sc.X[ind,1:2]
pm25.pp = spLM(y ~ 1+X, coords = coord, knots = sc.knots,
    n.samples = 10000, cov.model = "matern", modified.pp = FALSE,
    priors = list("beta.flat", "sigma.sq.ig" = c(1, 1), "tau.sq.ig" = c(1, 1),
        "phi.Unif" = c(1e-6, 100), "nu.Unif" = c(0.5, 7.5)),
    starting = list("beta" = coef(lm(y ~ X)), "sigma.sq" = 6, "tau.sq" = 0.07,
        "phi" = 7, "nu" = 1.4),
    tuning = list("sigma.sq" = 0.1, "tau.sq" = 0.5, "phi" = 0.1, "nu" = 0.1))
pm25.mod = spLM(y ~ 1+X, coords = coord, knots = sc.knots,
    n.samples = 10000, cov.model = "matern", modified.pp = TRUE,
    priors = list("beta.flat", "sigma.sq.ig" = c(1, 1), "tau.sq.ig" = c(1, 1),
        "phi.Unif" = c(1e-6, 100), "nu.Unif" = c(0.5, 7.5)),
    starting = list("beta" = coef(lm(y ~ X)), "sigma.sq" = 6, "tau.sq" = 0.07,
        "phi" = 7, "nu" = 1.4),
    tuning = list("sigma.sq" = 0.1, "tau.sq" = 0.5, "phi" = 0.1, "nu" = 0.1))
pred.X = sc.X
pred.X = cbind(1, pred.X[,1], pred.X[,2], pred.X[,3], pred.X[,2]^2)

pm25.pred.pp.samples = spPredict(pm25.pp, pred.coords = sc.X[,1:2], pred.covars = pred.X)
pm25.pred.mod.samples = spPredict(pm25.mod, pred.coords = sc.X[,1:2], pred.covars = pred.X)
pred.pm25.pp = apply(pm25.pred.pp.samples$p.y.predictive.samples, 1, mean)
pred.pm25.mod = apply(pm25.pred.mod.samples$p.y.predictive.samples, 1, mean)





pdf("./figs/ozone_pp.pdf", width = 8, height = 12)
par(mfrow = c(2, 1))
plot.nice(lonlat, pred.ozone.pp, zscale = range(pred.ozone.pp, pred.ozone.mod),
    nlevels = 20, main = "Ozone - Predictive Process")
plot.nice(lonlat, pred.ozone.mod, zscale = range(pred.ozone.pp, pred.ozone.mod),
    nlevels = 20, main = "Ozone - Modified Predictive Process")
par(mfrow = c(1, 1))
dev.off()

pdf("./figs/no2_pp.pdf", width = 8, height = 12)
par(mfrow = c(2, 1))
plot.nice(lonlat, pred.no2.pp, zscale = range(pred.no2.pp, pred.no2.mod,y ),
    nlevels = 20, main = "NO2 - Predictive Process")
plot.nice(lonlat, pred.no2.mod, zscale = range(pred.no2.pp, pred.no2.mod,y ),
    nlevels = 20, main = "NO2 - Modified Predictive Process")
par(mfrow = c(1, 1))
dev.off()

pdf("./figs/pm25_pp.pdf", width = 8, height = 12)
par(mfrow = c(2, 1))
plot.nice(lonlat, pred.pm25.pp, zscale = range(pred.pm25.pp, pred.pm25.mod),
    nlevels = 20, main = "PM2.5 - Predictive Process")
plot.nice(lonlat, pred.pm25.mod, zscale = range(pred.pm25.pp, pred.pm25.mod),
    nlevels = 20, main = "PM2.5 - Modified Predictive Process")
par(mfrow = c(1, 1))
dev.off()








### Form new observation matrix 
new.obs = obs
ind = which(is.na(new.obs$ozone))
new.obs$ozone[ind] = pred.ozone.mod[ind]
ind = which(is.na(new.obs$no2))
new.obs$no2[ind] = pred.no2.mod[ind]
ind = which(is.na(new.obs$pm25))
new.obs$pm25[ind] = pred.pm25.mod[ind]


### Multivariate fit
# Ozone: kmNorth, kmNorth^2, (kmEast * kmNorth * alt)
# NO2: kmNorth, alt, kmEast* kmNorth
# PM25: kmEast, kmNorth, alt, (kmNorth^2)
y.o = new.obs$ozone
y.n = new.obs$no2
y.p = new.obs$pm25
X = sc.X
X.o = cbind(X[,2], X[,2]^2, X[,1]*X[,2]*X[,3])
X.n = cbind(X[,2], X[,3], X[,1]*X[,2])
X.p = cbind(X[,1], X[,2], X[,3], X[,2]^2)

coord = sc.X[,1:2]
mvmod = spMvLM(list(y.o ~ 1 + X.o, y.n ~ 1 + X.n, y.p ~ 1 + X.p), n.samples = 5000,
    coords = coord, knots = sc.knots, modified.pp = TRUE, cov.model = "matern",
    priors = list("beta.flat", "K.iw" = list(5, 5*diag(3)),
        "phi.Unif" = list(rep(1e-6, 3), rep(10, 3)),
        "nu.Unif" = list(rep(0.5, 3), rep(7.5, 3)),
        "Psi.iw" = list(10, diag(3))),
    starting = list("beta" = c(coef(lm(y.o ~ 1 + X.o)), coef(lm(y.n ~ 1 + X.n)),
        coef(lm(y.p ~ 1 + X.p))),
        "A" = c(2.5, 0, 0, 2.4, 0, 0.7),
        "L" = c(1, -0.3, -0.1, 2.5, 0.3, 0.5),
        "phi" = c(6, 6, 6),
        "nu" = c(2.5, 2.0, 0.55)),
    tuning = list("A" = c(0.005, 0.01, 0.005, 0.01, 0.005, 0.005), "phi" = c(0.01, 0.01, 0.01),
        "nu" = c(0.01, 0.01, 0.01), "L" = c(0.01, 0.01, 0.005, 0.005, 0.005, 0.005)))

length(unique(mvmod$p.theta.samples[,1])) / nrow(mvmod$p.theta.samples)
dim(mvmod$p.theta.samples)
colMeans(mvmod$p.theta.samples)
tail(mvmod$p.theta.samples, 0)
plot(mvmod$p.theta.samples[,7])


mvpred.samples = spPredict(mvmod, pred.coord = coord, pred.covars = cbind(1, X.o, 1, X.n, 1, X.p))

mvpred.out = apply(mvpred.samples[[1]], 1, mean)

pred.ozone.mv = mvpred.out[1:158]
pred.no2.mv = mvpred.out[159:316]
pred.pm25.mv = mvpred.out[317:474]

ind.o = is.na(pred.ozone.mv) | (pred.ozone.mv < 0) | (pred.ozone.mv > 1000)
pred.ozone.mv[!ind]

ind.n = is.na(pred.no2.mv) 
# | (pred.no2.mv < 0) | (pred.no2.mv > 1000)
pred.no2.mv[!ind.n]

ind.p = is.na(pred.pm25.mv) | (pred.pm25.mv < -1000)
pred.pm25.mv[!ind.p]




pdf("./figs/knots.pdf", width = 9, height = 9)
plot(lonlat, col = 'blue', pch = 16, bty = 'n',
    main = "Locations and Knots", xlab = "Longitude", ylab = "Latitude")
points(knots[,1:2], col = 'orange', pch = 15)
map("state", add = TRUE)
dev.off()




pdf("./figs/ozone_pp.pdf", width = 9, height = 15)
par(mfrow = c(3, 1))
plot.nice(lonlat, pred.ozone.pp, zscale = range(pred.ozone.pp, pred.ozone.mod),
    nlevels = 20, main = "Ozone - PP")
plot.nice(lonlat, pred.ozone.mod, zscale = range(pred.ozone.pp, pred.ozone.mod),
    nlevels = 20, main = "Ozone - Modified PP")
plot.nice(lonlat[!ind.o,], pred.ozone.mv[!ind.o], 
    nlevels = 20, main = "Ozone - Multivariate Modified PP")
par(mfrow = c(1, 1))
dev.off()

pdf("./figs/no2_pp.pdf", width = 9, height = 15)
par(mfrow = c(3, 1))
plot.nice(lonlat, pred.no2.pp, zscale = range(pred.no2.pp, pred.no2.mod),
    nlevels = 20, main = "NO2 - PP")
plot.nice(lonlat, pred.no2.mod, zscale = range(pred.no2.pp, pred.no2.mod),
    nlevels = 20, main = "NO2 - Modified PP")
plot.nice(lonlat[!ind.n,], pred.no2.mv[!ind.n],
    nlevels = 20, main = "NO2 - Multivariate Modified PP")
par(mfrow = c(1, 1))
dev.off()

pdf("./figs/pm25_pp.pdf", width = 9, height = 15)
par(mfrow = c(3, 1))
plot.nice(lonlat, pred.pm25.pp, zscale = range(pred.pm25.pp, pred.pm25.mod),
    nlevels = 20, main = "PM2.5 - PP")
plot.nice(lonlat, pred.pm25.mod, zscale = range(pred.pm25.pp, pred.pm25.mod),
    nlevels = 20, main = "PM2.5 - Modified PP")
plot.nice(lonlat[!ind.p,], pred.pm25.mv[!ind.p],
    nlevels = 20, main = "PM2.5 - Multivariate Modified PP")
par(mfrow = c(1, 1))
dev.off()


pdf("./figs/fit.pdf", width = 9, height = 15)
par(mfrow = c(3, 1))
ind = which(!is.na(obs$ozone))
y = obs$ozone[ind]
pred = pred.ozone.mod[ind]
samples = ozone.pred.mod.samples[[1]][ind,]
qq = apply(samples, 1, quantile, c(0.025, 0.975))
plot(y, pred, ylim = range(qq), pch = 16, bty = 'n',
    xlab = "Data", ylab = "Prediction", main = "Ozone - Modified PP")
segments(y, qq[1,], y, qq[2,])
abline(0, 1)
 
ind = which(!is.na(obs$no2))
y = obs$no2[ind]
pred = pred.no2.mod[ind]
samples = no2.pred.mod.samples[[1]][ind,]
qq = apply(samples, 1, quantile, c(0.025, 0.975))
plot(y, pred, ylim = range(qq), pch = 16, bty = 'n',
    xlab = "Data", ylab = "Prediction", main = "NO2 - Modified PP")
segments(y, qq[1,], y, qq[2,])
abline(0, 1)

ind = which(!is.na(obs$pm25))
y = obs$pm25[ind]
pred = pred.pm25.mod[ind]
samples = pm25.pred.mod.samples[[1]][ind,]
qq = apply(samples, 1, quantile, c(0.025, 0.975))
plot(y, pred, ylim = range(qq), pch = 16, bty = 'n',
    xlab = "Data", ylab = "Prediction", main = "PM2.5 - Modified PP")
segments(y, qq[1,], y, qq[2,])
abline(0, 1)
par(mfrow = c(1, 1))
dev.off()


sum(!is.na(obs$ozone))
sum(!is.na(obs$no2))
sum(!(is.na(obs$ozone) | is.na(obs$no2)))
m(!is.na(obs$no2))
