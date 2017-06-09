library(mwBASE)
library(MASS)

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
bezier = function(s, phi, nu)
    ifelse(abs(s) < sqrt(phi), (1 - (s^2)/phi)^nu, 0)

obs = read.table("./obs.txt", header = TRUE)
knots = read.table("./knots.txt", header = TRUE)

obs = obs[!is.na(obs$ozone),-c(7,8)]

# Remove outliers and leverage points (make them NA)
obs = obs[-c(31, 72, 77, 103),]

sc.X = cbind(scale(obs$kmEast), scale(obs$kmNorth),
    scale(obs$alt))
sc.knots = cbind(scale(knots$kmEast), scale(knots$kmNorth))



m = nrow(obs)
p = nrow(knots)
# dd = t(tail(as.matrix(dist(rbind(sc.X[,1:2], sc.knots))), p)[,1:m])
dd = t(tail(as.matrix(dist(rbind(obs[,4:5], knots[,4:5]))), p)[,1:m])
colnames(dd) = NULL
rownames(dd) = NULL

sum(bezier(dd, 57, 1) > 0 )
bezier(dd[1,], 1)

# Ozone: kmNorth, kmNorth^2, (kmEast * kmNorth * alt)
mu.X = cbind(1, sc.X[,2], sc.X[,2]^2, sc.X[,1] * sc.X[,2] * sc.X[,3])

dat = list("y" = obs$ozone, "X" = mu.X, "knots" = knots)

calc.post = function(dat, param){
    y = dat$y
    X = dat$X
    n = length(y)
    p = nrow(dat$knots)
    k = ncol(X)
    beta = param[1:k]
    w = param[(k+1):(k+p)]
    phi = param[k+p+1]
    nu = param[k+p+2]
    sig2 = param[k+p+3]
    tau2 = param[k+p+4]

    if (nu <= 0 || sig2 <= 0 || tau2 <= 0 || phi <= 0)
        return (-Inf)

    K.psi = bezier(dd, phi, nu)

    if (sum(K.psi > 0) == 0)
        return (-Inf)

    # likelihood
    out = -n/2 * log(tau2) -sum((y - X %*% beta - K.psi %*% w)^2) / (2 * tau2)

    # priors
    out = out + sum(dnorm(beta, 0, sqrt(10), log = TRUE))
    out = out + sum(dnorm(w, 0, sqrt(sig2), log = TRUE))
#   out = out + dgamma(phi, 1, 0.1, log = TRUE)
    out = out - log(phi)
    out = out + dgamma(nu, 1, 1, log = TRUE)
    out = out + dgamma(sig2, 1, 0.1, log = TRUE)
    out = out + dgamma(tau2, 1, 0.1, log = TRUE)

    return (out)
    }

p = nrow(knots)
k = ncol(mu.X)

mcmc = mcmc_sampler(dat, calc.post, nparam = k + p + 4, nburn = 10000, nmcmc = 10000,
    chain_init = c(runif(k+p), 300, runif(3)))

mean(mcmc$accept)
plot(mcmc$params[,5], type = 'l')

ss = seq(0, 1, length = nrow(mcmc$params))
pairs(mcmc$params[,1:k], col = rgb(ss, 0, 0))

pairs(mcmc$params[,k+sample(p, 4)], col = rgb(ss, 0, 0))

pairs(mcmc$params[,(k+p+1):(k+p+4)], col = rgb(ss, 0, 0))

plot_hpd(mcmc$params[,k+p+1], col1 = 'dodgerblue')

plot_hpd(mcmc$params[,k+p+2], col1 = 'dodgerblue')

plot_hpd(mcmc$params[,k+p+3], col1 = 'dodgerblue')

plot_hpd(mcmc$params[,k+p+4], col1 = 'dodgerblue')


xx = sort(unique(knots[,1]))
yy = sort(unique(knots[,2]))

W = 0*diag(p)
for (i in 1:(p-1)){
    for (j in (i+1):p){
        tmpix = which(knots[i,1] == xx)
        tmpiy = which(knots[i,2] == yy)
        tmpjx = which(knots[j,1] == xx)
        tmpjy = which(knots[j,2] == yy)
        dx = tmpix - tmpjx
        dy = tmpiy - tmpjy
        plot(knots[,1:2])
#       points(xx[tmpix], yy[tmpiy], col = 'red', pch = 16)
#       points(xx[tmpjx], yy[tmpjy], col = 'blue', pch = 16)
        if (abs(dx) + abs(dy) <= 1){
            W[i,j] = -1
            W[j,i] = -1
            cat("Neighbors:", i, j, "\n")
            }
#       readline()
        }
    }

diag(W) = apply(W, 1, function(x) sum(x == -1))

W[1:6,1:6]


plot(knots[,1:2])
points(obs[,1:2], col = 'red')

plot(sc.knots)
points(sc.X[,1:2], col = 'red')
