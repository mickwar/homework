library(geoR)

set.seed(9)
type = c("spherical", "powered.exponential", "cauchy", "wave", "matern", "matern")
kappa = c(0.5, 2, 1, 0.5, 0.5, 1.5)
ranges = double(length(type))
mains = type
mains[2] = paste(mains[2], "nu=2")
mains[3] = "rational.quadratic"
mains[5] = paste(mains[5], "nu=1/2")
mains[6] = paste(mains[6], "nu=3/2")

# Solve for the range yielding 0.05 correlation at unit distance
for (i in 1:length(type)){
    ranges[i] = uniroot(f = function(x){
        cov.spatial(1, cov.model = type[i], cov.pars = c(1, x), kappa = kappa[i]) - 0.05},
        c(0.001, 10))$root
    }

# Covariance functions
pdf("./figs/covariance.pdf", width = 9, height = 12)
par(mfrow = c(3,2))
for (i in 1:length(type)){
    h = seq(0, 5, length = 1000)
    y = cov.spatial(h, cov.model = type[i], cov.pars = c(1, ranges[i]), kappa = kappa[i])
    plot(h, y, type = 'l', lwd = 3, bty = 'n', main = mains[i])
    legend("topright", legend = paste0("phi=",round(ranges[i], 3)), bty = 'n')
    }
par(mfrow = c(1,1))
dev.off()


# Semi-variograms
pdf("./figs/semi-variogram.pdf", width = 9, height = 12)
par(mfrow = c(3,2))
for (i in 1:length(type)){
    h = seq(0, 5, length = 1000)
    y = cov.spatial(h, cov.model = type[i], cov.pars = c(1, ranges[i]), kappa = kappa[i])
    plot(h, 1*(1-y), type = 'l', lwd = 3, bty = 'n', main = mains[i])
    legend("topleft", legend = paste0("phi=",round(ranges[i], 3)), bty = 'n')
    }
par(mfrow = c(1,1))
dev.off()

# Simulations
pdf("./figs/simulations.pdf", width = 9, height = 12)
par(mfrow = c(3,2))
for (i in 1:length(type)){
    xx = seq(-5, 5, length = 100)
    y = grf(grid = cbind(xx, 0), cov.model = type[i], cov.pars = c(1, ranges[i]), kappa = kappa[i])
    plot(xx, y$data, type='l', bty = 'n', main = mains[i])
    legend("topleft", legend = paste0("phi=",round(ranges[i], 3)), bty = 'n')
    }
dev.off()

