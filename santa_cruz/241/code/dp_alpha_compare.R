library(MCMCpack)
cdf = function(x, t){
    out = length(t)
    for (i in 1:length(t))
        out[i] = mean(x <= t[i])
    return (out)
    }


alpha_vec = c(1, 10, 100, 100)
n_vec = c(15, 15, 15, 1000)

#pdf("./dp_comp_all.pdf", height = 24, width = 12)
#par(mfrow = c(4,2), mar = c(4.1, 3.1, 3.1, 1.1))
for (i in 1:length(alpha_vec)){
    alpha = alpha_vec[i]
    n = n_vec[i]
#   pdf(paste0("./dp_comp_", alpha, "_", n, ".pdf"), height = 6, width = 12)
#   par(mfrow = c(1,2), mar = c(4.1, 3.1, 3.1, 1.1))

    set.seed(1)
    y = rnorm(n, 4, 2)

    B = 1000
    x = seq(-8, 8, length = 200)
    k = length(x)

    ### Prior
    dG0 = pnorm(x, 0, 1) - pnorm(c(-Inf, x[-k]), 0, 1)
    dG = rdirichlet(B, alpha * dG0)
    G = apply(dG, 1, cumsum)
    Gmean = apply(G, 1, mean)

    matplot(x, G, type = 'l', col = rgb(0.5, 0.5, 0.5, 0.3), lty = 1, main = "Prior")
    curve(pnorm(x, 0, 1), add = TRUE, lwd = 3)
    lines(x, Gmean, col = 'green', lty = 2, lwd = 3)
    #matplot(x, G[,1:5], type = 'l', col = rgb(0.5, 0.0, 0.5), lty = 1, add = TRUE)
    legend(min(x), 1, legend = c("Prior draws", "Prior mean", "True mean"),
        lwd = c(1, 3, 3), lty = c(1, 2, 1), col = c("gray50", "green", "black"))


    ### Posterior
    dG0 = pnorm(x, 0, 1) - pnorm(c(-Inf, x[-k]), 0, 1)
    dFn = cdf(y, x) - cdf(y, c(-Inf, x[-k]))
    dG0 = (alpha*dG0 + n*dFn) / (alpha + n)
    dG = rdirichlet(B, (alpha+n) * dG0)
    G = apply(dG, 1, cumsum)
    Gmean = apply(G, 1, mean)

    matplot(x, G, type = 'l', col = rgb(0.5, 0.5, 0.5, 0.3), lty = 1, main = "Posterior")
    curve(pnorm(x, 0, 1), add = TRUE, lwd = 3)
    lines(x, Gmean, col = 'green', lty = 2, lwd = 3)
    #matplot(x, G[,1:5], type = 'l', col = rgb(0.5, 0.0, 0.5), lty = 1, add = TRUE)
    points(x, cdf(y, x), col = 'red', pch = 20, cex = 0.5)
    legend(min(x), 1, legend = c("ecdf of data", "Posterior draws", "Posterior mean", "True mean"),
        lwd = c(NA, 1, 3, 3), lty = c(NA, 1, 2, 1), pch = c(20, NA, NA, NA), 
        col = c("red", "gray50", "green", "black"))
#   title(paste0("Alpha = ", alpha, " -- ", "nobs = ", n), outer = TRUE, line = -1.5,
#       cex.main = 2)
    title(paste0("Alpha = ", alpha, " -- ", "nobs = ", n), outer = TRUE, line = -1.5-45.2*(i-1),
        cex.main = 2)

#   dev.off()
    }
#dev.off()

fq = apply(G, 1, quantile, c(0.025, 0.5, 0.975))
par(mfrow = c(1,1), mar = c(5.1, 4.1, 4.1, 2.1))
plot(x, fq[2,], lwd = 3, ylim = range(fq), type = 'l')
lines(x, fq[1,], lwd = 1)
lines(x, fq[3,], lwd = 1)
plot(ecdf(y), col = 'red', add = TRUE)


fq = apply(dG, 2, quantile, c(0.025, 0.5, 0.975))
par(mfrow = c(1,1), mar = c(5.1, 4.1, 4.1, 2.1))
plot(x, fq[2,], lwd = 3, ylim = range(fq), type = 'l')
lines(x, fq[1,], lwd = 1)
lines(x, fq[3,], lwd = 1)
lines(density(y), col = 'red')
points(y, rep(0.1, n), col = 'blue', pch = 20)


