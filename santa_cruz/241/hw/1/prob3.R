library(MCMCpack)

cdf = function(x, t){
    out = length(t)
    for (i in 1:length(t))
        out[i] = mean(x <= t[i])
    return (out)
    }
gen.y1 = function(n)
    rnorm(n)
gen.y2 = function(n){
    out = double(n)
    mix = sample(1:3, n, replace = TRUE, prob = c(0.5, 0.3, 0.2))
    out = ifelse(mix == 1, rnorm(sum(mix == 1), -2.5, 0.5), out)
    out = ifelse(mix == 2, rnorm(sum(mix == 2), 0.5, 0.7), out)
    out = ifelse(mix == 3, rnorm(sum(mix == 3), 1.5, 2.0), out)
    }

B = 1000

#alpha = 100; m = -4; s = 2
#alpha = 100; m = -1; s = 0.1
 alpha = 1; m = -4; s = 2
G0.dist = function(t)
    pnorm(t, m, s)

G0.dist.post = function(t)
    alpha/(alpha+n)*G0.dist(t) + 1/(alpha+n)*sapply(t, function(t) sum(y <= t))

#tvec = seq(m-3*s, m+s*3, length = 1000)
tvec = seq(-9, 5, length = 1000)
n = c(20, 200, 2000)


pdf("./figs/prob3c.pdf", height = 12, width = 10)
par(mfrow = c(3,2), mar = c(4.1, 3.1, 3.1, 1.1), oma = c(0, 0, 3, 0))
for (n in c(20, 200, 2000)){
    ### Ferguson

    # posterior
    set.seed(1)
    y = gen.y1(n)

    dG0 = G0.dist.post(c(tvec[-1], Inf)) - G0.dist.post(tvec)
    dG.pre = rdirichlet(B, (alpha + n)*dG0)

    G.draws = apply(dG.pre, 1, cumsum)
    Gmean = apply(G.draws, 1, mean)
    Glines = apply(G.draws, 1, quantile, c(0.025, 0.975))

    plot(0, xlim = range(tvec), ylim = c(0, 1), cex.main = 3.0, xlab = "", type='n',
        main = bquote(alpha == .(alpha) ~ "," ~ n == .(n)))
    polygon(c(tvec, rev(tvec)), c(Glines[1,], rev(Glines[2,])), col = 'steelblue', border = NA)
    curve(G0.dist(x), col = rgb(1, 0, 0, 0.7), add = TRUE, lwd = 3)
    lines(tvec, Gmean, col = rgb(0, 0, 1, 0.7), lwd = 3)
    legend(min(tvec), 1, legend = c("eCDF of data", "Estimated posterior mean", 
        "95% p.w. credible intervals", paste0("CDF of baseline: N(", m, ",", s, ")")),
        lwd = c(1, 3, 3, 3), lty = c(1, 1, 1, 1), col = c(rgb(0,0,0), rgb(0,0,1,0.7),
        "steelblue", rgb(1,0,0,0.7)), pch = c(20, NA, NA, NA))
    lines(ecdf(y), col.01line = NA, verticals = TRUE)

    # posterior
    set.seed(1)
    y = gen.y2(n)

    dG0 = G0.dist.post(c(tvec[-1], Inf)) - G0.dist.post(tvec)
    dG.pre = rdirichlet(B, (alpha + n)* dG0)

    G.draws = apply(dG.pre, 1, cumsum)
    Gmean = apply(G.draws, 1, mean)
    Glines = apply(G.draws, 1, quantile, c(0.025, 0.975))

    plot(0, xlim = range(tvec), ylim = c(0, 1), cex.main = 3.0, xlab = "", type='n',
        main = bquote(alpha == .(alpha) ~ "," ~ n == .(n)))
    polygon(c(tvec, rev(tvec)), c(Glines[1,], rev(Glines[2,])), col = 'steelblue', border = NA)
    curve(G0.dist(x), col = rgb(1, 0, 0, 0.7), add = TRUE, lwd = 3)
    lines(tvec, Gmean, col = rgb(0, 0, 1, 0.7), lwd = 3)
    legend(min(tvec), 1, legend = c("eCDF of data", "Estimated posterior mean", 
        "95% p.w. credible intervals", paste0("CDF of baseline: N(", m, ",", s, ")")),
        lwd = c(1, 3, 3, 3), lty = c(1, 1, 1, 1), col = c(rgb(0,0,0), rgb(0,0,1,0.7),
        "steelblue", rgb(1,0,0,0.7)), pch = c(20, NA, NA, NA))
    lines(ecdf(y), col.01line = NA, verticals = TRUE)

    }
title("               Data: N(0,1)                     Data: Mixture of normals", outer = TRUE, cex.main = 3)
dev.off()

