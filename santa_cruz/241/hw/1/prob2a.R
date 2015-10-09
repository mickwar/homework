library(MCMCpack)

cdf = function(x, t){
    out = length(t)
    for (i in 1:length(t))
        out[i] = mean(x <= t[i])
    return (out)
    }
stick_break = function(n, alpha, Gdraw){
    get_pi = function(n, alpha){
        V = c(0, rbeta(n, 1, alpha))
        return (V[-1] * cumprod(1-V[-(n+1)]))
        }
    pi_vec = get_pi(n, alpha)
    theta = Gdraw(n)
    # Getting more draws from G0 than necessary

    y = sample(theta, replace = TRUE, prob = pi_vec)
    return (y)
    }

tl = 500
alpha_vec = c(0.1, 1, 10)

pdf("./figs/prob2a.pdf", height = 18, width = 15)
par(mfrow = c(3,2), mar = c(4.1, 3.1, 3.1, 1.1), oma = c(0, 0, 3, 0))
for (alpha in alpha_vec){

    ### Ferguson
    set.seed(1)
    B = 50
    x = seq(-8, 8, length = tl)
    k = length(x)

    dG0 = pnorm(x, 0, 1) - pnorm(c(-Inf, x[-k]), 0, 1)
    dG = rdirichlet(B, alpha * dG0)
    G = apply(dG, 1, cumsum)
    Gmean = apply(G, 1, mean)
    Gvar = apply(G, 1, var)

    matplot(x, G, type = 'l', col = rgb(0.5, 0.5, 0.5, 0.4), lty = 1,
        main = bquote(alpha == .(alpha)), cex.main = 3.0, xlab = "")
    curve(pnorm(x, 0, 1), col = rgb(1, 0, 0, 0.7), add = TRUE, lwd = 3)
    lines(x, Gmean, col = rgb(0, 0, 1, 0.7), lwd = 3)
    lines(x, Gvar, col = rgb(0, 1, 0, 0.7), lwd = 3)
    #matplot(x, G[,1:5], type = 'l', col = rgb(0.5, 0.0, 0.5), lty = 1, add = TRUE)
    legend(min(x), 1, legend = c("Prior draws", "Estimated prior mean", "Estimated prior variance",
        "CDF of N(0,1)"), lwd = c(1, 3, 3, 3), lty = c(1, 1, 1, 1),
        col = c("gray50", "blue", "green", "red"))

    ### Constructive
    n = 850
    set.seed(1)
    x = matrix(0, B, n)
    for (i in 1:B)
        x[i,] = stick_break(n, alpha, function(n) rnorm(n, 0, 1))

    tvec = seq(-8, 8, length = tl)
    q = apply(x, 1, function(x) cdf(x, tvec))
    qlines = apply(q, 1, quantile, c(0.025, 0.5, 0.975))
    qmean = apply(q, 1, mean)
    qvar = apply(q, 1, var)

    matplot(tvec, q, type = 'l', lty = 1, col = rgb(0.5, 0.5, 0.5, 0.4),
        main = bquote(alpha == .(alpha)), cex.main = 3.0, xlab = "")
    curve(pnorm(x, 0, 1), col = rgb(1, 0, 0, 0.7), add = TRUE, lwd = 3)
    lines(tvec, qmean, col = rgb(0, 0, 1, 0.7), lwd = 3)
    lines(tvec, qvar, col = rgb(0, 1, 0, 0.7), lwd = 3)
    legend(min(tvec), 1, legend = c("Prior draws", "Estimated prior mean",
        "Estimated prior variance", "CDF of N(0,1)"), lwd = c(1, 3, 3, 3), lty = c(1, 1, 1, 1),
        col = c("gray50", "blue", "green", "red"))


    }
title(main = "Ferguson's definition                               Stick breaking    ", outer = TRUE, cex.main = 4)
dev.off()

