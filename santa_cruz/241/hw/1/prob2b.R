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
#alpha

a_vec = c(1, 50, 1/5)
b_vec = c(10, 1, 1/100)

 pdf("./figs/prob2b.pdf", height = 18, width = 12)
par(mfrow = c(3,1), mar = c(4.1, 3.1, 3.1, 1.1), oma = c(0, 0, 3, 0))
for (j in 1:3){

    ### Constructive
    set.seed(1)
    B = 400
    n = 850
    x = matrix(0, B, n)
    alpha = rgamma(B, a_vec[j], b_vec[j])
    for (i in 1:B)
        x[i,] = stick_break(n, alpha[i], function(n) rnorm(n, 0, 1))

    tvec = seq(-8, 8, length = tl)
    q = apply(x, 1, function(x) cdf(x, tvec))
    qlines = apply(q, 1, quantile, c(0.025, 0.5, 0.975))
    qmean = apply(q, 1, mean)
    qvar = apply(q, 1, var)

    opts = par(no.readonly = TRUE)

    matplot(tvec, q, type = 'l', lty = 1, col = rgb(0.5, 0.5, 0.5, 0.4), cex.main = 3.0, xlab = "",
        main = bquote(alpha ~ "~" ~ Gamma(.(a_vec[j]), .(b_vec[j]))))
    curve(pnorm(x, 0, 1), col = rgb(1, 0, 0, 0.7), add = TRUE, lwd = 3)
    lines(tvec, qmean, col = rgb(0, 0, 1, 0.7), lwd = 3)
    lines(tvec, qvar, col = rgb(0, 1, 0, 0.7), lwd = 3)
    legend(min(tvec), 1, legend = c("Prior draws", "Estimated prior mean",
        "Estimated prior variance", "CDF of N(0,1)"), lwd = c(1, 3, 3, 3), lty = c(1, 1, 1, 1),
        col = c("gray50", "blue", "green", "red"))


    dens = density(alpha)
    rng.x = range(tvec)
    x.diff = diff(rng.x)
  
    left = rng.x[1] + x.diff*2/3
    right = rng.x[2]
    par(fig = c(grconvertX(c(left, right), from="user", to="ndc"),
        grconvertY(c(0.25, 0.65), from="user", to="ndc")),
        mar = c(0.1, 0.1, 1.0, 0.1), new = TRUE)
    #plot(density(x),col="blue",cex.main=.5,lwd=3)
    plot(dens, type="l", col="gray20", cex.main=1.5, axes=FALSE,
        main = expression(paste("Prior for ", alpha)))
    axis(1, cex.axis=1.2)
    par(opts)
    par("mfg" = c(min(j+1, 3), 1, 3, 1))

    }
 dev.off()

