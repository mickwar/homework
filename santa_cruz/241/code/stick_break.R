use = function(name, a, b){
    Gdraw <<- function(n)
        get(paste0("r", name))(n, a, b)
    Gp <<- function(n)
        get(paste0("p", name))(n, a, b)
    Gf <<- function(x)
        get(paste0("d", name))(x, a, b)
    }

#use("beta", 4, 1)
#use("gamma", 5, 1/30)
#use("t", 5, 0)
 use("norm", 0, 1)
#use("unif", 0, 1)


m = 1000
alpha = 100
n = 1000
tl = 500


cdf = function(x, t){
    out = length(t)
    for (i in 1:length(t))
        out[i] = mean(x <= t[i])
    return (out)
    }

dp_draw = function(n, alpha, Gdraw){
    get_pi = function(n, alpha){
        V = c(0, rbeta(n, 1, alpha))
        return (V[-1] * cumprod(1-V[-(n+1)]))
        }
    pi_vec = get_pi(n, alpha)
    theta = Gdraw(n)
    # Getting more draws from G0 than necessary

#   y = sort(sample(theta, replace = TRUE, prob = pi_vec))
    y = sample(theta, replace = TRUE, prob = pi_vec)
    return (y)
    }

x = matrix(0, m, n)


for (i in 1:m)
    x[i,] = dp_draw(n, alpha, Gdraw)
#x = t(replicate(m, dp_draw(n, alpha, Gdraw)))

tvec = seq(min(x), max(x), length = tl)
q = apply(x, 1, function(x) cdf(x, tvec))
qlines = apply(q, 1, quantile, c(0.025, 0.5, 0.975))
qmean = apply(q, 1, mean)

bw = 1.06*sd(x)/(n^(1/5))
f = apply(x, 1, function(x) density(x, bw = bw, from = min(tvec), to = max(tvec), n = length(tvec))$y)
flines = apply(f, 1, quantile, c(0.025, 0.5, 0.975))
fmean = apply(f, 1, mean)

par(mfrow = c(2, 1), mar = c(4.1, 3.1, 2.1, 1.1))
matplot(tvec, q[,1:100], type = 'l', lty = 1, col = rgb(0.5, 0.5, 0.5, 0.3))
curve(Gp(x), col = rgb(1, 0, 0, 0.7), add = TRUE, lwd = 3)
lines(tvec, qmean, col = rgb(0, 0, 1, 0.7), lwd = 3)
lines(tvec, qlines[2,], col = rgb(0, 0.4, 0, 0.7), lwd = 3)
lines(tvec, qlines[1,], col = rgb(0, 0.8, 0, 0.7), lwd = 2, lty = 2)
lines(tvec, qlines[3,], col = rgb(0, 0.8, 0, 0.7), lwd = 2, lty = 2)

matplot(tvec, f[,1:100], type = 'l', col = rgb(0.5, 0.5, 0.5, 0.3), lty = 1)
lines(tvec, Gf(tvec), col = rgb(1, 0, 0, 0.7), lwd = 3)
lines(tvec, fmean, col = rgb(0, 0, 0.55, 0.7), lwd = 3)
lines(tvec, flines[2,], col = rgb(1, 0.55, 0, 0.7), lwd = 3)
lines(tvec, flines[1,], col = rgb(1, 0.65, 0, 0.7), lwd = 2, lty = 2)
lines(tvec, flines[3,], col = rgb(1, 0.65, 0, 0.7), lwd = 2, lty = 2)

#plot(density(f[100,]))
