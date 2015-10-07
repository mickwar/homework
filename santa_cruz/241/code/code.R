# rdirichlet2 = function(n, alpha){
#     k = length(alpha)
# 
#     # get all the chisq draws
# #   x = matrix(rchisq(n*k, 2*alpha), n, k, byrow=TRUE)
# 
#     # alternate form (the scale parameter doesn't matter, each
#     # draw just needs to have the same scale)
#     x = matrix(rgamma(n*k, alpha, 1), n, k, byrow=TRUE)
# 
#     # normalize
#     x / apply(x, 1, sum)
#     }

cdf = function(x, t){
    out = length(t)
    for (i in 1:length(t))
        out[i] = mean(x <= t[i])
    return (out)
    }

#plot(rdirichlet(1000, c(1, 8, 1))[,2:3], pch = 20, xlim = c(0, 1), ylim = c(0,1))

par(mfrow = c(2,1), mar = c(4.1, 3.1, 3.1, 1.1))

N = 1000
library(MCMCpack)
alpha = 5
#x = seq(0.000, 1.000, by = 0.005)
#k = length(x)
#dG0 = punif(x[-1], 0, 1) - punif(x[-k], 0, 1)
#dG = cbind(0, rdirichlet(N, alpha * dG0))
#G = apply(dG, 1, cumsum)
#Gmean = apply(G, 1, mean)
#
#matplot(x, G, type = 'l', col = rgb(0.5, 0.5, 0.5, 0.3), lty = 1)
#curve(punif(x, 0, 1), add = TRUE, lwd = 3)
#lines(x, Gmean, col = 'green', lty = 2, lwd = 3)
#matplot(x, G[,1:5], type = 'l', col = rgb(0.5, 0.0, 0.5), lty = 1, add = TRUE)



x = seq(12, 42, length = 100)
k = length(x)
dG0 = pnorm(x, 29, 5) - pnorm(c(-Inf, x[-k]), 29, 5)
dG = rdirichlet(N, alpha * dG0)
G = apply(dG, 1, cumsum)
Gmean = apply(G, 1, mean)

matplot(x, G, type = 'l', col = rgb(0.5, 0.5, 0.5, 0.3), lty = 1)
curve(pnorm(x, 29, 5), add = TRUE, lwd = 3)
lines(x, Gmean, col = 'green', lty = 2, lwd = 3)
matplot(x, G[,1:5], type = 'l', col = rgb(0.5, 0.0, 0.5), lty = 1, add = TRUE)


#n = 15

dat = read.table("~/files/data/MPG_saturn.txt", header = TRUE)
y = dat$miles / dat$gallons
n = length(y)
#y = rnorm(n, 3, 2)

x = seq(12, 42, length = 500)
k = length(x)
dG0 = pnorm(x, 29, 5) - pnorm(c(-Inf, x[-k]), 29, 5)
dFn = cdf(y, x) - cdf(y, c(-Inf, x[-k]))
dG0 = (alpha*dG0 + n*dFn) / (alpha + n)
dG = rdirichlet(N, (alpha+n) * dG0)
G = apply(dG, 1, cumsum)
Gmean = apply(G, 1, mean)

matplot(x, G, type = 'l', col = rgb(0.5, 0.5, 0.5, 0.3), lty = 1)
curve(pnorm(x, 29, 5), add = TRUE, lwd = 3)
lines(x, Gmean, col = 'green', lty = 2, lwd = 3)
matplot(x, G[,1:5], type = 'l', col = rgb(0.5, 0.0, 0.5), lty = 1, add = TRUE)
points(x, cdf(y, x), col = 'red', pch = 20, cex = 0.5)
#plot(ecdf(y), add = TRUE)


### Polya urn
# N = 10000
# alpha = c(15, 2, 12)
# p = matrix(0, N, length(alpha))
# p[1,] = alpha
# for (i in 2:N){
#     j = sample(length(alpha), 1, prob = p[i-1,])
#     p[i,] = p[i-1,]
#     p[i,j] = p[i,j] + 1
#     }
# 
# #alpha / sum(alpha)
# p[N,] / (N + sum(alpha) - 1) # as N -> Inf, this becomes a draw from Dirichlet(alpha)
# rdirichlet(1, alpha)
