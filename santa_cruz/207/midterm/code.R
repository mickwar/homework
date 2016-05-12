dat = read.table("~/volcano.txt", header=TRUE)
y = log(dat[1:62,3])
n = length(y)

plot(density(y))

# Hyperpriors
m = 0   # mu mean
s2 = 10 # mu variance
a = 3   # sig2 alpha
b = 3   # sig2 beta
c = 3   # tau2 alpha
d = 3   # tau2 beta
df = 5  # df for the t, fixed

nburn = 10000
nmcmc = 50000
window = 500

params.lambda.i = matrix(0, nburn + nmcmc, n)
params.mu.i = matrix(0, nburn + nmcmc, n)
params.mu = double(nburn + nmcmc)
params.tau2 = double(nburn + nmcmc)
params.sig2 = double(nburn + nmcmc)

params.tau2[1] = 1
params.sig2[1] = 1

for (i in 2:(nburn + nmcmc)){
    if (floor(i / window) == i/window)
        cat(i, "/", nburn+nmcmc, "\r")
    # Update lambda's
    params.lambda.i[i,] = rgamma(n, (df + 1)/2, df/2 + 1/(2*params.sig2[i-1])*
        (y - params.mu.i[i-1,])^2)

    # Update mu.i's
    params.mu.i[i,] = rnorm(n,
        (params.tau2[i-1]*y*params.lambda.i[i,] + params.sig2[i-1]*params.mu[i-1]) / 
            (params.tau2[i-1]*params.lambda.i[i,] + params.sig2[i-1]),
        sqrt((params.tau2[i-1]*params.sig2[i-1]) /
            (params.tau2[i-1]*params.lambda.i[i,] + params.sig2[i-1])))

    # Update mu
    params.mu[i] = rnorm(1, 
        (n*s2*mean(params.mu.i[i,]) + params.tau2[i-1]*m) /
            (n*s2+params.tau2[i-1]),
        sqrt((s2 * params.tau2[i-1]) /
            (n*s2+params.tau2[i-1])))

    # Update tau2
    params.tau2[i] = 1/rgamma(1, c + n/2, d + 1/2*sum((params.mu.i[i,] - params.mu[i])^2))

    # Update sig2
    params.sig2[i] = 1/rgamma(1, a + n/2, b + 1/2*sum(params.lambda.i[i,] *
        (y - params.mu.i[i,])^2))

    if (i == (nburn + nmcmc))
        cat("\n")
    }

params.lambda.i = tail(params.lambda.i, nmcmc)
params.mu.i = tail(params.mu.i, nmcmc)
params.mu = tail(params.mu, nmcmc)
params.tau2 = tail(params.tau2, nmcmc)
params.sig2 = tail(params.sig2, nmcmc)

plot(params.sig2, type='l')
plot(params.tau2, type='l')
plot(params.mu, type='l')
#matplot(tail(params.mu.i, 1000), type='l')
#matplot(tail(params.lambda.i, 1000), type='l')

### Posterior predictive distribution (for each observation)
pred.y = matrix(0, nmcmc, n)
for (j in 1:n)
    pred.y[,j] = rnorm(nmcmc, params.mu.i[,j], sqrt(params.sig2 / params.lambda.i[,j]))
qq = apply(pred.y, 2, quantile, c(0.025, 0.975))

plot(y, apply(pred.y, 2, mean), pch = 20, ylim = range(qq))
segments(x0 = y, x1 = y, y0 = qq[1,], y1 = qq[2,])
abline(0, 1)

#f = function(x)
#    dt(x-7, 5)
#curve(f(x), col = 'blue', add=  TRUE)

#plot(density(params.sig2))
#plot(density(params.tau2))


### Posterior predictive for a new observation
pred.mu.0 = rnorm(nmcmc, params.mu, sqrt(params.tau2))
pred.lambda.0 = rgamma(nmcmc, df/2, df/2)
pred.y.0 = rnorm(nmcmc, pred.mu.0, sqrt(params.sig2 / pred.lambda.0))
plot(density(y))
lines(density(pred.y.0), col = 'darkgreen', lwd = 3)
for (j in 1:n){
    dens = density(pred.y[,j])
    lines(dens$x, 0.1*dens$y, col = rgb(0, 1, 0, 0.25))
    }
points(y, rep(0, n), pch = 20)

#lines(density(pred.mu.0), col = 'blue',lwd = 3)

#apply(params.lambda.i, 2, range)

#plot(density(pred.lambda.0), lwd = 2)
#for (j in 1:n)
#    lines(density(params.lambda.i[,j]), col = rgb(0, 1, 0, 0.25))
#lines(density(c(params.lambda.i)), col = 'darkgreen', lwd = 2)
#
