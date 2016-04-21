### Model 1: assume all points come from the same distribution
library(MASS)
rwishart = function(V, df){
    p = nrow(V)
    U = chol(V)
    A = diag(sqrt(rchisq(p, df-(1:p)+1)))
    A[upper.tri(A)] = rnorm(p*(p-1)/2)
    X = A %*% U
    return (t(X) %*% X)
    }
riwishart = function(V, df){
    p = nrow(V)
#   U = backsolve(chol(V), diag(p))
    U = chol(V)
    A = diag(sqrt(rchisq(p, df-(1:p)+1)))
    A[upper.tri(A)] = rnorm(p*(p-1)/2)
    A = backsolve(A, diag(p))
    X = U %*% A
    return (X %*% t(X))
    }


dat = read.csv("~/files/data/CAtheft.csv")

y = dat
#for (i in 1:nrow(y))
#    y[i,3:6] = y[i,3:6] / y[i,2]
z = log(y[,3:6])

n = nrow(z)
p = ncol(z)

pairs(z, pch = 20)
plot(y[,2], z[,4], pch = 20)

zbar = apply(z, 2, mean)
#S1 = matrix(0, p, p)
#for (i in 1:n)
#    S1 = S1 + matrix(as.numeric(z[i,] - zbar), 4, 1) %*%
#        matrix(as.numeric(z[i,] - zbar), 1, 4) 

## m = zbar
## C0 = diag(p)
## r = 5
## S0 = diag(p)
## 
## nburn = 100
## nmcmc = 1000
## params.mu = matrix(0, nburn + nmcmc, p)
## params.sigma = rep(list(matrix(0, p, p)), nburn + nmcmc)
## 
## params.mu[1,] = zbar
## 
## for (i in 2:(nburn + nmcmc)){
##     if (floor(i/100) == i/100)
##         cat(i, "/", nburn+nmcmc, "\r")
##     S1 = matrix(0, p, p)
##     for (j in 1:n)
##         S1 = S1 + matrix(as.numeric(z[j,] - params.mu[i-1,]), 4, 1) %*%
##             matrix(as.numeric(z[j,] - params.mu[i-1,]), 1, 4) 
##     params.sigma[[i]] = riwishart(solve(S0 + S1), r + n)
##     Vinv = solve(params.sigma[[i]])
##     C1 = solve(solve(C0) + n*Vinv)
##     m1 = C1 %*% (solve(C0) %*% m + n * Vinv %*% zbar)
##     params.mu[i,] = mvrnorm(1, m1, C1)
##     if (i == nburn + nmcmc)
##         cat("\n")
##     }
## 
## params.mu = tail(params.mu, nmcmc)
## params.sigma = tail(params.sigma, nmcmc)
## 
## apply(params.mu, 2, mean)
## apply(params.mu, 2, quantile, c(0.025, 0.5, 0.975))
## 
## pred.z = matrix(0, nmcmc, p)
## for (i in 1:nmcmc)
##     pred.z[i,] = mvrnorm(1, params.mu[i,], params.sigma[[i]])
## pairs(rbind(pred.z, as.matrix(z)), col = c(rep(2, nmcmc), rep(1, n)), pch = 20)


### Model 2: each point comes from it's own distribution
m = zbar
C0 = 1*diag(p)
C0inv = solve(C0)
r = 6
S0 = 0.25*diag(p)
k = 5
D0 = 3*diag(p)

nburn = 1000
nmcmc = 2000
params.mu.i = rep(list(matrix(0, nburn + nmcmc, p)), n)
params.mu = matrix(0, nburn + nmcmc, p)
params.sigma = rep(list(matrix(0, p, p)), nburn + nmcmc)
params.V = rep(list(matrix(0, p, p)), nburn + nmcmc)

params.mu[1,] = zbar
for (i in 1:n)
    params.mu.i[[i]][i,] = as.matrix(z[i,])

calc.post = function(mu.i, mu, sigma, V, j){
    out = -(n+p+r+1)/2*determinant(sigma[[j]])$modulus[1]
    out = out -(n+p+k+1)/2*determinant(V[[j]])$modulus[1]
    Sinv = solve(sigma[[j]])
    Vinv = solve(V[[j]])
    for (i in 1:n){
        out = out - 1/2*t(as.numeric(z[i,]) - mu.i[[i]][j,]) %*%
            Sinv %*% (as.numeric(z[i,]) - mu.i[[i]][j,])
        out = out - 1/2*t(mu.i[[i]][j,] - mu[j,]) %*%
            Vinv %*% (mu.i[[i]][j,] - mu[j,])
        }
    out = out - 1/2*t(mu[j,] - m) %*% solve(C0) %*% (mu[j,] - m)
    out = out - 1/2*sum(diag(S0 %*% Sinv))
    out = out - 1/2*sum(diag(D0 %*% Vinv))
    return (out)
    }


for (i in 2:(nburn + nmcmc)){
    if (floor(i/100) == i/100)
        cat(i, "/", nburn+nmcmc, "\r")
    # Update sigma
    S1 = matrix(0, p, p)
    for (j in 1:n)
        S1 = S1 + matrix(as.numeric(z[j,] - params.mu.i[[j]][i-1,]), 4, 1) %*%
            matrix(as.numeric(z[j,] - params.mu.i[[j]][i-1,]), 1, 4) 
#   params.sigma[[i]] = riwishart(solve(S0 + S1), r + n)
    params.sigma[[i]] = solve(rwishart(solve(S0 + S1), r + n))

    # Update V
    D1 = matrix(0, p, p)
    for (j in 1:n)
        D1 = D1 + matrix(params.mu.i[[j]][i-1,] - params.mu[i-1,], 4, 1) %*%
            matrix(params.mu.i[[j]][i-1,] - params.mu[i-1,], 1, 4) 
#   params.V[[i]] = riwishart(solve(D0 + D1), r + n)
    params.V[[i]] = solve(rwishart(solve(D0 + D1), r + n))
    rwishart(solve(D0 + D1), r + n)

    # Update mu
    Siginv = solve(params.sigma[[i]])
    Vinv = solve(params.V[[i]])
    mubar = apply(sapply(params.mu.i, function(x) x[i-1,]), 1, mean)
    C1 = solve(C0inv + n*Vinv)
    m1 = C1 %*% (C0inv %*% m + n * Vinv %*% mubar)
    params.mu[i,] = mvrnorm(1, m1, C1)

    # Update mu_i
    A1 = solve(Siginv + Vinv)
    for (j in 1:n){
        b1 = A1 %*% (Siginv %*% as.numeric(z[j,]) + Vinv %*% params.mu[i,])
        params.mu.i[[j]][i,] = mvrnorm(1, b1, A1)
        }

    if (i == nburn + nmcmc)
        cat("\n")
    }


# Burn-in
params.mu = tail(params.mu, nmcmc)
for (j in 1:n)
    params.mu.i[[j]] = tail(params.mu.i[[j]], nmcmc)
params.sigma = tail(params.sigma, nmcmc)
params.V = tail(params.V, nmcmc)

# Plot of the log-posterior
post = double(nmcmc)
for (i in 1:nmcmc)
    post[i] = calc.post(params.mu.i, params.mu, params.sigma, params.V, i)
plot(post, type='l')

par(mfrow = c(2,2), mar=c(4.1, 2.1, 2.1, 1.1))
for (j in 1:p)
    plot(params.mu[,j], type='l')
#plot(params.mu.i[[2]][,2], type='l')
par(mfrow = c(p,p), mar=c(4.1, 2.1, 2.1, 1.1))
for (i in 1:p)
    for (j in 1:p)
        plot(sapply(params.sigma, function(x) x[i,j]), type='l')
for (i in 1:p)
    for (j in 1:p)
        plot(sapply(params.V, function(x) x[i,j]), type='l')
par(mfrow = c(1,1), mar=c(5.1, 4.1, 4.1, 2.1))

pairs(rbind(as.matrix(params.mu.i[[31]]), as.matrix(params.mu.i[[32]]),
    as.matrix(z[31,]), as.matrix(z[32,])),
    col=c(rep(rgb(0.3,0,0), nmcmc), rep(rgb(0,0,0.3), nmcmc), rgb(0.9,0,0), rgb(0,0,0.9)), pch = 20)

pairs(params.mu, pch = 20)


### Make predictions for a specific county (not for a new county)
j = 31
pred.z = matrix(0, nmcmc, p)
for (i in 1:nmcmc)
    pred.z[i,] = mvrnorm(1, params.mu.i[[j]][i,], params.sigma[[i]])
pairs(rbind(pred.z, as.matrix(z[j,])), col = c(rep(2, nmcmc), rep(1, 1)), pch = 20)
mean(apply(exp(pred.z), 1, sum) >= 3000)

mean.preds = matrix(0, n, p)
qq.preds = rep(list(matrix(0, 2, 4)), n)
for (j in 1:n){
    pred.z = matrix(0, nmcmc, p)
    for (i in 1:nmcmc)
        pred.z[i,] = mvrnorm(1, params.mu.i[[j]][i,], params.sigma[[i]])
    mean.preds[j,] = apply(pred.z, 2, mean)
    qq.preds[[j]] = apply(pred.z, 2, quantile, c(0.025, 0.975))
    }

j = 23
x = 1+double(n)
x[j] = 2
pairs(mean.preds, pch = 20, col = x)
pairs(y[,3:6], pch = 20, col = col.vec)

#col.vec = rgb((y[,2] - min(y[,2])) / diff(range(y[,2])),
#    0, 1-(y[,2] - min(y[,2])) / diff(range(y[,2])))
col.vec = rgb((rank(y[,2]) - 1) / (n-1),
    0, 1-(rank(y[,2]) - 1) / (n-1))
par(mfrow = c(2,2), mar = c(5.1, 4.1, 2.1, 1.1))
for (j in 1:p){
    temp = sapply(qq.preds, function(x) x[,j])
    plot(z[,j], mean.preds[,j], pch = 20, ylim = range(temp), main = names(z)[j],
        xlab = "Observed log(Theft)", ylab = "Predicted log(Theft)")
    abline(0, 1)
    segments(x0 = z[,j], y0 = temp[1,], y1 = temp[2,], col = col.vec)
    }
par(mfrow = c(1,1), mar = c(5.1, 4.1, 4.1, 2.1))


apply(exp(pred.z), 2, mean) - y[j,3:6]
#exp(apply(pred.z, 2, mean)) - y[j,3:6]

