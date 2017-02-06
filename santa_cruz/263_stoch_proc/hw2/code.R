### No 1
sig2 = 1/2
phi = 1/2
alpha = 0.5

n = 300
ss = seq(-3, 3, length = n)

sigma = sig2*exp(-phi*abs(as.matrix(dist(ss, diag = TRUE, upper = TRUE)))^alpha)
L = t(chol(sigma))

z = rnorm(n)
plot(ss, L %*% z, pch = 16, type='o')



### No 2
library(MASS)
set.seed(1)
n = 1000
x = rnorm(n)
sig = sqrt(0.01)
y = 0.3 + 0.4*x + 0.5*sin(2.7*x) + 1.1/(1+x^2) + rnorm(n, 0, sig)
plot(x, y)

# GP with exponential correlation function
B = 10000
post.sig2 = double(B)+1
post.tau2 = double(B)+1
post.phi = double(B)+1
accept.phi = double(B)
post.mu = double(B)
post.f = matrix(0, B, n)
Dmat = as.matrix(dist(x, diag = TRUE, upper = TRUE))
H.phi = exp(-post.phi[1]*Dmat)
U.phi = chol(H.phi)
U.inv = backsolve(U.phi, diag(n))
f.phi.post = function(phi, f, mu, U.inv)
    sum(log(diag(U.inv)))-1/2*sum((t(f - mu) %*% U.inv)^2)

for (b in 5:B){

    ll = f.phi.post(post.phi[b-1], post.f[b-1,], post.mu[b-1], U.inv)
    cand.phi = rnorm(1, post.phi[b-1], 0.5)
    if (cand.phi > 0 && cand.phi <= 50){
        cand.H.phi = exp(-cand.phi*Dmat)
        cand.U.phi = chol(cand.H.phi)
        cand.U.inv = backsolve(cand.U.phi, diag(n))
        ll.cand = f.phi.post(cand.phi, post.f[b-1,], post.mu[b-1], cand.U.inv)
        if (log(runif(1)) <= ll.cand - ll){
            post.phi[b] = cand.phi
            accept.phi[b] = 1
            H.phi = cand.H.phi
            U.inv = cand.U.inv
        } else {
            post.phi[b] = post.phi[b-1]
            }
    } else {
        post.phi[b] = post.phi[b-1]
        }

    post.sig2[b] = 1/rgamma(1, 1 + n/2, 1 + 1/2*sum((y - post.f[b-1,])^2))
    post.tau2[b] = 1/rgamma(1, 1 + n/2, 1 + 1/2*sum((t(post.f[b-1,] - post.mu[b-1]) %*% U.inv)^2))

    vm = (1/1^2 + sum(t(U.inv) %*% U.inv) / post.tau2[b])^(-1)
    mm = vm * (0/1^2 + sum(t(U.inv) %*% U.inv %*% post.f[b-1,]) / post.tau2[b])
    post.mu[b] = rnorm(1, mm, sqrt(vm))

    A = solve(diag(n)/post.sig2[b] + t(U.inv) %*% U.inv / post.tau2[b])
    M = A %*% (y /  post.sig2[b] + post.mu[b]/post.tau2[b] * (t(U.inv) %*% U.inv) %*% (double(n)+1))
    post.f[b,] = mvrnorm(1, M, A)

    plot(x, post.f[b,])
    points(x, y, col= 'red', pch = 16)

    }


mean(accept.phi)

pp = apply(post.f, 2, mean)
qq = apply(post.f, 2, quantile, c(0.025, 0.5, 0.975))

plot(x, y, pch = 16, ylim = range(qq))
lines(sort(x), pp[order(x)], col = 'blue', lwd = 3)
lines(sort(x), qq[1, order(x)], col = 'blue', lwd = 1)
lines(sort(x), qq[3, order(x)], col = 'blue', lwd = 1)

dim(t(post.f))
length(post.sig2)

params = cbind(post.sig2, post.tau2, post.phi, post.mu, t(post.f))
out = list("params"=params, "accept.phi" = accept.phi, "x" = x, "y" = y)
save(out, file = "./numb2.RData")



m = 100
z = seq(min(x) + 0.01, max(x) - 0.01, length = m)

