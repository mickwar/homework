# Monthly values of the The Southern Oscillation Index (SOI) during 
# 1950-1995. This series consists of 540 observations on the SOI, 
# computed as the ``difference of the departure from the long-term 
# monthly mean sea level pressures'' at Tahiti in the South Pacific 
# and Darwin in Northern Australia. The index is one measure of the 
# so-called "El Nino-Southern Oscillation" -- an event of critical 
# importance and interest in climatological studies in recent decades. 
# The fact that most of the observations in the last part of the series 
# take negative values is related to a recent warming in the tropical 
# Pacific, and this pattern has been of key interest in recent studies. 
# A key question of interest is to determine just how unusual this event 
# is, and if it can reasonably be explained by standard "stationary" 
# time series models, or requires models that include drifts/trends 
# that may be related to global climatic change. 

library(MASS)

# tmp = read.table("~/files/repos/data/soi.txt")[,1]
# tmp = lh
tmp = diff(c(co2))

dat = NULL
dat$tt = seq(1, length(tmp), by = 1)
dat$y = tmp

y = dat$y
tt = dat$tt
n = length(y)

plot(dat$tt, dat$y, type='l', bty='n')

calc.post = function(omega, param){
    beta = param[1:2]
    v = param[3]

    Fm = cbind(cos(omega * tt), sin(omega * tt)) # = F'
    exp(sum(dnorm(y, Fm %*% beta, v, log = TRUE)))

#   (-0.5)*determinant(Fm %*% t(Fm))$modulus[1] +
#       ((2-n)/2)*log(1 - sum((t(beta.hat) %*% Fm)^2) / sum(y^2))

    }

nmcmc = 10000
nburn = 10000

post.omega = double(nburn + nmcmc)
post.beta = matrix(0, nburn + nmcmc, 2)
post.v = double(nburn + nmcmc)
omega.grid = 2*pi*(1:floor(n/2 - 0.5))/n

post.omega[1] = omega.grid[39]
post.v[1] = 1

for (i in 2:(nburn + nmcmc)){
    if (floor(i/1000) == i/1000)
        cat(i, "/", nburn+nmcmc, "\r")

    Fm = rbind(cos(post.omega[i-1] * tt), sin(post.omega[i-1] * tt))
    post.beta[i,] = mvrnorm(1, solve(Fm %*% t(Fm)) %*% Fm %*% y, post.v[i-1]*solve(Fm %*% t(Fm)))
    post.v[i] = 1/rgamma(1, n/2, 1/2*(sum((y - t(Fm) %*% post.beta[i,])^2)))

    pp = sapply(omega.grid, function(x) calc.post(x, c(post.beta[i,], post.v[i])))
    post.omega[i] = sample(omega.grid, 1, prob = pp)

    }

plot(head(post.v, i-1))
plot(head(post.beta[,2], i-1))

post.omega = tail(post.omega, nmcmc)
post.v = tail(post.v, nmcmc)
post.beta = tail(post.beta, nmcmc)

plot(density(post.omega))
plot(density(post.v))
plot(density(post.beta[,1]))
plot(density(post.beta[,2]))

plot(density(2*pi/post.omega))


plot(post.omega, type='l')
plot(post.v, type='l')
plot(post.beta[,1], type='l')
plot(post.beta[,2], type='l')

plot(post.omega, post.beta[,1], pch = 16)
plot(post.omega, post.beta[,2], pch = 16)

omega = seq(0.001, pi - 0.001, length = 1000)
lambda = 2*pi / omega

tp = sapply(omega, function(x) calc.post(x, dat))

# tp = sapply(lambda, function(x) calc.post(2*pi/x, dat)-2*log(x))

plot(omega, tp, type='l')

# plot(lambda, tp, type='l', xlim = c(2, 100))

calc.post(0.001, dat)
