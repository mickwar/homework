tmp = read.table("~/files/repos/data/soi.txt")[,1]
#tmp = lh
#tmp = diff(c(co2))

dat = NULL
dat$tt = seq(1, length(tmp), by = 1)
dat$y = tmp

y = dat$y
tt = dat$tt
n = length(y)

plot(dat$tt, dat$y, type='l', bty='n')

k = 1:floor(n/2-1/2)

omega = expand.grid(2*pi*k/n, 2*pi*k/n)

omega.1 = 2*pi*k/n
omega.2 = 2*pi*l/n

a.hat.1 = (2/n)*sapply(omega[,1], function(w) sum(y * cos(w * tt)))
b.hat.1 = (2/n)*sapply(omega[,1], function(w) sum(y * sin(w * tt)))
a.hat.2 = (2/n)*sapply(omega[,2], function(w) sum(y * cos(w * tt)))
b.hat.2 = (2/n)*sapply(omega[,2], function(w) sum(y * sin(w * tt)))

Iw = (n/2)*(a.hat.1^2 + b.hat.1^2 + a.hat.2^2 +b.hat.2^2)

plot(omega[,1], (1 - Iw/sum(y^2))^(2 - n/2))

library(rgl)
plot3d(omega[,1], omega[,2], log(Iw, 10))
plot3d(omega[,1], omega[,2], log((1 - Iw/sum(y^2))^(2 - n/2), 10))

calc.post = function(omega, param){
    beta = param[1:2]
    v = param[3]

    Fm = cbind(cos(omega * tt), sin(omega * tt)) # = F'
    exp(sum(dnorm(y, Fm %*% beta, v, log = TRUE)))

#   (-0.5)*determinant(Fm %*% t(Fm))$modulus[1] +
#       ((2-n)/2)*log(1 - sum((t(beta.hat) %*% Fm)^2) / sum(y^2))

