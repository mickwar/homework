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

tmp = read.table("~/files/repos/data/soi.txt")[,1]
tmp = diff(c(co2))

dat = NULL
dat$tt = seq(1, length(tmp), by = 1)
dat$y = tmp

plot(dat$tt, dat$y, type='l', bty='n')

calc.post = function(omega, dat){
    tt = dat$tt
    y = dat$y
    n = length(y)
    Fm = rbind(cos(omega * tt), sin(omega * tt)) # = F'
    beta.hat = solve(Fm %*% t(Fm)) %*% Fm %*% y

#   exp((-0.5)*determinant(Fm %*% t(Fm))$modulus[1]) *
#       (1 - (t(beta.hat) %*% Fm %*% t(Fm) %*% beta.hat) / sum(y^2))^((2 - n)/2)

    (-0.5)*determinant(Fm %*% t(Fm))$modulus[1] +
        ((2-n)/2)*log(1 - (t(beta.hat) %*% Fm %*% t(Fm) %*% beta.hat) / sum(y^2))

    }

omega = seq(0.001, pi - 0.001, length = 1000)
lambda = 2*pi / omega

tp = sapply(omega, function(x) calc.post(x, dat))

tp = sapply(lambda, function(x) calc.post(2*pi/x, dat)-2*log(x))

plot(omega, tp, type='l')

plot(lambda, tp, type='l', xlim = c(2, 100))

calc.post(0.001, dat)
