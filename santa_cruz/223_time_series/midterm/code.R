### Part 1
dat = read.csv("./googletrendsUCSC.csv")[,2]
times = read.csv("./googletrendsUCSC.csv")[,1]
as.numeric(times)

y = diff(dat)
n = length(y)
tt = 1:n

plot(y, type='l', bty = 'n')
# axis(2)
# ind = round(seq(1, length(y), length = 6))
# axis(1, at = ind, labels = times[ind])

k = 1:floor(n/2-1/2)
omega = 2*pi*k/n

# omega = seq(0.001, pi-0.001, length = 1000)

a.hat = (2/n)*sapply(omega, function(w) sum(y * cos(w * tt)))
b.hat = (2/n)*sapply(omega, function(w) sum(y * sin(w * tt)))
Iw = (n/2)*(a.hat^2 + b.hat^2)
pw = (1 - Iw/sum(y^2))^(1 - n/2)

lambda = 2*pi/omega
pl = pw * (2*pi/(lambda^2))^0

omega[which.max(pw)]
lambda[which.max(pl)]

pdf("./figs/spectral.pdf", height = 8, width = 12)
par(mfrow = c(1,2), mar = c(4.1, 4.1, 2.1, 1.1))
plot(omega, log(pw), type='l', bty='n', main = "Frequency", ylab = "log posterior",
    cex.main = 2.0, cex.lab = 1.5)
plot(lambda, log(pl), type='l', bty='n', main = "Period (in months)", ylab = "log posterior",
    cex.main = 2.0, cex.lab = 1.5, xlab = "lambda", xlim = c(2, 16))
dev.off()


## pp = pw / sum(pw)
## cc = pl / sum(pl)
## bb = sample(omega, 100000, replace = TRUE, prob = pp)
## ll = 2*pi/bb
## 
## 
## #pdf("./sp_lh.pdf", height = 8, width = 13)
## par(mfrow = c(1,2), mar = c(4.1, 4.1, 2.1, 1.1))
## yy = table(bb) / length(bb)
## xx = as.numeric(names(yy))
## plot(xx, as.numeric(yy), xlim = c(0, pi), type='h', bty='n',
##     xlab = "omega", ylab = "posterior", main = "Frequency",
##     cex.main = 2, cex.lab = 1.5)
## lines(omega, pp, col='red')
## 
## yy = table(ll) / length(ll)
## xx = as.numeric(names(yy))
## plot(xx, as.numeric(yy), type='h', bty='n',
##     xlab = "lambda", ylab = "posterior", main = "Period (in 10 minutes)",
##     cex.main = 2, cex.lab = 1.5, bty='n')
## lines(lambda, cc, col='red')
## #dev.off()

### Part 2

get.vals = function(y, delta, m0, C0, n0, d0){
    n = length(y)
    Rt = double(n+1)
    Qt = double(n+1)
    At = double(n+1)
    ft = double(n+1)
    et = double(n+1)
    mt = double(n+1)
    nt = double(n+1)
    dt = double(n+1)
    St = double(n+1)
    Ct = double(n+1)

    mt[1] = m0
    Ct[1] = C0
    nt[1] = n0
    dt[1] = d0

    for (i in 2:(n+1)){
        Rt[i] = Ct[i-1] / delta # Ct[i-1] + Wt[i]
        Qt[i] = Rt[i] + St[i-1]
        At[i] = Rt[i]/Qt[i]   
        ft[i] = mt[i-1]
        et[i] = y[i-1] - ft[i]
        mt[i] = mt[i-1] + At[i]*et[i]
        nt[i] = nt[i-1] + 1
        dt[i] = dt[i-1] + St[i-1]*(et[i]^2) / Qt[i]
        St[i] = dt[i] / nt[i]
        Ct[i] = At[i] * St[i]
        }
    
    return (list("Rt"=Rt, "Qt"=Qt, "At"=At,
        "ft"=ft, "et"=et, "mt"=mt, "nt"=nt,
        "dt"=dt, "St"=St, "Ct"=Ct))   
    }
obs.pred.dens = function(dat, params){
    df = params$nt[-length(params$nt)]
    mu = params$ft[-1]
    sig2 = params$Qt[-1]
    out = lgamma((df+1)/2) - lgamma(df/2) - 1/2*log(df*pi*sig2) -
        (df+1)/2*log(1+(dat - mu)^2 / (df*sig2))
    return (sum(out))
    }

par(mfrow = c(1,1), mar = c(5.1, 4.1, 4.1, 2.1))
plot(dat, bty='n', type='l')
out = get.vals(dat, delta = 0.8, m0 = mean(dat), C0 = var(dat), n0 = 5, d0 = 100)


plot(dat, bty='n', type='l', ylim = c(0, 100))
lines(out$mt[-1], col = 'blue', lwd = 3)
lines(out$mt[-1] + sqrt(out$Qt[-1])*qt(0.975, out$nt[-(n+1)]), col = 'green')
lines(out$mt[-1] - sqrt(out$Qt[-1])*qt(0.975, out$nt[-(n+1)]), col = 'green')

opd = double(11)
del = seq(0.5, 1, length = 11)
for (i in 1:length(del)){
    out = get.vals(dat, delta = del[i], m0 = mean(dat), C0 = var(dat), n0 = 5, d0 = 100)
    opd[i] = obs.pred.dens(dat, out)
    }
