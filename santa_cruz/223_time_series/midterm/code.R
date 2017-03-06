### Part 1
dat = read.csv("./googletrendsUCSC.csv")[,2]
times = read.csv("./googletrendsUCSC.csv")[,1]
x.ind = round(seq(1, length(dat), length = 13))
as.numeric(times)

y = diff(dat)
n = length(y)
tt = 1:n

plot(y, type='l', bty = 'n', axes = FALSE, ylab = "Differences",
    xlab = "Time", cex.lab = 1.4, cex.main = 2, main = "First-order difference of search interest")
axis(2)
axis(1, at = x.ind, times[x.ind])
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


### Part 2
get.vals = function(y, delta, m0, C0, n0, d0){
    n = length(y)
    Rt = double(n+1);    Qt = double(n+1);    At = double(n+1)
    ft = double(n+1);    et = double(n+1);    mt = double(n+1)
    nt = double(n+1);    dt = double(n+1);    St = double(n+1)
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
    
    return (list("Rt"=Rt[-1], "Qt"=Qt[-1], "At"=At[-1],
        "ft"=ft[-1], "et"=et[-1], "mt"=mt[-1], "nt"=nt[-1],
        "dt"=dt[-1], "St"=St[-1], "Ct"=Ct[-1], "delta"=delta,
        "m0"=m0, "C0"=C0, "n0"=n0, "d0"=d0))   
    }
obs.pred.dens = function(dat, params){
    df = c(params$n0, params$nt[-length(params$nt)])
    mu = params$ft
    sig2 = params$Qt
    out = lgamma((df+1)/2) - lgamma(df/2) - 1/2*log(df*pi*sig2) -
        (df+1)/2*log(1+(dat - mu)^2 / (df*sig2))
    return (sum(out))
    }
smooth = function(dat, params){
    n = length(dat)
    out.at = double(n)
    out.Rt = double(n)
    out.at[n] = params$mt[n]
    out.Rt[n] = params$Ct[n]
    for (i in (n-1):1){
        out.at[i] = params$mt[i] - params$delta*(params$mt[i] - out.at[i+1])
        out.Rt[i] = params$Ct[i] - params$delta^2*(params$Ct[i] - out.Rt[i+1])
        }
    return (list("at.s"=out.at, "Rt.s"=out.Rt))
    }
forecast = function(dat, params, h = 12){
    n = length(dat)
    delta = params$delta
    at.p = rep(params$mt[n], h)
    Rt.p = double(h+1)
    Rt.p[1] = params$Ct[n]
    for (i in 2:(h-1))
        Rt.p[i] = Rt.p[i-1] + 
    Rt.p

    ft.p = rep(params$mt[n], h)

    C_T + (1 - d)/d * C_T = 1/d C_T
    1/d C_T + 1/d C_T - d/d C_T
    (2- d)/d C_T + 1/d C_T - C_T
    
    }

# Optimal discount factor
n = length(dat)
del = seq(0.501, 1, by = 0.001)
m0 = dat[1]
C0 = var(dat)
n0 = 1
d0 = 100
opd = double(length(del))
for (i in 1:length(del)){
    out = get.vals(dat, delta = del[i], m0 = m0, C0 = C0, n0 = n0, d0 = d0)
    opd[i] = obs.pred.dens(dat, out)
    }

pdf("./figs/discount.pdf", height = 8, width = 8)
par(mfrow = c(1,1), mar = c(5.1, 4.1, 4.1, 2.1))
plot(del, opd, type='l', ylab = "log observed predictive density", bty = 'n',
    xlab = "Discount Factor" ~ delta, cex.lab = 1.4, cex.main = 2)
abline(v = del[which.max(opd)], lty = 3, col = 'gray', lwd = 3)
dev.off()
(d.max = del[which.max(opd)])




out = get.vals(dat, delta = d.max, m0 = m0, C0 = C0, n0 = n0, d0 = d0)
par(mfrow = c(1,1), mar = c(5.1, 4.1, 4.1, 2.1))

plot(dat, bty='n', type='o', ylim = c(0, 120), col = 'gray60',
    xlab = "Time", ylab = "Search interest", cex.lab = 1.4, cex.main = 2,
    axes = FALSE, main = ~ y[t] ~ "|" ~ D[t-1] )
axis(2)
axis(1, at = x.ind, labels = times[x.ind])
lines(out$mt, col = 'blue', lwd = 4)
lines(out$mt + sqrt(out$Qt)*qt(0.975, out$nt-1), col = 'lightblue', lwd = 2)
lines(out$mt - sqrt(out$Qt)*qt(0.975, out$nt-1), col = 'lightblue', lwd = 2)

plot(dat, bty='n', type='o', ylim = c(0, 120), col = 'gray60',
    xlab = "Time", ylab = "Search interest", cex.lab = 1.4, cex.main = 2,
    axes = FALSE, main = ~ theta[t] ~ "|" ~ D[t])
axis(2)
axis(1, at = x.ind, labels = times[x.ind])
lines(out$mt, col = 'green', lwd = 3)
lines(out$mt + sqrt(out$Ct)*qt(0.975, out$nt), col = 'lightgreen', lwd = 2)
lines(out$mt - sqrt(out$Ct)*qt(0.975, out$nt), col = 'lightgreen', lwd = 2)

spar = smooth(dat, out)
plot(dat, bty='n', type='o', ylim = c(0, 120), col = 'gray60',
    cex.lab = 1.4, cex.main = 2, ylab = "Search Interest", axes = FALSE,
    xlab = "Time", main = ~ theta[t] ~ "|" ~ D[T])
axis(2)
axis(1, at = x.ind, labels = times[x.ind])
lines(spar$at.s, col = 'red', lwd = 3)
lines(spar$at.s + sqrt(spar$Rt.s*out$St[n]/out$St)*qt(0.975, out$nt[n]), col = 'pink', lwd = 2)
lines(spar$at.s - sqrt(spar$Rt.s*out$St[n]/out$St)*qt(0.975, out$nt[n]), col = 'pink', lwd = 2)

