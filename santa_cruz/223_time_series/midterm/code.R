### Part 1
dat = read.csv("./googletrendsUCSC.csv")[,2]
times = read.csv("./googletrendsUCSC.csv")[,1]
x.ind = round(seq(1, length(dat), length = 13))
as.numeric(times)

y = diff(dat)
n = length(y)
tt = 1:n

# Plot of data
pdf("./figs/data.pdf", width = 14, height = 8)
par(mfrow = c(1,2), mar = c(5.1, 4.1, 4.1, 2.1))
plot(dat, type='o', bty = 'n', axes = FALSE, ylab = "Search Interest",
    xlab = "Time", cex.lab = 1.4, cex.main = 2, main = "Original Data",
    col = 'gray60')
axis(2)
axis(1, at = x.ind, times[x.ind])
plot(y, type='o', bty = 'n', axes = FALSE, ylab = "Differences",
    xlab = "Time", cex.lab = 1.4, cex.main = 2, col = 'gray60',
    main = "First-order differences",
    )
axis(2)
axis(1, at = x.ind, times[x.ind])
dev.off()



# Spectrum
k = 1:floor(n/2-1/2)
omega = 2*pi*k/n

a.hat = (2/n)*sapply(omega, function(w) sum(y * cos(w * tt)))
b.hat = (2/n)*sapply(omega, function(w) sum(y * sin(w * tt)))
Iw = (n/2)*(a.hat^2 + b.hat^2)
pw = (1 - Iw/sum(y^2))^(1 - n/2)

lambda = 2*pi/omega
pl = pw * (2*pi/(lambda^2))^0

omega[which.max(pw)]
lambda[which.max(pl)]

pdf("./figs/spectral.pdf", width = 14, height = 8)
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
    Rt.p = (1/delta^(1:h))*params$Ct[n]
    ft.p = rep(params$mt[n], h)
    Qt.p = Rt.p + params$St[n]
    return (list("at.p"=at.p, "Rt.p"=Rt.p,
        "ft.p"=ft.p, "Qt.p"=Qt.p))
    }

# Optimal discount factor
n = length(dat)
del = seq(0.500, 1, by = 0.001)
m0 = 60
C0 = 200
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
abline(v = del[which.max(opd)], lty = 2, col = 'gray', lwd = 3)
dev.off()
(d.max = del[which.max(opd)])




out = get.vals(dat, delta = d.max, m0 = m0, C0 = C0, n0 = n0, d0 = d0)
pdf("./figs/current.pdf", width = 14, height = 8)
par(mfrow = c(1,2), mar = c(5.1, 4.1, 4.1, 2.1))
plot(dat, bty='n', type='o', ylim = c(0, 120), col = 'gray60',
    xlab = "Time", ylab = "Search interest", cex.lab = 1.4, cex.main = 2,
    axes = FALSE, main = ~ theta[t] ~ "|" ~ D[t])
axis(2)
axis(1, at = x.ind, labels = times[x.ind])
lines(out$mt, col = 'green', lwd = 3)
lines(out$mt + sqrt(out$Ct)*qt(0.975, out$nt), col = 'lightgreen', lwd = 2)
lines(out$mt - sqrt(out$Ct)*qt(0.975, out$nt), col = 'lightgreen', lwd = 2)

plot(dat, bty='n', type='o', ylim = c(0, 120), col = 'gray60',
    xlab = "Time", ylab = "Search interest", cex.lab = 1.4, cex.main = 2,
    axes = FALSE, main = ~ y[t] ~ "|" ~ D[t-1] )
axis(2)
axis(1, at = x.ind, labels = times[x.ind])
lines(out$mt, col = 'blue', lwd = 4)
lines(out$mt + sqrt(out$Qt)*qt(0.975, out$nt-1), col = 'lightblue', lwd = 2)
lines(out$mt - sqrt(out$Qt)*qt(0.975, out$nt-1), col = 'lightblue', lwd = 2)
dev.off()

spar = smooth(dat, out)
pdf("./figs/smooth.pdf", width = 8, height = 8)
par(mfrow = c(1,1), mar = c(5.1, 4.1, 4.1, 2.1))
plot(dat, bty='n', type='o', ylim = c(0, 120), col = 'gray60',
    cex.lab = 1.4, cex.main = 2, ylab = "Search Interest", axes = FALSE,
    xlab = "Time", main = ~ theta[t] ~ "|" ~ D[T])
axis(2)
axis(1, at = x.ind, labels = times[x.ind])
lines(spar$at.s, col = 'red', lwd = 3)
lines(spar$at.s + sqrt(spar$Rt.s*out$St[n]/out$St)*qt(0.975, out$nt[n]), col = 'pink', lwd = 2)
lines(spar$at.s - sqrt(spar$Rt.s*out$St[n]/out$St)*qt(0.975, out$nt[n]), col = 'pink', lwd = 2)
dev.off()


h = 12
ppar = forecast(dat, out, h)
pdf("./figs/forecast.pdf", width = 14, height = 8)
par(mfrow = c(1,2), mar = c(5.1, 4.1, 4.1, 2.1))
plot(1:(n+h), c(dat, rep(NA, h)), bty='n', type='o', ylim = c(0, 120), col = 'gray60',
    cex.lab = 1.4, cex.main = 2, ylab = "Search Interest", axes = FALSE,
    xlab = "Time", main = ~ theta[T+h] ~ "|" ~ D[T])
axis(2)
axis(1, at = x.ind, labels = times[x.ind])
lines(out$mt, col = 'green', lwd = 3)
lines(out$mt + sqrt(out$Ct)*qt(0.975, out$nt), col = 'lightgreen', lwd = 2)
lines(out$mt - sqrt(out$Ct)*qt(0.975, out$nt), col = 'lightgreen', lwd = 2)
lines((n+1):(n+h), ppar$at.p, col = 'darkgreen', lwd = 3)
lines((n+1):(n+h), ppar$at.p + sqrt(ppar$Rt.p)*qt(0.975, out$nt[n]),
    col = 'darkgreen', lwd = 2)
lines((n+1):(n+h), ppar$at.p - sqrt(ppar$Rt.p)*qt(0.975, out$nt[n]),
    col = 'darkgreen', lwd = 2)

plot(1:(n+h), c(dat, rep(NA, h)), bty='n', type='o', ylim = c(0, 120), col = 'gray60',
    xlab = "Time", ylab = "Search interest", cex.lab = 1.4, cex.main = 2,
    axes = FALSE, main = ~ y[T+h] ~ "|" ~ D[T] )
axis(2)
axis(1, at = x.ind, labels = times[x.ind])
lines(out$mt, col = 'blue', lwd = 4)
lines(out$mt + sqrt(out$Qt)*qt(0.975, out$nt-1), col = 'lightblue', lwd = 2)
lines(out$mt - sqrt(out$Qt)*qt(0.975, out$nt-1), col = 'lightblue', lwd = 2)
lines((n+1):(n+h), ppar$ft.p, col = 'darkblue', lwd = 3)
lines((n+1):(n+h), ppar$ft.p + sqrt(ppar$Qt.p)*qt(0.975, out$nt[n]),
    col = 'darkblue', lwd = 2)
lines((n+1):(n+h), ppar$at.p - sqrt(ppar$Qt.p)*qt(0.975, out$nt[n]),
    col = 'darkblue', lwd = 2)
dev.off()

