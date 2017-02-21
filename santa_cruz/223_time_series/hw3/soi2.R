tmp = read.table("~/files/repos/data/soi.txt")[,1]
#tmp = lh
#tmp = diff(c(co2))

#tmp = rnorm(4000)

#tmp = sin(seq(0, 1, length=1000)*2*pi * 4*10) +
#    sin(seq(0, 1, length=1000)*2*pi * 7*10) +
#    sin(seq(0, 1, length=1000)*2*pi * 18*10) +
#    sin(seq(0, 1, length=1000)*2*pi * 29*10)

dat = NULL
dat$tt = seq(1, length(tmp), by = 1)
dat$y = tmp

y = dat$y
tt = dat$tt
n = length(y)

# pdf("./dat_soi.pdf", width = 8, height = 8)
# par(mfrow = c(1,1), mar = c(5.1, 4.1, 4.1, 2.1))
# plot(dat$tt / 12, dat$y, type='l', bty='n', cex.main = 2, cex.lab = 1.5,
#     main = "SOI", ylab = "Index", xlab = "Time (in yearss)")
# dev.off()

k = 1:floor(n/2-1/2)
omega = 2*pi*k/n

# omega = seq(0.001, pi-0.001, length = 1000)

a.hat = (2/n)*sapply(omega, function(w) sum(y * cos(w * tt)))
b.hat = (2/n)*sapply(omega, function(w) sum(y * sin(w * tt)))

Iw = (n/2)*(a.hat^2 + b.hat^2)
pw = (1 - Iw/sum(y^2))^(1 - n/2)


lambda = 2*pi/omega
pl = pw * 2*pi/(lambda^2)


omega[which.max(pw)]
lambda[which.max(pl)]

#pdf("./fp_lh.pdf", height = 8, width = 13)
par(mfrow = c(1,2), mar = c(4.1, 4.1, 2.1, 1.1))
plot(omega, log(pw), type='l', bty='n', main = "Frequency", ylab = "log posterior",
    cex.main = 2.0, cex.lab = 1.5)
plot(lambda, log(pl), type='l', bty='n', main = "Period", ylab = "log posterior",
    cex.main = 2.0, cex.lab = 1.5)
#plot(lambda/12, log(pl), type='l', bty='n', main = "Period (in years)", ylab = "log posterior",
#    cex.main = 2.0, cex.lab = 1.5)
#dev.off()

# par(mfrow = c(1,1), mar = c(5.1, 4.1, 4.1, 2.1))
# plot(omega, pw, type='l', bty='n', main = "Frequency", ylab = "posterior")
# plot(lambda, pl, type='l', bty='n', main = "Period", ylab = "posterior")

# plot(omega, pw, type='l')

pp = pw / sum(pw)
cc = pl / sum(pl)
bb = sample(omega, 100000, replace = TRUE, prob = pp)
ll = 2*pi/bb


plot(table(bb) / length(bb), xlim = c(0, 0.2), bty='n')
lines(omega, pp, col='red')

plot(table(ll) / length(ll), bty='n')
lines(lambda, cc, col='red')

mean(bb)
range(bb)
var(bb)

as.numeric(names(which.max(table(bb))))
table(bb)
