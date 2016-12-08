### Settings
set.seed(35)
Nmax = 75*1
k = 4
niter = ceiling(Nmax / k)
nmcmc = 10000
nsims = 5
delta.cr = 0.20
delta.tox = 0.05
l.cr = 0.02
u.tox = 0.80

### Initialize
stopped = matrix(0, niter, nsims)
pi.e.cr = matrix(0, niter, nsims)
pi.e.tox = matrix(0, niter, nsims)


 probs = c(0.2, 0.4) # scenario 1
#probs = c(0.4, 0.2) # scenario 2
pvec = double(4)

pvec[1] = probs[1] * (1 - probs[2])        # CR, ~TOX,     A1
pvec[2] = probs[1] * probs[2]              # CR, TOX,      A2
pvec[3] = (1 - probs[1]) * (1 - probs[2])  # ~CR, ~TOX,    A3
pvec[4] = (1 - probs[1]) *  probs[2]       # ~CR, TOX,     A4

for (l in 1:nsims){
    theta.s = c(2.037, 6.111, 30.555, 2.037)

    theta.e = theta.s
    theta.e = theta.e * 4 / sum(theta.e)

    eta.s.cr = rbeta(nmcmc, theta.s[1] + theta.s[2], theta.s[3] + theta.s[4])
    eta.s.tox = rbeta(nmcmc, theta.s[2] + theta.s[4], theta.s[1] + theta.s[3])

    for (i in 1:niter){
        # Sample a new cohort
        y = sample(1:k, k, replace = TRUE, prob = pvec)

        # Update theta.e
        for (j in 1:4)
            theta.e[j] = theta.e[j] + sum(y == j)

        # Sample from eta.e
        eta.e.cr = rbeta(nmcmc, theta.e[1] + theta.e[2], theta.e[3] + theta.e[4])
        eta.e.tox = rbeta(nmcmc, theta.e[2] + theta.e[4], theta.e[1] + theta.e[3])

        # Estimate pi.e.cr and pi.e.tox
        pi.e.cr[i,l] = mean(eta.e.cr > (eta.s.cr + delta.cr))
        pi.e.tox[i,l] = mean(eta.e.tox > (eta.s.tox + delta.tox))

        # Update stop vector
        if (pi.e.cr[i,l] < l.cr)
            stopped[i,l] = 1
        if (pi.e.tox[i,l] > u.tox)
            stopped[i,l] = 2
        if ((pi.e.cr[i,l] < l.cr) && (pi.e.tox[i,l] > u.tox))
            stopped[i,l] = 3
        if ((i > 1) && (stopped[i-1,l] != 0))
            stopped[i,l] = stopped[i-1,l]
        }
    }

pdf("./figs/scen1.pdf", width = 8, height = 9)
par(mfrow = c(3,2), mar = c(4.1, 2.1, 2.1, 1.1))
for (l in 1:nsims){
    plot(pi.e.cr[,l], ylim = c(0, 1), col = 'green', type='b', axes = FALSE,
        ylab = "", xlab = "n")
    axis(2)
    axis(1, at = (1:niter), labels = (1:niter)*k)
    points(pi.e.tox[,l], col = 'red', type='b')
    abline(h = u.tox, lty = 2)
    abline(h = l.cr, lty = 2)
    if (max(stopped[,l]) != 0){
        if (max(stopped[,l]) == 1)
            col = 'green'
        if (max(stopped[,l]) == 2)
            col = 'red'
        if (max(stopped[,l]) == 3)
            col = 'blue'
        abline(v = which.max(stopped[,l]), lty = 2, col =col)
        }
    }
plot(0, type='n', axes = FALSE, xlim = c(0,1), ylim = c(0,1), xlab = "", ylab ="")
text(0.1, 0.8, expression(pi[n](TOX)), cex = 3, adj = 0, col = 'red')
text(0.1, 0.5, expression(pi[n](CR)), cex = 3, adj = 0, col = 'green')
par(mfrow = c(1,1), mar = c(5.1, 4.1, 4.1, 2.1))
dev.off()

