library(MASS)
library(coda)

y = read.table("./dataexam2th.txt", header = TRUE)
y = head(y, 20)
matplot(y, type='l', lty=1)

n = NROW(y) # T
p = NCOL(y) # I

qstar = function(phi, i)
    y[1,i]^2 * (1 - phi^2) + sum((y[-1,i] - phi*y[-n,i])^2)

autotune = function(accept, target = 0.25, k = 2.5)
    (1+(cosh(accept-target)-1)*(k-1)/(cosh(target-
        ceiling(accept-target))-1))^sign(accept-target)


# Part 2
tau2 = 0.1
nburn = 5000
nmcmc = 10000

#mark1 = 1000
#mark2 = 3000
window = 200

nu = double(nburn + nmcmc)
phi_i = matrix(0, nburn + nmcmc, p)
phi = double(nburn + nmcmc)

accept = matrix(0, nburn+nmcmc, p)

#cand.sig = 0.001*diag(p)
cand.sig = rep(0.1, p)


# Initial values
nu[1] = 1


for (b in 2:(nburn + nmcmc)){

    # Update nu
    temp = 0
    for (i in 1:p)
        temp = temp + qstar(phi_i[b-1,i], i)

    nu[b] = 1/rgamma(1, n*p/2, temp/2)

    # Update phi_i
    for (i in 1:p){
        phi_i[b,i] = phi_i[b-1,i]
        cand = rnorm(1, phi_i[b,i], cand.sig[i])
        if (cand < 1 && cand > -1){
            cand.post = -0.5/nu[b] * qstar(cand,i) +
                0.5 * log(1-cand^2) -
                0.5/tau2 * (cand - phi[b-1])^2
            curr.post = -0.5/nu[b] * qstar(phi_i[b,i],i) +
                0.5 * log(1-phi_i[b,i]^2) -
                0.5/tau2 * (phi_i[b,i] - phi[b-1])^2
            if (log(runif(1)) < cand.post - curr.post){
                phi_i[b,i] = cand
                accept[b,i] = 1
                }
            }
        if ((b/window == floor(b / window)) && (b <= nburn))
            cand.sig[i] = cand.sig[i] * autotune(mean(accept[(b-window+1):b,i]),
                target = 0.25, k = max(window / 50, 1.5))
        }

#   # Tune it
#   if (b == mark2)
#       cand.sig = cov(phi_i[(mark1+1):mark2,])

#   if (b == nburn)
#       accept = 0

#   if ((b > mark2) && ((b / window) == floor(b / window)) && (b <= nburn))
#       cand.sig = cand.sig * autotune(sum(accept)/(b-mark2))
    
    # Update phi
    phi[b] = rnorm(1, mean(phi_i[b,]), sqrt(tau2/p))
    }

# Burn-in
nu = tail(nu, nmcmc)
phi_i = tail(phi_i, nmcmc)
phi = tail(phi, nmcmc)
accept = tail(accept, nmcmc)

apply(accept, 2, mean)

pairs(phi_i)

### Not really necessary for Gibbs samples
### Trace plots
plot(tail(nu, 1000), type='l')
for (i in 1:p){
    plot(tail(phi_i[,i], 1000), type='l')
    readline()
    }
plot(tail(phi, 1000), type='l')

### ACF
acf(nu)
for (i in 1:p){
    acf(phi_i[,i])
    readline()
    }
acf(phi)

apply(phi_i, 2, effectiveSize)


### Posterior densities
#par(mfrow = c(2,1), mar = c(4.1, 2.1, 2.1, 1.1), oma = c(0,0,4,0))
par(mfrow = c(2,1), mar = c(4.1, 2.1, 2.1, 1.1))
plot(density(nu), lwd = 3, main = expression(nu), xlab="", ylab="")

dens = NULL
for (i in 1:p)
    dens[[length(dens)+1]] = density(phi_i[,i])
dens[[length(dens)+1]] = density(phi)

plot(dens[[p+1]], lwd = 3, ylim = c(0, max(sapply(dens, function(x) max(x$y)))),
    xlim = range(sapply(dens, function(x) range(x$x))), xlab="", ylab="",
    main = expression(phi~","~phi[1]~","~phi[2]~","~phi[3]~","~phi[4]~","~phi[5]))
for (i in 1:p)
    lines(dens[[i]], col = i+1, lwd = 2)
legend("topright", box.lty = 0, legend = c(expression(phi), expression(phi[1]),
    expression(phi[2]), expression(phi[3]), expression(phi[4]), expression(phi[5])),
    col = 1:6, lwd = c(3,rep(2,p)))
par(mfrow = c(1,1), mar = c(5.1, 4.1, 4.1, 2.1))
#mtext("Posterior distributions", outer = TRUE, at = 0.05, cex = 1.5, line = 2, adj = 0)
#mtext("Based on conditional likelihood", outer = TRUE, at = 0.05, cex = 1.2, line = 0.8, adj = 0)

### Summary statistics
params = cbind(nu, phi_i, phi)
colnames(params) = c("nu", paste("phi_", 1:p, sep=""), "phi")

cbind("mean"=apply(params, 2, mean),
    "var"=apply(params, 2, var),
    t(apply(params, 2, quantile, c(0.025, 0.975))))

### Predictions
# k = 5
# pred = matrix(0, nmcmc, n)
# pred[,1] = y[1,k]
# for (i in 1:nmcmc)
#     for (t in 2:n)
#         pred[i,t] = phi_i[i,k] * pred[i,t-1] + rnorm(1, 0, sqrt(nu[i]))
# mean.pred = apply(pred, 2, mean)
# q.pred = apply(pred, 2, quantile, c(0.025, 0.975))
# 
# plot(y[,k], type='l', ylim = range(q.pred))
# lines(mean.pred, col = 'red')
# lines(q.pred[1,], col = 'red')
# lines(q.pred[2,], col = 'red')


# pred_y = matrix(0, 5000, n)
# pred_y[,1] = 0
# for (i in 1:nrow(pred_y))
#     for (t in 2:n)
#         pred_y[i,t] = 0.5 * pred_y[i,t-1] + rnorm(1, 0, sqrt(1.06))
# 
# mpred = apply(pred_y, 2, mean)
# qpred = apply(pred_y, 2, quantile, c(0.025, 0.975))
# plot(mpred, type='l', ylim = range(qpred))
# lines(qpred[1,], col = 'red')
# lines(qpred[2,], col = 'red')

