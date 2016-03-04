y = read.table("./dataexam2th.txt", header = TRUE)
matplot(y, type='l')

n = NROW(y) # T
p = NCOL(y) # I

qphi = function(phi_i, i)
    sum((y[-1,i] - phi_i[1]*y[-n,i])^2)

# Speed it up
sum_cross = double(p)
sum_square = double(p)
for (i in 1:p){
    sum_cross[i] = sum(y[-1,i]*y[-n,i])
    sum_square[i] = sum(y[-n,i]^2)
    }

# Part 1
tau2 = 0.1
nburn = 2000
nmcmc = 10000

nu = double(nburn + nmcmc)
phi_i = matrix(0, nburn + nmcmc, p)
phi = double(nburn + nmcmc)

nu[1] = 1

for (b in 2:(nburn + nmcmc)){

    # Update nu
    temp = 0
    for (i in 1:p)
        temp = temp + qphi(phi_i[b-1,i], i)
    nu[b] = 1/rgamma(1, n*p/2, temp/2)

    # Update phi_i
    for (i in 1:p)
        phi_i[b,i] = rnorm(1, (tau2*sum_cross[i] + nu[b]*phi[b-1]) / (tau2*sum_square[i] + nu[b]), 
            sqrt((nu[b]*tau2) / (tau2*sum_square[i] + nu[b])))

    # Update phi
    phi[b] = rnorm(1, mean(phi_i[b,]), sqrt(tau2/p))
    }

# Burn-in
nu = tail(nu, nmcmc)
phi_i = tail(phi_i, nmcmc)
phi = tail(phi, nmcmc)

### Not really necessary for Gibbs samples
# ### Trace plots
# plot(tail(nu, 1000), type='l')
# for (i in 1:p){
#     plot(tail(phi_i[,i], 1000), type='l')
#     readline()
#     }
# plot(tail(phi, 1000), type='l')
# 
# ### ACF
# acf(nu)
# for (i in 1:p){
#     acf(phi_i[,i])
#     readline()
#     }
# acf(phi)


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
pred = matrix(0, nmcmc, n)
pred[,1] = y[1,1]
for (i in 1:nmcmc){
    for (t in 2:n){
        pred[i,t] = phi_i[i,1] * pred[i,t-1] + rnorm(1, 0, sqrt(nu[i]))
        }
    }
mean.pred = apply(pred, 2, mean)

plot(y[,1], type='l')
lines(pred[4,], col = 'red')
lines(mean.pred, col = 'red')

y[1,1]*phi_i[1,1]
