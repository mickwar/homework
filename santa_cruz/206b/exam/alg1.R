library(xtable)

y = read.table("./dataexam2th.txt", header = TRUE)
#y = head(y, 20)
#pdf("./figs/data.pdf", height = 8, width = 8)
#par(mar = c(4.1, 2.1, 2.1, 1.1))
matplot(y, type='l', lty=1, col = 2:6)
#dev.off()

n = NROW(y) # T
p = NCOL(y) # I

# Speed it up
sum_cross = double(p)
sum_square = double(p)
for (i in 1:p){
    sum_cross[i] = sum(y[-1,i]*y[-n,i])
    sum_square[i] = sum(y[-n,i]^2)
    }

# Part 1
tau2 = 0.1
nburn = 5000
nmcmc = 20000

nu = double(nburn + nmcmc)
phi_i = matrix(0, nburn + nmcmc, p)
phi = double(nburn + nmcmc)

nu[1] = 1

for (b in 2:(nburn + nmcmc)){

    # Update nu
    temp = 0
    for (i in 1:p)
        temp = temp + sum((y[-1,i] - phi_i[b-1,i]*y[-n,i])^2)
    nu[b] = 1/rgamma(1, n*p/2, temp/2)

    # Update phi_i
    for (i in 1:p)
        phi_i[b,i] = rnorm(1, (tau2*sum_cross[i] + nu[b]*phi[b-1]) /
            (tau2*sum_square[i] + nu[b]), 
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
pdf("./figs/posts_1.pdf", width = 8, height = 8)
par(mfrow = c(2,1), mar = c(3.1, 2.1, 2.1, 1.1))
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
    col = 1:6, lwd = c(3,rep(2,p)), cex = 1.5)
par(mfrow = c(1,1), mar = c(5.1, 4.1, 4.1, 2.1))
dev.off()

### Summary statistics
params = cbind(nu, phi_i, phi)
colnames(params) = c("nu", paste("phi_", 1:p, sep=""), "phi")

xtable(cbind("mean"=apply(params, 2, mean),
    "sd"=apply(params, 2, sd),
    t(apply(params, 2, quantile, c(0.025, 0.975)))))


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

