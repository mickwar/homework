library(splines)
source("~/files/R/mcmc/bayes_functions.R")

dat = read.table("~/files/data/eu_poll.txt", header = TRUE)
n = NROW(dat)

#key_events = c("5-Jan-2016", "2-Feb-2016", "20-Feb-2016", "21-Feb-2016",
#    "20-Mar-2016", "15-Apr-2016", "18-Apr-2016", "22-Apr-2016")
#key_time = as.numeric(as.Date(key_events, "%d-%b-%Y")) - 16678
#key_time = key_time[c(1, 2, 5, 6, 8)]

logit = function(p) log(p/(1-p))
ilogit = function(x) 1/(1+exp(-x))

table(dat[,5:6])

### Re-order (so earliest poll is first) 
dat = dat[n:1,]
dat[,5] = as.character(dat[,5])
dat[,6] = as.character(dat[,6])

### Combine the companies which did few polls
dat[dat[,5] == "BMG" | dat[,5] == "Opinium" | dat[,5] == "Panelbase", 5] = "Other"
#dat[dat[,5] == "BMG" | dat[,5] == "Opinium" | dat[,5] == "Panelbase" |
#    dat[,5] == "IpsosMori" | dat[,5] == "TNS", 5] = "Other"

### Orginaize data
# Observations
y = logit(dat[,2] / (dat[,2] + dat[,3]))

# Time vector
time = as.numeric(as.Date(dat[,1], "%d-%b-%Y")) - 16678

# Design matrices
# Treat Other and Online as baselines
company = matrix(model.matrix(lm(y ~ 0 + factor(dat[,5]))), nrow=n)
colnames(company) = c(sort(unique(dat[,5])))
company = company[,-which(colnames(company) == "Other")]
type = matrix(model.matrix(lm(y ~ 0 + factor(dat[,6]))), nrow=n)
colnames(type) = c(sort(unique(dat[,6])))
type = type[,-which(colnames(type) == "online")]

quantile(time, seq(0,1,length=5))

(max(time) - min(time)) / 4
time

pdf("figs/data.pdf", width = 6, height = 6)
plot(time, y, type='l', xlab = "", axes = FALSE, ylab = expression("Logit ("~y[t]~")"),
    cex.lab = 1.2, main = "Voting preference data")
axis(1, at = time[c(1, 17, 34, 61, 100)], labels = dat[c(1, 17, 34, 61, 100),1])
axis(2)
points(time, y, col = 3+type, pch = 20, cex = 1.5)
dev.off()

#####
nburn = 200
nmcmc = 10000

design.matrices = list(
    matrix(rep(1, n), n, 1, dimnames = list(NULL, "Intercept")),
    cbind("Intercept"=1, time),
    cbind("Intercept"=1, company),
    cbind("Intercept"=1, type),
    cbind("Intercept"=1, time, company),
    cbind("Intercept"=1, time, type),
    cbind("Intercept"=1, company, type),
    cbind("Intercept"=1, time, company, type))

model.DIC = double(length(design.matrices))
model.gof = matrix(0, nmcmc, length(design.matrices))
model.ppl = matrix(0, 2, length(design.matrices))


for (model.n in 1:length(design.matrices)){
    X = design.matrices[[model.n]]

    ### Gibbs sampler
    p = NCOL(X)
    prior.m = rep(0, p)
    #prior.m = solve(t(X) %*% X) %*% t(X) %*% y
    prior.g = n
    prior.a = 0
    prior.b = 0

    A = t(X) %*% X
    chol.A = t(chol(solve(A)))
    post.mean = 1/(prior.g + 1) * (prior.g * solve(t(X) %*% X) %*% t(X) %*% y + prior.m)

    param.beta = matrix(0, nburn + nmcmc, p)
    param.sig2 = double(nburn + nmcmc)
    param.sig2[1] = 1

    for (i in 2:(nburn + nmcmc)){
        param.beta[i,] = post.mean +
            sqrt(prior.g * param.sig2[i-1] / (prior.g + 1)) * chol.A %*% rnorm(p, 0, 1)
        param.sig2[i] = 1/rgamma(1,
            prior.a + n/2 + p/2, 
            prior.b + 0.5*sum((y - X %*% param.beta[i,])^2) +
                0.5/prior.g * t(param.beta[i,] - prior.m) %*% A %*% (param.beta[i,] - prior.m))
        }

    # Truncate
    param.beta = tail(param.beta, nmcmc)
    param.sig2 = tail(param.sig2, nmcmc)

    ### Posterior predictions
    pred.y = matrix(0, n, nmcmc)
    for (i in 1:nmcmc)
        pred.y[,i] = rnorm(n, X %*% param.beta[i,], 1*sqrt(param.sig2[i]))

    ### DIC
    dtheta = matrix(0, nmcmc)
    for (i in 1:nmcmc)
        dtheta[i] = -2*sum(dnorm(y, X %*% param.beta[i,], sqrt(param.sig2[i]), log = TRUE))
    model.DIC[model.n] = mean(dtheta) + var(dtheta)/2

    ### Bayes GOF
    model.gof[,model.n] = bayes.gof(y, cbind(param.beta, param.sig2),
        function(y, param) pnorm(y, X %*% param[1:NCOL(X)], sqrt(param[NCOL(X)+1])),
        every = nmcmc + 1)

    ### PPL
    model.ppl[1, model.n] = sum((y - apply(pred.y, 1, mean))^2)
    model.ppl[2, model.n] = sum(apply(pred.y, 1, var))

    cat("Complete:", model.n, "/", length(design.matrices), "\n")
    }

pdf("figs/gof.pdf", width = 6, height = 6)
plot(0, type='n', xlim = c(-0.1, 1.1), ylim = c(0, 4), axes = FALSE, ylab = "Density",
    main = "Bayesian goodness-of-fit p-values", xlab = "p-value")
for (j in 1:8)
    lines(density(model.gof[,j]), col = j, lwd = 2)
legend("topright", legend = paste0("M", 1:8), col = 1:8, box.lty = 0, lwd = 3, cex = 1.2)
axis(1)
axis(2)
dev.off()

cbind(model.DIC, t(model.ppl), apply(model.ppl, 2, sum))


### Plots
mm = apply(pred.y, 1, mean)
qq = apply(pred.y, 1, quantile, c(0.025, 0.975))
pdf("figs/pred_obsfit.pdf", width = 6, height = 6)
plot(y, mm, pch = 20, ylim = range(qq), col = 'darkblue', xlab = "Observed",
    ylab = "Predicted", main = "M8 model predictions", cex.lab = 1.3)
abline(0, 1)
segments(x0 = y, y0 = qq[1,], x1 = y, y1 = qq[2,], col = 'dodgerblue')
points(y, mm, pch = 20, col = 'darkblue')
dev.off()

hpd.beta = apply(param.beta, 2, hpd.uni)
hpd.sig2 = hpd.uni(param.sig2)

colnames(X)[9] = "Phone"
hpd.plot(density(param.sig2), hpd.sig2, main = expression("Model Variance"),
    xlab = expression(sigma^2))
pdf("figs/post.pdf", width = 9, height = 9)
par(mfrow = c(ceiling(sqrt(p)), ceiling(sqrt(p))), mar = c(4.1, 2.1, 2.1, 1.1))
for (i in 1:p){
    hpd.plot(density(param.beta[,i]), hpd.beta[,i], main = colnames(X)[i],
        xlab = bquote(beta[.(i-1)]))
    abline(v = 0, col = 'red')
    }
par(mfrow = c(1,1), mar = c(5.1, 4.1, 4.1, 2.1))
dev.off()


# plot(time, y, type='l', ylim = range(qq))
# points(time, y, pch = 20)
# lines(time, mm, col = 'red')
# lines(time, qq[1,], col = 'blue')
# lines(time, qq[2,], col = 'blue')
# abline(v = key_time, col = 'green')

pdf("figs/pred_time.pdf", width = 10, height = 6)
plot(time, ilogit(y), type='l', ylim = range(ilogit(qq)), axes = FALSE,
    xlab = "", ylab = "Probability of Leaving", lwd = 3, main = "M8 predictions over time")
points(time, ilogit(y), pch = 20)
lines(time, ilogit(mm), col = 'dodgerblue', lwd = 2)
lines(time, ilogit(qq[1,]), col = col.mult('dodgerblue', 'gray50'), lwd = 1)
lines(time, ilogit(qq[2,]), col = col.mult('dodgerblue', 'gray50'), lwd = 1)
axis(1, at = time[c(1, 17, 34, 61, 100)], labels = dat[c(1, 17, 34, 61, 100),1])
axis(2)
dev.off()
#lines(time, ilogit(y), pch = 20, lwd = 3)
#abline(h = 0.5)
#abline(v = key_time, col = 'green')
    
