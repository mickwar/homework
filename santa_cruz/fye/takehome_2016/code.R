source("~/files/R/mcmc/bayes_functions.R")
dat = read.table("~/files/data/fye_2016_data.txt", header = TRUE)
#dat[,1] = log(dat[,1]+1)
n = NROW(dat)

summary(dat)
#pdf("figs/pairs.pdf", width = 6, height = 6)
pairs(dat, pch = 20,
#   upper.panel = function(x, y, ...){ points(x, y, ...) },
    lower.panel = function(x, y, ...){
        text(x = sum(range(x))/2, y = sum(range(y))/2,
            label = as.character(round(cor(x,y), 3)), cex=2.5)
        }
    )
#dev.off()

library(rgl)
plot3d(dat[,3], dat[,4], dat[,1], col = dat[,2]+1, box = FALSE, axes = FALSE,
    xlab = "", ylab = "", zlab = "")
for (i in 1:11){
    lines3d(rbind(c(min(dat[,3]), min(dat[,4])+(i-1)*diff(range(dat[,4]))/10, 0),
        c(max(dat[,3]), min(dat[,4])+(i-1)*diff(range(dat[,4]))/10, 0)),
        col = 'gray75')
    lines3d(rbind(c(min(dat[,3])+(i-1)*diff(range(dat[,3]))/10, min(dat[,4]), 0),
        c(min(dat[,3])+(i-1)*diff(range(dat[,3]))/10, max(dat[,4]), 0)),
        col = 'gray75')
    }
z = double(2*n)
z[seq(2, 2*n, by = 2)] = dat[,1]
points3d(dat[,3], dat[,4], dat[,1], col = dat[,2]+1, size = 5)
segments3d(rep(dat[,3], each = 2), rep(dat[,4], each = 2), z,
    col = "gray30")
axis3d('x-'); mtext3d(colnames(dat)[3], 'x-', 2.5)
axis3d('y-'); mtext3d(colnames(dat)[4], 'y-', 2.5)
axis3d('z+'); mtext3d(colnames(dat)[1], 'z+', 2.5)
#rgl.postscript("figs/data.pdf", "pdf")

###
y = dat[,1]
rad = dat[,2]
temp = dat[,3]
wind = dat[,4]
n = length(y)

library(splines)
design.matrices = list(
    matrix(rep(1, n), n, 1, dimnames = list(NULL, "Intercept")),
    cbind("Intercept"=1, rad, temp, wind),
    cbind("Intercept"=1, temp, wind, temp*wind),
    cbind("Intercept"=1, temp, wind, temp^2, wind^2, temp*wind),
    cbind("Intercept"=1, rad, temp, wind, temp*wind),
    cbind("Intercept"=1, rad, ns(temp, knots = c(75)), ns(wind, knots = c(10)),
        temp*wind),
    cbind("Intercept"=1, rad, ns(temp, knots = c(70, 80)), ns(wind, knots = c(5, 12)),
        temp*wind),
    cbind("Intercept"=1, rad, ns(temp, knots = c(70, 80)), ns(wind, knots = c(5, 12)),
        ns(temp, knots = c(70, 80))*ns(wind, knots = c(5, 12))))#,
#   cbind("Intercept"=1, rad, ns(temp, knots = c(70, 80)), ns(wind, knots = c(8, 12)),
#       rad*temp, rad*wind, temp*wind, rad*temp,wind))
qnum = length(design.matrices)

nburn = 200
nmcmc = 10000

model.DIC1 = double(qnum)
model.DIC2 = double(qnum)
model.gof = matrix(0, nmcmc, qnum)
model.ppl = matrix(0, 2, qnum)

#model.DIC = c(model.DIC, 0)
#model.gof = cbind(model.gof, 0)
#model.ppl = cbind(model.ppl, 0)

for (model.n in 1:qnum){
    X = design.matrices[[model.n]]


    p = NCOL(X)
    prior.m = rep(0, p)
    prior.g = n
    prior.a = 0
    prior.b = 0

    A = t(X) %*% X
    max(eigen(A)$values) / min(eigen(A)$values)
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
    model.DIC1[model.n] = mean(dtheta) + var(dtheta)/2
    model.DIC2[model.n] = 2*mean(dtheta) + 2*sum(dnorm(y,
        X %*% apply(param.beta, 2, mean), sqrt(mean(param.sig2))))

    ### Bayes GOF
    model.gof[,model.n] = bayes.gof(y, cbind(param.beta, param.sig2),
        function(y, param) pnorm(y, X %*% param[1:NCOL(X)], sqrt(param[NCOL(X)+1])),
        every = nmcmc + 1)

    ### PPL
    model.ppl[1, model.n] = sum((y - apply(pred.y, 1, mean))^2)
    model.ppl[2, model.n] = sum(apply(pred.y, 1, var))

    cat("Complete:", model.n, "/", length(design.matrices), "\n")
    }



plot(0, type='n', xlim = c(-0.1, 1.1), ylim = c(0, 20), axes = FALSE, ylab = "Density",
    main = "Bayesian goodness-of-fit p-values", xlab = "p-value")
for (j in 1:qnum)
    lines(density(model.gof[,j]), col = j, lwd = 2)
legend("topright", legend = paste0("M", 1:qnum), col = 1:qnum, box.lty = 0, lwd = 3, cex = 1.2)
axis(1)
axis(2)
apply(model.gof, 2, function(x) mean(x > 0.05))

cbind(model.DIC1, model.DIC2, t(model.ppl), apply(model.ppl, 2, sum))

hpd.beta = apply(param.beta, 2, hpd.uni)
hpd.sig2 = hpd.uni(param.sig2)

hpd.plot(density(param.sig2), hpd.sig2, main = expression("Model Variance"),
    xlab = expression(sigma^2))
#pdf("figs/post.pdf", width = 9, height = 9)
par(mfrow = c(ceiling(sqrt(p)), ceiling(sqrt(p))), mar = c(4.1, 2.1, 2.1, 1.1))
for (i in 1:p){
    hpd.plot(density(param.beta[,i]), hpd.beta[,i], main = colnames(X)[i],
        xlab = bquote(beta[.(i-1)]))
    abline(v = 0, col = 'red')
    }
par(mfrow = c(1,1), mar = c(5.1, 4.1, 4.1, 2.1))
#dev.off()


mm = apply(pred.y, 1, mean)
qq = apply(pred.y, 1, quantile, c(0.025, 0.975))
#pdf("figs/pred_obsfit.pdf", width = 6, height = 6)
plot(y, mm, pch = 20, ylim = range(qq), col = 'darkblue', xlab = "Observed",
    ylab = "Predicted", main = "M8 model predictions", cex.lab = 1.3)
abline(0, 1)
segments(x0 = y, y0 = qq[1,], x1 = y, y1 = qq[2,], col = 'dodgerblue')
points(y, mm, pch = 20, col = 'darkblue')
#dev.off()

plot(temp, y, type='n')
points(temp, mm, pch = 20, col = 'darkblue')
segments(x0 = temp, y0 = qq[1,], x1 = temp, y1 = qq[2,], col = 'dodgerblue')
points(temp, y, pch = 20, cex = 2)

plot(wind, y, type='n')
points(wind, mm, pch = 20, col = 'darkblue')
segments(x0 = wind, y0 = qq[1,], x1 = wind, y1 = qq[2,], col = 'dodgerblue')
points(wind, y, pch = 20, cex = 2)

plot3d(temp, wind, y, size = 5)
points3d(temp, wind, mm, col = 'darkblue')

plot3d(temp, wind, mm, col = 'darkblue', size = 5)
