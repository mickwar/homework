source("~/files/R/mcmc/bayes_functions.R")
dat = read.table("~/files/data/fye_2016_data.txt", header = TRUE)

### Add 1 (to handle any zeroes if any) and take log transform
dat[,1] = log(dat[,1]+1)
names(dat)[1]="log(ozone+1)"

###
y = dat[,1]
rad = dat[,2]
temp = dat[,3]
wind = dat[,4]
n = length(y)

summary(dat)

#pdf("figs/pairs.pdf", width = 6, height = 6)
pairs(dat, pch = 20,
#   upper.panel = function(x, y, ...){ points(x, y, ...) },
    lower.panel = function(x, y, ...){
        cc = round(cor(x, y), 3)
        text(x = sum(range(x))/2, y = sum(range(y))/2,
            label = as.character(cc), cex=1+2.0*abs(cc))
        },
    diag.panel = function(x, ...){
        usr = par("usr"); on.exit(par(usr))
        par(usr = c(usr[1:2], 0, 1.5) )
        h = hist(x, plot = FALSE)
        breaks = h$breaks; nB <- length(breaks)
        y = h$counts; y = y/max(y)
        rect(breaks[-nB], 0, breaks[-1], y, col = "gray50", border= "white", ...)
        box()
        }
    )
#dev.off()

library(rgl)
plot3d(temp, wind, y, col = dat[,2]+1, box = FALSE, axes = FALSE,
    xlab = "", ylab = "", zlab = "", zlim = c(0, 6))
for (i in 1:11){
    lines3d(rbind(c(min(temp), min(wind)+(i-1)*diff(range(wind))/10, 0),
        c(max(temp), min(wind)+(i-1)*diff(range(wind))/10, 0)),
        col = 'gray75')
    lines3d(rbind(c(min(temp)+(i-1)*diff(range(temp))/10, min(wind), 0),
        c(min(temp)+(i-1)*diff(range(temp))/10, max(wind), 0)),
        col = 'gray75')
    }
points3d(temp, wind, y, col = dat[,2]+1, size = 5)
segments3d(rep(temp, each = 2), rep(wind, each = 2), c(rbind(0, y)),
    col = "gray30")
axis3d('x-'); mtext3d(colnames(dat)[3], 'x-', 2.5)
axis3d('y-'); mtext3d(colnames(dat)[4], 'y-', 2.5)
axis3d('z+'); mtext3d("log(ozone+1)", 'z+', 1.5, at = 3.5)
#rgl.postscript("figs/log_data.pdf", "pdf")

plot3d(temp, wind, exp(y)-1, col = dat[,2]+1, box = FALSE, axes = FALSE,
    xlab = "", ylab = "", zlab = "")
for (i in 1:11){
    lines3d(rbind(c(min(temp), min(wind)+(i-1)*diff(range(wind))/10, 0),
        c(max(temp), min(wind)+(i-1)*diff(range(wind))/10, 0)),
        col = 'gray75')
    lines3d(rbind(c(min(temp)+(i-1)*diff(range(temp))/10, min(wind), 0),
        c(min(temp)+(i-1)*diff(range(temp))/10, max(wind), 0)),
        col = 'gray75')
    }
points3d(temp, wind, exp(y)-1, col = dat[,2]+1, size = 5)
segments3d(rep(temp, each = 2), rep(wind, each = 2), c(rbind(0, exp(y)-1)),
    col = "gray30")
axis3d('x-'); mtext3d(colnames(dat)[3], 'x-', 2.5)
axis3d('y-'); mtext3d(colnames(dat)[4], 'y-', 2.5)
axis3d('z+'); mtext3d(colnames(dat)[1], 'z+', 1.5)
#rgl.postscript("figs/data.pdf", "pdf")


design.matrices = list(
    matrix(rep(1, n), n, 1, dimnames = list(NULL, "Intercept")),
    cbind("Intercept"=1, rad),
    cbind("Intercept"=1, temp),
    cbind("Intercept"=1, wind),
    cbind("Intercept"=1, rad, temp),
    cbind("Intercept"=1, rad, temp, rad*temp),
    cbind("Intercept"=1, rad, wind),
    cbind("Intercept"=1, rad, wind, rad*wind),
    cbind("Intercept"=1, temp, wind),
    cbind("Intercept"=1, temp, wind, temp*wind),
    cbind("Intercept"=1, rad, temp, wind),
    cbind("Intercept"=1, rad, temp, wind, rad*temp),
    cbind("Intercept"=1, rad, temp, wind, rad*wind),
    cbind("Intercept"=1, rad, temp, wind, temp*wind),
    cbind("Intercept"=1, rad, temp, wind, rad*temp, rad*wind),
    cbind("Intercept"=1, rad, temp, wind, rad*temp, temp*wind),
    cbind("Intercept"=1, rad, temp, wind, rad*wind, temp*wind),
    cbind("Intercept"=1, rad, temp, wind, rad*temp, rad*wind, temp*wind))
qnum = length(design.matrices)

picked.model = c(10, 14)
picked.col = c('firebrick', 'dodgerblue')




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

    cat("Complete:", model.n, "/", qnum, "\r")
    }
cat("\n")

### Indexes for M1 and M2 (with colors)
picked = c(10, 14)
cols.A = c('firebrick', 'dodgerblue')
cols.B = c('darkred', 'darkblue')

### Matrix for DIC and PPL
cbind(model.DIC1, t(model.ppl), apply(model.ppl, 2, sum))[picked,]



### Plot of GOF p-values
pdf("figs/gof.pdf", width = 6, height = 6)
plot(0, type='n', xlim = c(-0.1, 1.1), ylim = c(0, 10), axes = FALSE, ylab = "Density",
    main = "Bayesian goodness-of-fit p-values", xlab = "p-value")
for (j in 1:2)
    lines(density(model.gof[,picked[j]]), col = cols.A[j], lwd = 3)
legend("topright", legend = c("M1", "M2"), col = cols.A,
    box.lty = 0, lwd = 3, cex = 1.2)
axis(1)
axis(2)
apply(model.gof[,picked], 2, function(x) mean(x > 0.05))
dev.off()

