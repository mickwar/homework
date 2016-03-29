orbits = read.table("~/files/data/205b/orbits.txt", header = TRUE)
n = NROW(orbits)

f1 = function(p, log = TRUE){
    alpha = p[1]
    beta = p[2]
    sig2 = p[3]
    out = -n/2 * log(2*pi*sig2) - 1/(2*sig2)*sum((orbits$Y - alpha -
        beta*cos(omega.hat * orbits$X + delta.hat))^2)
    if (log)
        return (out)
    return (exp(out))
    }
f2 = function(p, log = TRUE){
    omega = p[1]
    delta = p[2]
    out = -n/2 * log(2*pi*sig2.hat) - 1/(2*sig2.hat)*sum((orbits$Y - alpha.hat -
        beta.hat*cos(omega * orbits$X + delta))^2)
    if (log)
        return (out)
    return (exp(out))
    }
f3 = function(p, log = TRUE){
    alpha = p[1]
    beta = p[2]
    omega = p[3]
    delta = p[4]
    sig2 = p[5]
    out = -n/2 * log(2*pi*sig2) - 1/(2*sig2)*sum((orbits$Y - alpha -
        beta*cos(omega * orbits$X + delta))^2)
    if (log)
        return (out)
    return (exp(out))
    }

### (a)
### Profile likelihood
# Initial values for the "fixed" parameters
alpha.hat = mean(c(range(orbits$Y)))    # Estimate for the intercept 
beta.hat = diff(range(orbits$Y))/2      # Estimate for amplitude
sig2.hat = 0.1                          # Guess for variance

# The starting values for the optimizer
xy.1 = expand.grid(
    mean(c(range(orbits$Y))) + seq(-0.2, 0.2, length = 10),
    diff(range(orbits$Y))/2 + seq(-0.2, 0.2, length = 10),
    seq(0.005, 0.1, length = 10))
xy.2 = expand.grid(seq(0, 3, length = 40), seq(-pi, pi, length = 40))

# Iterate through the optimization 5 times (though only 2 may be really necessary)
for (j in 1:5){
    temp.2 = Inf
    for (i in 1:nrow(xy.2)){
        # Optimize treating alpha, beta, sig^2 as fixed (i.e. using alpha.hat,
        # beta.hat, and sig2.hat as the fixed values)
        temp = optim(as.double(xy.2[i,]), function(x) -f2(x, log = TRUE),
            method = "L-BFGS-B", lower = c(0, -pi), upper = c(3, pi))
        # If the i'th starting points produced a better mode, update the parameters
        if (temp$value < temp.2){
            temp.2 = temp$value
            omega.hat = temp$par[1]
            delta.hat = temp$par[2]
            }
        }
    temp.2 = Inf
    for (i in 1:nrow(xy.1)){
        # Optimize treating omega, delta as fixed (using omega.hat and delta.hat)
        temp = optim(as.double(xy.1[i,]), function(x) -f1(x, log = TRUE),
            method = "L-BFGS-B", lower = c(min(orbits$Y), min(orbits$Y), 0.001),
            upper = c(max(orbits$Y), max(orbits$Y), var(orbits$Y)))
        # If the i'th starting points produced a better mode, update the parameters
        if (temp$value < temp.2){
            temp.2 = temp$value
            alpha.hat = temp$par[1]
            beta.hat = temp$par[2]
            sig2.hat = temp$par[3]
            }
        }
    }


f3(c(alpha.hat, beta.hat, omega.hat, delta.hat, sig2.hat))
c(alpha.hat, beta.hat, omega.hat, delta.hat, sig2.hat)

# Plots
xx = seq(-1, 13, length = 1000)
plot(orbits, pch = 20, ylim = range(orbits$Y) + c(-2,2)*sqrt(sig2.hat))
lines(xx, alpha.hat+beta.hat*cos(omega.hat*xx+delta.hat), type='l')
lines(xx, alpha.hat+beta.hat*cos(omega.hat*xx+delta.hat)+2*sqrt(sig2.hat), type='l', col = 'red')
lines(xx, alpha.hat+beta.hat*cos(omega.hat*xx+delta.hat)-2*sqrt(sig2.hat), type='l', col = 'red')

# Approximate predictive coverage
mean(orbits$Y <= alpha.hat+beta.hat*cos(omega.hat*orbits$X+delta.hat)+2*sqrt(sig2.hat) &
    orbits$Y >= alpha.hat+beta.hat*cos(omega.hat*orbits$X+delta.hat)-2*sqrt(sig2.hat))

### (b)
### Wald-like
I.omega = beta.hat / sig2.hat * sum(
    ((orbits$Y - alpha.hat - beta.hat*cos(omega.hat*orbits$X + delta.hat)) *
    cos(omega.hat*orbits$X + delta.hat)*orbits$X^2) +
    (beta.hat * orbits$X^2 * sin(omega.hat * orbits$X + delta.hat)^2))

I.delta = beta.hat /sig2.hat * sum(
    ((orbits$Y - alpha.hat - beta.hat*cos(omega.hat*orbits$X + delta.hat)) *
    cos(omega.hat * orbits$X + delta.hat)) +
    beta.hat*sin(omega.hat*orbits$X + delta.hat)^2)

# Approx conf int for omega and delta
omega.hat + 1/sqrt(I.omega) * qnorm(0.975) * c(-1, 1)
delta.hat + 1/sqrt(I.delta) * qnorm(0.975) * c(-1, 1)

### (c)
### Bootstrap confidence intervals
B = 5000
boot.par = matrix(0, B, 2) # 2 columns for omega and delta
for (b in 1:B){
    # Generate a bootstrap sample (replace orbits$Y because of how I coded f3)
    orbits$Y = rnorm(n, alpha.hat + beta.hat * cos(omega.hat * orbits$X + delta.hat),
        sqrt(sig2.hat))
    # Compute mles for the bootstrap sample, with the actual mles as starting points
    # Not doing the iterative profile likelihood as before since it takes too long
    # and we are already near the correct values so there shouldn't be any issues
    temp = optim(c(alpha.hat, beta.hat, omega.hat, delta.hat, sig2.hat),
        function(x) -f3(x, log = TRUE), method = "L-BFGS-B",
        lower = c(min(orbits$Y), min(orbits$Y), 0, -pi, 0.001),
        upper = c(max(orbits$Y), max(orbits$Y), 3, pi, var(orbits$Y)))
    boot.par[b,] = temp$par[c(3,4)]
    }

quantile(boot.par[,1], c(0.025, 0.975))
quantile(boot.par[,2], c(0.025, 0.975))

plot(density(boot.par[,1]), main = expression(omega), xlab = "")
plot(density(boot.par[,2]), main = expression(delta), xlab = "")

# Comparison
rbind(omega.hat + 1/sqrt(I.omega) * qnorm(0.975) * c(-1, 1),
    quantile(boot.par[,1], c(0.025, 0.975)))

rbind(delta.hat + 1/sqrt(I.delta) * qnorm(0.975) * c(-1, 1),
    quantile(boot.par[,2], c(0.025, 0.975)))


### Problem 2
### (b)
life = read.table("~/files/data/205b/lifetime.txt", header = TRUE)
n = nrow(life)
tobs = mean(life$x) - mean(life$y)
lhat = 0.5*(mean(life$x) + mean(life$y))
I.lambda = -(2*n/(lhat^2) - 2/(lhat^3) * (sum(life$x) + sum(life$y)))
pnorm(tobs / (1 / sqrt(I.lambda)),0,1,lower.tail = FALSE)

### (c)
### Permutation test
t.perm = double(10000)
w = as.double(unlist(life))
for (i in 1:length(t.perm)){
    s = sample(w)
    t.perm[i] = mean(head(s, 20)) - mean(tail(s, 20))
    }
mean(t.perm >= tobs)


### Problem 3
n = c(10, 100, 1000)
alpha = 0.1

pdf("./power.pdf", width = 8, height = 12)
par(mfrow = c(3,1), mar=c(4.1, 4.1, 2.1, 1.1))
for (i in 1:3){
    qq = seq(0, 2, length = 100)    # (theta_0 - theta)/S
    plot(qq, pt(-sqrt(n[i])*qq + qt(alpha, n[i]-1, lower.tail = FALSE), n[i]-1, lower.tail = FALSE), type='l',
        xlab = expression((theta - theta[0])/S),
        ylab = expression(beta(theta)), lwd = 1, main = paste0("n=",n[i]))
    lines(qq, pnorm(-sqrt(n[i])*qq + qnorm(alpha, lower.tail = FALSE), lower.tail = FALSE), type='l', col = 'red', lwd = 1)
    legend("topleft", box.lty = 0, legend=c(expression(beta[1](theta)), expression(beta[2](theta))), col = c(1,2), lwd = 2)
    }
dev.off()

