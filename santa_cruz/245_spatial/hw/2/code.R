library(geoR)

### Karhunen-Loeve representation
get.lambda.psi = function(J, L){
    vl = double(J)
    vl[1] = uniroot(function(v) 1-v*tan(v*L), interval = c(0, 0.5*pi/L-0.0001))$root
    for (i in 2:J)
        vl[i] = uniroot(function(v) 1-v*tan(v*L),
            interval = c((i-1.5)*pi/L+0.0001, (i-0.5)*pi/L-0.0001))$root


    wk = double(J)
    wk[1] = uniroot(function(w) w+tan(w*L), interval = c(0, 0.5*pi/L-0.0001))$root
    for (i in 1:J)
        wk[i] = uniroot(function(w) w+tan(w*L),
            interval = c((i-1.5)*pi/L+0.0001, (i-0.5)*pi/L-0.0001))$root

    lambda = double(J)
    for (i in 1:J){
        if (i %% 2 == 1){   # odd
            lambda[i] = 2/(1 + vl[i]^2)
        } else {            # even
            lambda[i] = 2/(1 + wk[i]^2)
            }
        }
    psi = function(s, i){
        out = matrix(0, length(s), length(i))
        for (j in 1:length(s)){
            out[j,] = ifelse(i %% 2 == 1,
                cos(vl[i]*s[j])/sqrt(L + sin(2*vl[i]*L)/(2*vl[i])),
                sin(wk[i]*s[j])/sqrt(L - sin(2*wk[i]*L)/(2*wk[i])))
            }
        out = ifelse(is.na(out), 0, out)
        return (out)
        }

    return (list("lambda" = lambda, "psi" = psi))
    }
gen.Xs = function(ss, lambda, psi){
    J = length(lambda)
    Z = rnorm(J)
    parts = matrix(rep(sqrt(lambda), each = length(ss)), ncol = J) *
        psi(ss, 1:J) *
        matrix(rep(Z, each = length(ss)), ncol = J)
    out = apply(parts, 1, sum)
    return (out)
    }

pdf("./figs/no3.pdf", width = 12, height = 10)
par(mfrow = c(3,2), mar = c(4.1, 2.1, 2.1, 1.1))
for (j in c(5, 10, 20, 50)){
    set.seed(1)
    tmp = get.lambda.psi(J = j, L = 2*pi)
    lambda = tmp$lambda
    psi = tmp$psi
    ss = seq(-2*L, 2*L, length = 200)
    plot(ss, gen.Xs(ss, lambda, psi), type='l', ylim = c(-3, 3),
        bty = 'n', main = paste0("J = ",j), xlab = "", ylab = "")
    }
plot(ss, grf(n = length(ss), grid = cbind(ss, 0), cov.model = "matern", kappa = 1/2,
    cov.pars = c(1, 1))$data, type = 'l', ylim = c(-3, 3), bty = 'n', main = "Exact")
plot(ss, grf(n = length(ss), grid = cbind(ss, 0), cov.model = "matern", kappa = 1/2,
    cov.pars = c(1, 1))$data, type = 'l', ylim = c(-3, 3), bty = 'n', main = "Exact")
par(mfrow = c(1,1), mar = c(5.1, 4.1, 4.1, 2.1))
dev.off()

### Page 13 approx
L = 2*pi
f = function(k) 1/(1 + k^2)
get.psi2 = function(s, j){
    out = matrix(0, length(s), j)
    for (i in 1:length(s))
        out[i,] = cos(1:j * pi * s[i] / (2*L))
    return (out)
    }

pdf("./figs/no4.pdf", width = 12, height = 10)
par(mfrow = c(3,2), mar = c(4.1, 2.1, 2.1, 1.1))
for (j in c(5, 10, 20, 50)){
    set.seed(1)
    ss = seq(-2*L, 2*L, length = 200)
    lam2 = f(1:j * pi / (2 * L))
    psi2 = get.psi2(ss, j)
    plot(ss, gen.Xs(ss, lam2, psi2), type='l', ylim = c(-3, 3),
        bty = 'n', main = paste0("J = ",j), xlab = "", ylab = "")
    }
plot(ss, grf(n = length(ss), grid = cbind(ss, 0), cov.model = "matern", kappa = 1/2,
    cov.pars = c(1, 1))$data, type = 'l', ylim = c(-3, 3), bty = 'n', main = "Exact")
plot(ss, grf(n = length(ss), grid = cbind(ss, 0), cov.model = "matern", kappa = 1/2,
    cov.pars = c(1, 1))$data, type = 'l', ylim = c(-3, 3), bty = 'n', main = "Exact")
par(mfrow = c(1,1), mar = c(5.1, 4.1, 4.1, 2.1))
dev.off()

### Comparison to empirical estimates for 100 realizations
gp = grf(n = length(ss), grid = cbind(ss, 0), cov.model = "matern", kappa = 1/2,
    cov.pars = c(1, 1), nsim = 100)$data

E = svd(gp)

J = 50
E$d = E$d[1:J]^2
E$u = E$u[,1:J]
psi1 = psi(ss, 1:J)


pdf("./figs/eval.pdf", width = 9, height = 9)
par(mfrow = c(1,1), mar = c(5.1, 5.1, 4.1, 2.1))
plot(1:50, lambda/max(lambda), pch = 16, bty = 'n', main = "(scaled) Eigenvalues",
    ylab = expression(lambda), xlab = "", cex.main = 2, cex.lab = 2)
points(1:50+0.1, lam2/max(lam2), pch = 16, col = 'blue')
points(1:50-0.1, E$d/max(E$d), pch = 16, col = 'green')
text(24+8, 1, "Max Eigenvalue", cex = 2, adj = 0)
text(24+8, 0.95, "K-L", cex = 1.8, col = 'black', adj = 0)
text(24+8, 0.91, "Spectral", cex = 1.8, col = 'blue', adj = 0)
text(24+8, 0.87, "Empirical", cex = 1.8, col = 'green', adj = 0)
text(34+8, 0.95, round(max(lambda), 2), cex = 1.8, col = 'black', adj = 0)
text(34+8, 0.91, round(max(lam2), 2), cex = 1.8, col = 'blue', adj = 0)
text(34+8, 0.87, round(max(E$d), 2), cex = 1.8, col = 'green', adj = 0)
par(mfrow = c(1,1), mar = c(5.1, 4.1, 4.1, 2.1))
dev.off()


pdf("./figs/evec.pdf", width = 12, height = 16)
par(mfrow = c(3,2), mar = c(5.1, 5.1, 4.1, 2.1))
for (i in 1:6){
    plot(ss, psi1[,i], ylim = c(-1, 1), type='l', bty = 'n', main = paste0("Eigenvector ", i),
        xlab = "", ylab = expression(psi), cex.lab = 2, cex.main = 2)
    lines(ss, psi2[,i], col = 'blue')
    lines(ss, E$u[,i], col = 'green')
    }
par(mfrow = c(1,1), mar = c(5.1, 4.1, 4.1, 2.1))
dev.off()
