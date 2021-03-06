#c(0.8893,
#0.8893 + 1.2564,
#0.8893 + 1.9860,
#0.8893 + 1.9860 + 1.2564 - 1.6631,
#0.8893 + 1.1889,
#0.8893 + 1.1889 + 1.2564 - 1.6130)

# Model y = XB + eps, eps ~ N(0, sig2)
# x_i = (x_i1, x_i2) (no intercept)
y = c(82, 79, 74, 83, 80, 81, 84, 81)
X = cbind(c(10, 9, 9, 11, 11, 10, 10, 12), c(15, 14, 13, 15, 14, 14, 16, 13))
n = length(y)
r = NCOL(X)

# Parameter estimates
(betahat = solve(t(X) %*% X) %*% t(X) %*% y)
(sig2hat = t(y - X %*% betahat) %*% (y - X %*% betahat) / (n - r))

summary(lm(y ~ 0 + X)) # Built-in check

# Each are estimable
lambda1 = c(1, 0)
lambda2 = c(1, 1)

# 95% confidence intervals
t(lambda1) %*% betahat + c(-1, 1) * qt(0.975, n - r) *
    sqrt(sig2hat * t(lambda1) %*% solve(t(X) %*% X) %*% lambda1)
t(lambda2) %*% betahat + c(-1, 1) * qt(0.975, n - r) *
    sqrt(sig2hat * t(lambda2) %*% solve(t(X) %*% X) %*% lambda2)

# t-test for H_0: beta_2 = 3
lambda3 = c(0, 1)
2*pt(abs((t(lambda3) %*% betahat - 3) / sqrt(sig2hat * t(lambda3) %*% solve(t(X) %*% X) %*% lambda3)), n-r, lower.tail = FALSE)

# t-test for H_0: beta_1 - beta_2 = 0
lambda4 = c(1, -1)
2*pt(abs((t(lambda4) %*% betahat - 0) / sqrt(sig2hat * t(lambda4) %*% solve(t(X) %*% X) %*% lambda4)), n-r, lower.tail = FALSE)
((t(lambda4) %*% betahat - 0) / sqrt(sig2hat * t(lambda4) %*% solve(t(X) %*% X) %*% lambda4))




Z = cbind(1, X)
summary(lm(y ~ 0 + Z))
q = NCOL(Z)
bhat = solve(t(Z) %*% Z) %*% t(Z) %*% y

Pz = Z %*% solve(t(Z) %*% Z) %*% t(Z)
P1 = matrix(1/n, n, n)
sd2 = function(x) sqrt(mean((x - mean(x))^2))

(t(y) %*% (Pz - P1) %*% y) / (t(y) %*% (diag(n) - P1) %*% y) # R^2
1-(t(y) %*% (diag(n) - Pz) %*% y) / (t(y) %*% (diag(n) - P1) %*% y) # R^2
cor(y, Pz %*% y)^2
cor(y, predict(lm(y ~ 1 + X)))^2
mean((y - mean(y)) / sd2(y) * (yhat - mean(yhat)) / sd2(yhat))^2 # Correlation squared

(t(y - mean(y)) %*% (yhat - mean(yhat)) / (sd2(y) * sd2(yhat)) / n ) ^2

1-t(y - yhat) %*% (y - yhat) / (t(y) %*% y - n *mean(y)^2)


(t(y - mean(y)) %*% ((Pz %*% y) - mean(Pz %*% y))) / (var(y) * var(Pz %*% y) * (n-1)^2/(n-0)^2)

sum((y - (Pz %*% y))^2) / (n-q)

t(y - mean(y)) %*% (y - mean(y)) / (n - 1)

yhat = Pz %*% y
mean(Pz %*% y)
t(Pz %*% y) %*% rep(1, n) / n 

(t(y - (Pz %*% y)) %*% (y - (Pz %*% y))) / t(y) %*% y

