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




Z = cbind(1, X)
summary(lm(y ~ 0 + Z))

Pz = Z %*% solve(t(Z) %*% Z) %*% t(Z)
P1 = matrix(1/n, n, n)

sqrt((t(y) %*% (Pz - P1) %*% y) / (t(y) %*% (diag(n) - P1) %*% y))
cor(y, Pz %*% y)
cor(y, predict(lm(y ~ 1 + X)))

