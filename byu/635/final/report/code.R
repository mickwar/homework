### simulation study (section 5)
library(MASS)

# functions
logit = function(p)
    log(p / (1-p))

logistic = function(x)
    1/(1+exp(-x))

# simulation parameters
K = 100
nk = 7

# fixed effects
X = matrix(1, K * nk, 4)
tl = seq(-3, 3)
X[,2] = rep(tl, K)
X[351:700,3] = 0
X[,4] = X[,2] * X[,3]
alpha = c(-2.5, 1, -1, -0.5)

# random effects
D1 = matrix(c(1,0,0,0), 2, 2, byrow=TRUE)
D2 = matrix(c(0.5,0,0,0.25), 2, 2, byrow=TRUE)
b = mvrnorm(K, c(0, 0), D1)
c1 = 16*sqrt(3)/(15*pi)


det(c1^2 * D1 %*% t(b) %*% b + diag(2))
ci = double(K)
for (i in 1:K)
    ci[i] =  (1 + c1^2 * t(b[i,]) %*% D1 %*% b[i,])^(-1/2)

b = matrix(rep(b, each = nk), K*nk, 2, byrow = F) * cbind(1, rep(tl, K))

y1 = rbinom(K*nk, 1, prob = logistic(X %*% alpha + apply(b, 1, sum)))

