library(MASS)
### glmm on epileptic data
dat = read.table("~/files/data/epileptics.txt", header = TRUE)

n = nrow(dat)
p = 4 # number of visits (and measurements)

# set up the X design matrix
X = matrix(0, n*p, 5)
X[,1] = 1                               # intercept
X[,2] = rep(log(dat$Base/4), each=p)    # Base
X[,3] = rep(dat$Trt, each=p)            # Trt
X[,4] = X[,2] * X[,3]                   # Base * Trt interaction
X[,5] = rep(log(dat$Age), each=p)       # Age
#X[,6] = rep(c(0,0,0,1), n)              # Indicator for Visit 4

# set up the Z matrix
# Z = cbind(kronecker(diag(n), rep(1, p)),
#     kronecker(diag(n), seq(-3, 3, by = 2)/10),
#     diag(n*p))
# dim(Z)
A = cbind(1, seq(-3, 3, by=2)/10)
#A = rep(1, p)
Z = kronecker(diag(n), A)



# set up Y
Y = as.numeric(t(dat[,2:5]))
Y1 = Y + 1 # for the gamma function instead of factorial

### priors
#betacovprior = 1e4 * solve(t(X) %*% X)
betacovprior = diag(NCOL(X))
betacovinv = solve(betacovprior)
betacovdet = determinant(betacovprior)$modulus[1]


calc.G = function(g){
    H = matrix(c(g[1], g[3]*sqrt(g[1]*g[2]), g[3]*sqrt(g[1]*g[2]), g[2]), 2, 2)
#   H = matrix(g[1], 1, 1)
#   Hinv = solve(H)
#   return (list("full" = kronecker(diag(n), H), "inv" = kronecker(diag(n), Hinv),
#       "H"=H))
    return (list("H"=H))
    }
calc.R = function(r){
    r[1] = 1e-6
    r[2] = 0
    K = r[1] * ((1 - r[2])*diag(p) + r[2])
#   Kinv = solve(K)
#   return (list("full" = kronecker(diag(n), K), "inv" = kronecker(diag(n), Kinv),
#       "K"=K))
    return (list("K"=K))
    }

# G = calc.G(g)
# R = calc.R(r)
# V = Z %*% calc.G(g)$full %*% t(Z) + calc.R(r)$full
# 
# V[1:4,1:4]
# W = A %*% G$H %*% t(A) + R$K
# 
# system.time(solve(V))
# system.time(kronecker(diag(n), solve(W)))



calc.post = function(params){
    # fixed effects
    at = 1
    beta = params[at:NCOL(X)]
    # parameters in G matrix (random effects)
    at = NCOL(X) + at
    g = params[at:(at+ngparam-1)]
    # parameters in R matrix (covariance of g(mu))
    at = ngparam + at
    r = params[at:(at+nrparam-1)]
    
    G = calc.G(g) # function that returns the full G and its inverse
    R = calc.R(r) # similar to calc.G

    W = A %*% G$H %*% t(A) + R$K
    Vinv = kronecker(diag(n), solve(W))
    # V = ZGZ' + R = I_n kron W

    # function g of mu (the link function), no random effects since g(mu)
    # is normally distribution (but takes into account the variances)
    gmu = X %*% beta
    mu = exp(gmu)

    # likelihood
    out = sum(Y * gmu - mu - lgamma(Y1))
    
    # priors 
    # gmu
#   out = out - n*0.5*determinant(W)$modulus[1] - (0.5*t(gmu) %*% Vinv %*% gmu)
    out = out - n*0.5*determinant(W)$modulus[1] - 0.5*t(gmu - X %*% beta) %*% Vinv %*% (gmu - X %*% beta)

    # beta
    out = out - 0.5*betacovdet - 0.5*t(beta) %*% betacovinv %*% beta

    # gparams
    out = out + dgamma(g[1], 2, 1, log = TRUE)
    out = out + dgamma(g[2], 2, 1, log = TRUE)
    out = out + dunif(g[3], -1, 1, log = TRUE)

    # rparams
#   out = out + dgamma(r[1], 2, 1, log = TRUE)
#   out = out + dbeta(0.5*r[2]+0.5, 10, 10, log = TRUE)
    
    return (out)
    }


### mcmc settings
autotune = function(accept, target = 0.25, k = 2.5)
    (1+(cosh(accept-target)-1)*(k-1)/(cosh(target-
        ceiling(accept-target))-1))^sign(accept-target)

ngparam = 3
nrparam = 2
nparams = NCOL(X) + ngparam + nrparam

lower = c(rep(-Inf, NCOL(X)), 0, 0, -1, 0, -1)
upper = c(rep(Inf, NCOL(X)), Inf, Inf, 1, Inf, 1)

nburn = 10000
nmcmc = 20000
params = matrix(0, nburn+nmcmc, nparams)

# starting values
#params[1,] = 0.5
params[1,] = c(-1.27, 0.87, -0.91, 0.33, 0.46, 1, 1, 0.5, 1, 0.5)
#params[1,] = c(-1.27, 0.87, -0.91, 0.33, 0.46, 1, 1, 0.5)
#params[1,] = old.params
#old.params = params[nrow(params),]

sigs = rep(1, nparams)
#sigs = old.sigs
#old.sigs = sigs

# if burn in has already been done (these should
# be the most up-to-date values)
#sigs = read.table("./mcmc_cand_sigmas.txt")[,1]
#params[1,] = read.table("./mcmc_init_params.txt")[,1]

# initialize log posterior value
post = calc.post(params[1,])
cand.param = params[1,]
keep.post = double(nburn+nmcmc)
keep.post[1] = post

# confidence intervals on acceptance rates?
accept = matrix(0, nburn+nmcmc, nparams)
window = 200

# mcmc loop
for (i in 2:(nburn+nmcmc)){
    cat("\rIteration",i,"/",nburn+nmcmc)
    params[i,] = params[i-1,]
    for (j in 1:nparams){
        cand = rnorm(1, params[i,j], sigs[j])
        if (cand >= lower[j] && cand <= upper[j]){
            cand.param[j] = cand
            cand.post = calc.post(cand.param)
            # check whether to accept draw or not
            if (log(runif(1)) < cand.post - post){
                post = cand.post
                params[i,j] = cand
                accept[i,j] = 1
            } else {
                cand.param[j] = params[i,j]
                }
        } else {
            cand.param[j] = params[i,j]
            }
        }
        # adjust the candidate sigma
    if (floor(i/window) == i/window && i <= nburn)
        sigs = sigs*autotune(apply(accept[(i-window+1):i,], 2,
            mean), k = max(window/50, 1.1))
    if (i == (nburn+nmcmc))
        cat("\n")
    keep.post[i] = post
    }

params = params[(nburn+1):(nburn+nmcmc),]
accept = accept[(nburn+1):(nburn+nmcmc),]
keep.post = keep.post[(nburn+1):(nburn+nmcmc)]

plot(keep.post, type='l')

apply(accept, 2, mean)

apply(params, 2, mean)
apply(params, 2, quantile, c(0.025, 0.50, 0.975))

# names = c("Intercept", "Base", "Trt", "Base*Trt", "Age", "V4",
#     "Gvar1", "Gvar2", "Gcor", "Rvar", "Rrho")
names = c("Intercept", "Base", "Trt", "Base*Trt", "Age",
    "Gvar1", "Gvar2", "Gcor", "Rvar", "Rrho")
for (j in 1:nparams){
    plot(params[,j], type='l', main=names[j])
    readline()
    }
