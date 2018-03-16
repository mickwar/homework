library(MASS)
library(mwBASE)

DATA_DIR = "~/files/tmp/268_data/"

### Inverse gamma density
digamma = function(x, a, b, log = TRUE){
    out = a*log(b) -lgamma(a)- (a+1)*log(x) - b/x
    if (log)
        return (out)
    return (exp(out))
    }
dmvnorm = function(x, mu, sigma, log = TRUE){
    out = -length(x)/2*log(2*pi) - 0.5*determinant(sigma)$modulus[1] -
        0.5*t(x-mu)%*%solve(sigma)%*%(x-mu)
    if (log)
        return (out)
    return (exp(out))
    }
### Normal Inverse Gamma density
dnig = function(beta, sig2, mu, lambda, a, b, log = TRUE){
    # lambda is precision
    k = length(mu)
    U = chol(lambda)
    Uinv = t(backsolve(U, diag(k)))

    t(beta-mu) %*% lambda %*% (beta-mu)
    t(beta-mu) %*% t(U) %*% U %*% (beta-mu)
    sum((U %*% (beta-mu))^2)
    determinant(U)
    sum(log(diag(U)))

    t(beta-mu) %*% solve(lambda) %*% (beta-mu)
    t(beta-mu) %*% t(Uinv) %*% Uinv %*% (beta-mu)
    sum((Uinv %*% (beta-mu))^2)
    determinant(Uinv)
    sum(log(diag(Uinv)))

    out = -k/2*log(2*pi) + a*log(b) - lgamma(a) -(a+k/2+1)*log(sig2) -
        (1/sig2)*(b+1/2*sum((U %*% (beta - mu))^2)) -sum(log(diag(Uinv)))
    if (log)
        return (out)
    return (exp(out))
    }
### Samples for NIG
rnig = function(n, mu, lambda, a, b){
    k = length(mu)
    U = backsolve(chol(lambda), diag(k))

    sigma2 = 1/rgamma(n, a, b)
#   mu + L %*% rnorm(k, 0, 1)
    tapply(sigma2, 1:n, function(s) t(mu + (sqrt(s)*U) %*% rnorm(k)))
    beta = t(sapply(sigma2, function(s) mu + (sqrt(s)*U) %*% rnorm(k)))
    return (list("beta"=beta, "sigma2"=sigma2))
    }
### Random Inverse Gaussian
rigauss = function(n, mu, lambda){
    nu = rnorm(n)
    y = nu^2
    x = mu + y*mu^2/(2*lambda) - mu/(2*lambda) *
        sqrt(4*mu*lambda*y + mu^2*y^2)
    # account of when mu is large (becomes a levy) (happens when x = 0)
    x = ifelse(x < 0, 0, x)
    ind = which(x == 0)
    k = length(ind)

    z = rnorm(n-k)
    if (k > 0){
        x[ind] = lambda[ind] / qnorm(1 - runif(k)/2)^2
        out = x
        out[-ind] = ifelse(z <= mu[-ind] / (mu[-ind] + x[-ind]),
            x[-ind], mu[-ind]^2 / x[-ind])
    } else {
        out = ifelse(z <= mu / (mu + x), x, mu^2 / x)
        }
    return (out)
    }
### NIG summation
nig_sum = function(mu1, lambda1, a1, b1, mu2, lambda2, a2, b2){
    k = length(mu1)
    mu = solve(lambda1 + lambda2) %*% (lambda1 %*% mu1 + lambda2 %*% mu2)
    a = a1 + a2 + k/2
    b = b1 + b2 + 1/2*t(mu1 - mu) %*% lambda1 %*% (mu1 - mu) +
        1/2*t(mu2 - mu) %*% lambda2 %*% (mu2 - mu)
    lambda = lambda1 + lambda2
    return (list("mu"=mu, "lambda"=lambda, "a"=a, "b"=as.numeric(b)))
    }

### Generate and write the data to m files
gen_data = function(n, m = 1, k, rho, sigma, beta,
    write.out = FALSE, data_dir){
    # Generate m files containing n rows and k+1 columns
    # First column, response y
    # Remaing k columns, covariates

    # x's are generated from mean-zero multivariate normal where all
    # variables have rho correlation, variance is 1.

    require(MASS)

    ### For generating X
    make.sigma = function(k, rho){
        out = matrix(rho, k, k)
        diag(out) = 1
        return (out)
        }

    if (missing(beta))
        beta = c(seq(1, 0.1, by = -0.1), rep(0, k-10))

    sig.mat = make.sigma(k, rho)

    if (!write.out){
        if (m > 1)
            message("Note: setting m = 1 to return one X, y pair")
        m = 1
        }

    for (i in 1:m){
        X = mvrnorm(n, rep(0, k), sig.mat)
        y = rnorm(n, mean = X %*% beta, sd = sigma)
        if (write.out){
            fname = paste0(c(data_dir, "dat_", rho, "_", sigma, "_",
                rep("0", nchar(m) - nchar(i)), i, ".txt"), collapse = "")
            write.table(round(data.frame(y, X), 3), quote = FALSE, row.names = FALSE, file = fname)
            rm(X, y)
            gc()
            }
        }
    if (!write.out)
        return (list("X" = X, "y" = y))
    }

### Read in data and calculate posterior distributions
### for each data subset
read_dat_and_get_params = function(rho, sigma, data_dir){
    files = list.files(data_dir, paste0("dat_", rho, "_", sigma, "_"), full.names = TRUE)
    m = length(files)

    if (m == 0)
        stop("Data not found. Generate some data homie.")

    z = read.table(files[1], header = TRUE)
    n = NROW(z)
    k = NCOL(z) - 1
    rm(z)

    mu.s = matrix(0, k, m)
    lambda.s = array(0, c(m, k, k))
    a.s = double(m)
    b.s = double(m)

    for (i in 1:m){
        z = read.table(files[i], header = TRUE)
        y = as.vector(z[,1])
        X = as.matrix(z[,-1])

        xty = t(X) %*% y
        yty = sum(y^2)

        lambda.s[i,,] = t(X) %*% X
        lam.inv = solve(t(X) %*% X)
        mu.s[,i] = lam.inv %*% xty
        a.s[i] = (n - k) / 2
        b.s[i] = 1/2 * (yty - t(xty) %*% lam.inv %*% xty)

        rm(X, y, z, xty, yty, lam.inv)
        gc()
        }

    return (list("mu.s" = mu.s, "lambda.s" = lambda.s, "a.s" = a.s, "b.s" = b.s))
    }


neach = 10000
m = 100
k = 40
beta = c(seq(1, 0.1, by = -0.1), rep(0, k - 10))
n = neach * m

# Init
MM = array(0, c(4, 4, k))
QQ = array(0, c(4, 4, 2, k))
mse = array(0, c(4, 4, 2))
mlen = array(0, c(4,4,2))
coverage = array(0, c(4, 4))
mgamma = array(0, c(4, 10, 2))


for (i1 in 1:2){
    for (i2 in 1:2){
        sigma = c(1, 10)[i1]
        rho = c(0, 0.99)[i2]

        l = 2 * ((i1 - 1) %% 2) + i2

        ### Read in data
#       set.seed(1)
#       gen_data(neach, m, k, rho, sigma, beta,
#           write.out = TRUE, data_dir)

        tmp = read_dat_and_get_params(rho, sigma, DATA_DIR)
        mu.s = tmp$mu.s
        lambda.s = tmp$lambda.s
        a.s = tmp$a.s
        b.s = tmp$b.s
        rm(tmp)

        ### Get data for predictions
        nx = 1000
        set.seed(1)
        tmp = gen_data(nx, 1, k, rho, sigma)
        new.X = tmp$X
        new.y = tmp$y
        rm(tmp)

        ### Calculate full posterior (sum of the subsets), OLS
        # (Algorithm 1)
        mu = mu.s[,1]
        lambda = lambda.s[1,,]
        a = a.s[1]
        b = b.s[1]

        for (i in 2:m){
            tmp = nig_sum(mu, lambda, a, b, mu.s[,i], lambda.s[i,,], a.s[i], b.s[i])
            mu = tmp$mu
            lambda = tmp$lambda
            a = tmp$a
            b = tmp$b
            rm(tmp)
            }
        # Non-informative prior
        alg1 = rnig(5000, mu, lambda, a, b)


        ### LASSO
        # (Algorithm 2, start with Algorithm 1)

        # Penalty terms
        pen.vec = rep(10, k)

        # Gibbs set up
        nburn = 2000
        nmcmc = 5000
        B = nburn + nmcmc
        draws.beta = matrix(0, B, k)
        draws.sig2 = double(B)
        draws.psi = matrix(0, B, k)

        draws.sig2[1] = 1
        draws.psi[1,] = rep(1, k)

        for (i in 2:B){
            lambda.psi = diag(1/draws.psi[i-1,])

            # Get NIG parameters
            tmp = nig_sum(mu, lambda, a, b, rep(0, k), lambda.psi, 0, 0)
            mu.bar = tmp$mu
            lambda.bar = tmp$lambda
            a.bar = tmp$a
            b.bar = tmp$b

            # Update sigma^2
            draws.sig2[i] = 1/rgamma(1, a.bar, b.bar)

            # Update beta
            draws.beta[i,] = mvrnorm(1, mu.bar, draws.sig2[i] * solve(lambda.bar))

            # Update psi
            draws.psi[i,] = 1/rigauss(k, abs(pen.vec * draws.sig2[i] / draws.beta[i,]), pen.vec^2)
            }

        draws.beta = tail(draws.beta, nmcmc)
        draws.sig2 = tail(draws.sig2, nmcmc)
        draws.psi = tail(draws.psi, nmcmc)

        alg2 = list("beta"=draws.beta, "sigma2"=draws.sig2, "psi"=draws.psi)


        ### SSVS
        # (Algorithm 3)
        # Notes: It seems the spike part (d) needs to have a really small variance
        # (0.001^2) for this data set. This is probably because the observation variance
        # is included when calculated the bernoulli probability for gamma.
        # In HW1, that variance was not included. Was this a mistake from Qian
        # or Raj? I want to say there is a mistake in Qian, since including
        # sigma^2 in calculating the probabilities must be dealt with when
        # selecting d and D. But maybe not, the variance of beta.hat in linear regression
        # does include sigma^2. So...
        # How to choose reasonable d and D?
        #
        # The above comment was written after running the model on data with no 
        # correlation between covariates. When there was strong correlation (0.99),
        # then every variable was selected by SSVS, the spread of the each beta
        # was much larger, those (non-zero) betas closer to 0 had credible intervals crossing 0.
        # Basically, this method, as it is, is not sufficient when dealing with
        # correlated predictors.

        # Variance for the spike (d) and slab (D) components
        #d = 0.01^2
        #D = 100
        #D = 100
        #d = 1

        g = 100
        tau = sqrt(0.01)

        # Beta prior for w (mixture weight)
        #a.w = rep(1, k)
        #b.w = rep(1, k)

        # Fixed w
#       w = 1 / (1 + 1/sqrt(1+D) * exp(log(n)*D/(2*(1+D))))     # BIC
        w = 1 / (1 + 1/sqrt(1+D) * exp(log(n)*g/(2*(1+g))))     # BIC
        w.vec = rep(w, k)


        # Gibbs set up
        nburn = 2000
        nmcmc = 5000
        B = nburn + nmcmc
        draws.beta = matrix(0, B, k)
        draws.sig2 = double(B)
        draws.gamma = matrix(0, B, k)   # Indicators for selected variables
        #draws.w = matrix(0, B, k)       # Probabilities for indicators

        draws.sig2[1] = 1
        draws.gamma[1,] = rep(0, k)     # Start with no selected variables
        #draws.w[1,] = rep(0.5, k)

        for (i in 2:B){
            lambda.gamma = diag(1/ifelse(draws.gamma[i-1,] == 1, g^2*tau^2, tau^2))

            # Get NIG parameters
            tmp = nig_sum(mu, lambda, a, b, rep(0, k), lambda.gamma, 0, 0)
            mu.bar = tmp$mu
            lambda.bar = tmp$lambda
            a.bar = tmp$a
            b.bar = tmp$b

            # Update sigma^2
            draws.sig2[i] = 1/rgamma(1, a.bar, b.bar)

            # Update beta
            draws.beta[i,] = mvrnorm(1, mu.bar, draws.sig2[i] * solve(lambda.bar))

            # Update gamma
        #   h1 = draws.w[i-1,] * dnorm(draws.beta[i,], 0, sqrt(D))
        #   h2 = (1 - draws.w[i-1,]) * dnorm(draws.beta[i,], 0, sqrt(d))
        #   h1 = w.vec * dnorm(draws.beta[i,], 0, sqrt(D))
        #   h2 = (1 - w.vec) * dnorm(draws.beta[i,], 0, sqrt(d))
        #   h1 = w.vec * dnorm(draws.beta[i,], 0, sqrt(draws.sig2[i]*D))
        #   h2 = (1 - w.vec) * dnorm(draws.beta[i,], 0, sqrt(draws.sig2[i]*d))
            h1 = w.vec * dnorm(draws.beta[i,], 0, g*tau)
            h2 = (1 - w.vec) * dnorm(draws.beta[i,], 0, tau)
            draws.gamma[i,] = rbinom(k, 1, h1 / (h1 + h2))

            # Update w
        #   draws.w[i,] = rbeta(k, a.w + draws.gamma[i,], b.w + 1 - draws.gamma[i,]) 
            }

        draws.beta = tail(draws.beta, nmcmc)
        draws.sig2 = tail(draws.sig2, nmcmc)
        draws.gamma = tail(draws.gamma, nmcmc)
        #draws.w = tail(draws.w, nmcmc)

        #alg3 = list("beta"=draws.beta, "sigma2"=draws.sig2, "gamma"=draws.gamma, "w"=draws.w)
        alg3 = list("beta"=draws.beta, "sigma2"=draws.sig2, "gamma"=draws.gamma)


        ### MC^3
        # (Algorithm 4)
        calc.post = function(gamma, mu, lambda, a, b){
            vg = (gamma == 1)
            vo = (gamma == 0)
            kg = sum(vg)
            ko = sum(vo)

            if (kg == 0 || ko == 0)
                return (-Inf)

            sig2 = b/(a-1)

            mu.g = mu[vg]
            mu.o = mu[vo]

            lam.gg = matrix(lambda[vg, vg], kg, kg)
            lam.go = matrix(lambda[vg, vo], kg, ko)
            lam.og = matrix(lambda[vo, vg], ko, kg)
            lam.oo = matrix(lambda[vo, vo], ko, ko)

            # gamma subscripts
             mu.squig = mu.g + solve(lam.gg) %*% lam.go %*% mu.o
            lam.squig = lam.gg
        #     a.squig = a + ko / 2
              a.squig = a + (1-m) * ko / 2
              b.squig = b + 1/2 * t(mu.o) %*% (lam.oo - lam.og %*% solve(lam.gg) %*% lam.go) %*% mu.o
            ### Not sure which a.squig is correct. The paper gives a + ko / 2, but when
            ### doing the exact calculation (using only the selected variables), then
            ### it's a + (1-m) * ko / 2

            mu.prior = rep(0, kg)
        #   lam.prior = 1/D * diag(kg)
            lam.prior = 1/D * lambda[vg, vg]
            a.prior = 1
            b.prior = 1 
            
            tmp = nig_sum(mu.squig, lam.squig, a.squig, b.squig,
                mu.prior, lam.prior, a.prior, b.prior)
            mu.bar = tmp$mu
            lam.bar = tmp$lambda
            a.bar = tmp$a
            b.bar = tmp$b
            
            return (dnig(rep(0, kg), sig2, mu.prior, lam.prior, a.prior, b.prior) -
                dnig(rep(0, kg), sig2, mu.bar, lam.bar, a.bar, b.bar) +
                sum(dbinom(gamma, 1, prob = w.vec, log = TRUE)))
        #   return (-log(sig2) - dnig(rep(0, kg), sig2, mu.bar, lam.bar, a.bar, b.bar) +
        #       sum(dbinom(gamma, 1, prob = w.vec, log = TRUE)))
            }

        alg4.draws = function(gammas, mu, lambda, a, b){
            B = NROW(gammas)
            k = NCOL(gammas)
            draws.beta = matrix(0, B, k)
            draws.sig2 = matrix(0, B)
            for (j in 1:B){
                vg = (gammas[j,] == 1)
                vo = (gammas[j,] == 0)
                kg = sum(vg)
                ko = sum(vo)

                if (kg == 0 || ko == 0)
                    return (-Inf)

                sig2 = b/(a-1)

                mu.g = mu[vg]
                mu.o = mu[vo]

                lam.gg = matrix(lambda[vg, vg], kg, kg)
                lam.go = matrix(lambda[vg, vo], kg, ko)
                lam.og = matrix(lambda[vo, vg], ko, kg)
                lam.oo = matrix(lambda[vo, vo], ko, ko)

            #   # gamma subscripts
                 mu.squig = mu.g + solve(lam.gg) %*% lam.go %*% mu.o
                lam.squig = lam.gg
            #     a.squig = a + ko / 2
                  a.squig = a + (1-m) * ko / 2
                  b.squig = b + 1/2 * t(mu.o) %*% (lam.oo - lam.og %*% solve(lam.gg) %*% lam.go) %*% mu.o

                mu.prior = rep(0, kg)
            #   lam.prior = 1/D * diag(kg)
                lam.prior = 1/D * lambda[vg, vg]
                a.prior = 1
                b.prior = 1 

                tmp = nig_sum(mu.squig, lam.squig, a.squig, b.squig,
                    mu.prior, lam.prior, a.prior, b.prior)
                mu.bar = tmp$mu
                lam.bar = tmp$lambda
                a.bar = tmp$a
                b.bar = tmp$b

                tmp = rnig(1, mu.bar, lam.bar, a.bar, b.bar)
                draws.beta[j,vg] = tmp$beta
                draws.sig2[j] = tmp$sigma2
                }
            return (list("beta"=draws.beta, "sigma2"=draws.sig2, "gamma"=gammas))
            }
            
        # d = 0.01^2
        D = 100

        # Gibbs set up
        nburn = 5000
        nmcmc = 20000
        B = nburn + nmcmc
        draws.gamma = matrix(0, B, k)   # Indicators for selected variables
        accept = double(B)

        draws.gamma[1,] = rbinom(k, 1, prob = 0.5)     # Start with randomly selected variables

        #w = 0.5
        #w = 1 / (1 + 1/sqrt(1+c) * exp(2*c/(2*(1+c))))         # AIC
        w = 1 / (1 + 1/sqrt(1+D) * exp(log(n)*D/(2*(1+D))))     # BIC
        w.vec = rep(w, k)

        curr.post = calc.post(draws.gamma[1,], mu, lambda, a, b)

        for (i in 2:B){
            draws.gamma[i,] = draws.gamma[i-1,]
            if (i %% 100 == 0)
                cat(i, "/", B, "\r")

            cand.gamma = draws.gamma[i-1,]
            tmp.ind = sample(k, 1)
            cand.gamma[tmp.ind] = 1 - cand.gamma[tmp.ind]

            cand.post = calc.post(cand.gamma, mu, lambda, a, b)

            if (log(runif(1)) <= cand.post - curr.post){
                draws.gamma[i,] = cand.gamma
                curr.post = cand.post
                accept[i] = 1
                }

            if (i == B)
                cat("\n")
            }

        draws.gamma = tail(draws.gamma, nmcmc)
        accept = tail(accept, nmcmc)

        alg4 = alg4.draws(draws.gamma, mu, lambda, a, b)

        ### Statistics for posterior betas
        MM[l,1,] = colMeans(alg1$beta)
        QQ[l,1,,] = apply(alg1$beta, 2, hpd_mult, force_uni = TRUE)
        MM[l,2,] = colMeans(alg2$beta)
        QQ[l,2,,] = apply(alg2$beta, 2, hpd_mult, force_uni = TRUE)
        MM[l,3,] = colMeans(alg3$beta)
        QQ[l,3,,] = apply(alg3$beta, 2, hpd_mult, force_uni = TRUE)
        MM[l,4,] = colMeans(alg4$beta)
        QQ[l,4,,] = apply(alg4$beta, 2, hpd_mult, force_uni = TRUE)

        ### MSE
        # MSE for beta
        mse[l,1,1] = mean((alg1$beta - t(matrix(beta, k, nrow(alg1$beta))))^2)
        mse[l,2,1] = mean((alg2$beta - t(matrix(beta, k, nrow(alg2$beta))))^2)
        mse[l,3,1] = mean((alg3$beta - t(matrix(beta, k, nrow(alg3$beta))))^2)
        mse[l,4,1] = mean((alg4$beta - t(matrix(beta, k, nrow(alg4$beta))))^2)

        # Prediction MSE
        mse[l,1,2] = mean((new.y - rowMeans(new.X %*% t(alg1$beta)))^2)
        mse[l,2,2] = mean((new.y - rowMeans(new.X %*% t(alg2$beta)))^2)
        mse[l,3,2] = mean((new.y - rowMeans(new.X %*% t(alg3$beta)))^2)
        mse[l,4,2] = mean((new.y - rowMeans(new.X %*% t(alg4$beta)))^2)

        ### Mean posterior interval lengths
        ind1 = 1:10
        ind2 = 11:k

        mlen[l,1,1] = mean(apply(apply(alg1$beta, 2, quantile, c(0.025, 0.975)), 2, diff)[ind1])
        mlen[l,1,2] = mean(apply(apply(alg1$beta, 2, quantile, c(0.025, 0.975)), 2, diff)[ind2])
        mlen[l,2,1] = mean(apply(apply(alg2$beta, 2, quantile, c(0.025, 0.975)), 2, diff)[ind1])
        mlen[l,2,2] = mean(apply(apply(alg2$beta, 2, quantile, c(0.025, 0.975)), 2, diff)[ind2])
        mlen[l,3,1] = mean(apply(apply(alg3$beta, 2, quantile, c(0.025, 0.975)), 2, diff)[ind1])
        mlen[l,3,2] = mean(apply(apply(alg3$beta, 2, quantile, c(0.025, 0.975)), 2, diff)[ind2])
        mlen[l,4,1] = mean(apply(apply(alg4$beta, 2, quantile, c(0.025, 0.975)), 2, diff)[ind1])
        mlen[l,4,2] = mean(apply(apply(alg4$beta, 2, quantile, c(0.025, 0.975)), 2, diff)[ind2])

        ### Coverage
        tx = nrow(alg1$beta)       
        pred.y = alg1$beta %*% t(new.X) + matrix(rnorm(nx * tx, 0, sqrt(alg1$sigma2)), tx, nx)
        qq = apply(pred.y, 2, quantile, c(0.025, 0.975))
        coverage[l,1] = mean(qq[1,] < new.y & qq[2,] > new.y)

        tx = nrow(alg2$beta)       
        pred.y = alg2$beta %*% t(new.X) + matrix(rnorm(nx * tx, 0, sqrt(alg2$sigma2)), tx, nx)
        qq = apply(pred.y, 2, quantile, c(0.025, 0.975))
        coverage[l,2] = mean(qq[1,] < new.y & qq[2,] > new.y)

        tx = nrow(alg3$beta)       
        pred.y = alg3$beta %*% t(new.X) + matrix(rnorm(nx * tx, 0, sqrt(alg3$sigma2)), tx, nx)
        qq = apply(pred.y, 2, quantile, c(0.025, 0.975))
        coverage[l,3] = mean(qq[1,] < new.y & qq[2,] > new.y)

        tx = nrow(alg4$beta)       
        pred.y = alg4$beta %*% t(new.X) + matrix(rnorm(nx * tx, 0, sqrt(alg4$sigma2)), tx, nx)
        qq = apply(pred.y, 2, quantile, c(0.025, 0.975))
        coverage[l,4] = mean(qq[1,] < new.y & qq[2,] > new.y)

        matplot(t(MM[l,,]), type = 'l', xlab = "Beta Index", ylab = "Coefficient",
            axes = FALSE, main = "Regression coefficient estimates", lwd = 2)
        axis(1); axis(2)
        legend("topright", bty = 'n', legend = c("OLS", "LASSO", "SSVS", "MC3", "Truth"),
            col = c(1:4, 1), lty = c(1:4, NA), cex = 1.5, lwd = 3, pch = c(rep(NA, 4), 4))
        points(beta, pch = 4, lwd = 3)

        ### Mean gammas
        mgamma[l,,] = head(cbind(colMeans(alg3$gamma), colMeans(alg4$gamma)), 10)
        }
    }

library(xtable)
xtable(log(mse[,,1]))

big_mlen = matrix(0, 8, 4)
big_mlen[1,] = t(mlen[1,,1])
big_mlen[2,] = t(mlen[2,,1])
big_mlen[3,] = t(mlen[3,,1])
big_mlen[4,] = t(mlen[4,,1])
big_mlen[5,] = t(mlen[1,,2])
big_mlen[6,] = t(mlen[2,,2])
big_mlen[7,] = t(mlen[3,,2])
big_mlen[8,] = t(mlen[4,,2])
xtable(big_mlen, digits = 3)

coverage
xtable(coverage)

big_mgamma = matrix(0, 8, 10)
big_mgamma[1:2,] = t(mgamma[1,,])
big_mgamma[3:4,] = t(mgamma[2,,])
big_mgamma[5:6,] = t(mgamma[3,,])
big_mgamma[7:8,] = t(mgamma[4,,])
xtable(big_mgamma)


pdf("figs/four_estimates.pdf", width = 12, height = 12)
par(mfrow = c(2,2))
for (i1 in 1:2){
    for (i2 in 1:2){
        sigma = c(1, 10)[i1]
        rho = c(0, 0.99)[i2]

        l = 2 * ((i1 - 1) %% 2) + i2
        matplot(t(MM[l,,]), type = 'l', xlab = "Beta Index", ylab = "Coefficient",
            axes = FALSE, main = bquote(rho ~ "=" ~ .(rho) ~~~~~ sigma ~ "=" ~ .(sigma)),
            lwd = 2)
        axis(1); axis(2)
        legend("topright", bty = 'n', legend = c("OLS", "LASSO", "SSVS", "MC3", "Truth"),
            col = c(1:4, 1), lty = c(1:4, NA), cex = 1.5, lwd = 3, pch = c(rep(NA, 4), 4))
        points(beta, pch = 4, lwd = 3)
        }
    }
dev.off()


# plot(new.y, colMeans(pred.y), ylim = range(qq))
# segments(new.y, qq[1,], new.y, qq[2,])
# abline(0,1)
# 
# plot(density(alg1$sigma2))




# pdf("slides/figs/estimates.pdf", width = 8, height = 8)
# matplot(t(MM[4,,]), type = 'l', xlab = "Beta Index", ylab = "Coefficient",
#     axes = FALSE, main = "Regression coefficient estimates", lwd = 2)
# axis(1); axis(2)
# legend("topright", bty = 'n', legend = c("OLS", "LASSO", "SSVS", "MC3", "Truth"),
#     col = c(1:4, 1), lty = c(1:4, NA), cex = 1.5, lwd = 3, pch = c(rep(NA, 4), 4))
# points(beta, pch = 4, lwd = 3)
# dev.off()
# 
# ### Gammas
# library(xtable)
# xtable(head(cbind(colMeans(alg3$gamma), colMeans(alg4$gamma)), 10))
# xtable(t(head(cbind(colMeans(alg3$gamma), colMeans(alg4$gamma)), 10)))
