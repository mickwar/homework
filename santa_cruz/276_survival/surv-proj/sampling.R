get_samples = function(nmcmc = 10000, nburn = 10000, window = 200, k,
    chain_init, sig_eta, sig_beta, sig_alpha, display = 1000){

    require(MASS)

    trailing = function(x, digits = 4)
        formatC(x, digits=digits, format="f")
    nice_time = function(seconds){
        # floor() or round() would work as well
        seconds = ceiling(seconds)
        days = seconds %/% (60*60*24)
        seconds = seconds %% (60*60*24)
        hours = seconds %/% (60*60)
        seconds = seconds %% (60*60)
        minutes = seconds %/% (60)
        seconds = seconds %% (60)
        out = ""
        if (days > 0)
            out = paste0(out, days, "d ", hours, "h ", minutes, "m ", seconds, "s")
        if (days == 0 && hours > 0)
            out = paste0(out, hours, "h ", minutes, "m ", seconds, "s")
        if (days == 0 && hours == 0 && minutes > 0)
            out = paste0(out, minutes, "m ", seconds, "s")
        if (days == 0 && hours == 0 && minutes == 0)
            out = paste0(out, seconds, "s")
        return (out)
        }

    calc_post_eta = function(eta, w)
        return (sum(dgamma(w, eta, eta, log = TRUE)) + dgamma(eta, 0.001, 0.001, log = TRUE))

    calc_post_beta = function(beta, time, cluster, x, nu, gamma, alpha, w)
        return (sum(x[nu == 1,] %*% beta) - 
            sum(gamma*w[cluster]*(time^alpha)*exp(x %*% beta)) +
            sum(dnorm(beta, 0, sqrt(10^3), log = TRUE)))

    calc_post_alpha = function(alpha, time, cluster, x, nu, beta, gamma, w)
        return (sum(log(alpha) + alpha*log(time[nu == 1])) -
            sum(gamma*w[cluster]*(time^alpha)*exp(x %*% beta)) +
            dgamma(alpha, 0.001, 0.001, log = TRUE))


    # Try to get what nparam should be
    nparam = 38 + 5

    if (missing(k))
        k = window / 50

    if (missing(chain_init))
        chain_init = runif(nparam)

    if (missing(sig_eta))
        sig_eta = 0.1

    if (missing(sig_beta))
        sig_beta = 0.1*diag(2)

    if (missing(sig_alpha))
        sig_alpha = 0.1

    params = matrix(0, nburn + nmcmc, nparam)
    accept_eta = double(nburn + nmcmc)
    accept_beta = double(nburn + nmcmc)
    accept_alpha = double(nburn + nmcmc)
    params[1,] = chain_init

    cluster = dat[,1]
    time = dat[,2]
    nu = dat[,3]
    snu = apply(matrix(nu, ncol = 2), 1, sum)
    x = dat[,4:5]

    begin_time = as.numeric(Sys.time())
    for (i in 2:(nburn + nmcmc)){
        # Display loop
        if (floor(i/display) == i/display && display > 0){
            curr_time = as.numeric(Sys.time()) - begin_time
            cat("\r   ", i, " / ", nburn+nmcmc, " -- ",
                trailing(100*i/(nburn+nmcmc), 2),"% -- Remaining: ",
                nice_time(curr_time*(nmcmc+nburn-i)/(i-1)), "            ", sep = "")
            }

        # Current state is previous state
        params[i,] = params[i-1,]

        eta = params[i, 1]
        beta = params[i, 2:3]
        gamma = params[i, 4]
        alpha = params[i, 5]
        w = params[i, 6:43]

        # Update gamma
        exb = exp(x %*% beta)
        gamma = rgamma(1, 0.001 + sum(snu),
            0.001 + sum(w[cluster] * (time^alpha) * exb))

        # Update w
        w = rgamma(38, eta + snu,
            eta + gamma * apply(matrix(time^alpha * exb, ncol = 2), 1, sum))


        # Update eta
        cand_eta = rnorm(1, eta, sig_eta)
        if (cand_eta > 0){
            curr_post = calc_post_eta(eta, w)
            cand_post = calc_post_eta(cand_eta, w)
            if (log(runif(1)) <= cand_post - curr_post){
                eta = cand_eta
                accept_eta[i] = 1               
                }
            }

        # Update beta
        cand_beta = mvrnorm(1, beta, sig_beta)
        if (cand_eta > 0){
            curr_post = calc_post_beta(beta, time, cluster, x, nu, gamma, alpha, w)
            cand_post = calc_post_beta(cand_beta, time, cluster, x, nu, gamma, alpha, w)
            if (log(runif(1)) <= cand_post - curr_post){
                beta = cand_beta
                accept_beta[i] = 1               
                }
            }

        # Update alpha
        cand_alpha = rnorm(1, alpha, sig_alpha)
        if (cand_eta > 0){
            curr_post = calc_post_alpha(alpha, time, cluster, x, nu, beta, gamma, w)
            cand_post = calc_post_alpha(cand_alpha, time, cluster, x, nu, beta, gamma, w)
            if (log(runif(1)) <= cand_post - curr_post){
                alpha = cand_alpha
                accept_alpha[i] = 1               
                }
            }

        # Replace new state in params matrix
        params[i, 1] = eta
        params[i, 2:3] = beta
        params[i, 4] = gamma
        params[i, 5] = alpha
        params[i, 6:43] = w
        
        # Tune candidate sigmas
        if ((floor(i/window) == i/window) && (i <= nburn)){
            sig_eta = sig_eta * autotune(mean(accept_eta[(i-window+1):i]),
                target = 0.234, k = k)
            sig_beta = autotune(mean(accept_beta[(i-window+1):i]),
                target = 0.234, k = k) * (sig_beta + window *
                var(params[(i-window+1):i, 2:3]) / i)
            sig_alpha = sig_alpha * autotune(mean(accept_alpha[(i-window+1):i]),
                target = 0.234, k = k)
            }


        # End display
        if (i == (nburn + nmcmc) && display > 0){
            curr_time = as.numeric(Sys.time()) - begin_time
            cat("\r   ", i, " / ", nburn+nmcmc, " -- ",
                trailing(100, 2),"% -- Elapsed: ",
                nice_time(curr_time), "            \n", sep = "")
            }
        }

    # Discard the burn-in
    params = tail(params, nmcmc)
    accept = cbind(accept_eta, accept_beta, accept_alpha)
    accept = tail(accept, nmcmc)

    return (list(
        "params" = params,
        "accept" = accept,
        "sig_eta" = sig_eta,
        "sig_beta" = sig_beta,
        "sig_alpha" = sig_alpha))
    }
