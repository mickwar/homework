get_samples = function(nmcmc = 10000, nburn = 10000, window = 200, k,
    chain_init, cand_sig, display = 1000){

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

    # Try to get what nparam should be
    nparam = 38 + 5

    if (missing(k))
        k = window / 50

    if (missing(chain_init))
        chain_init = runif(nparam)

    if (missing(cand_sig))
        cand_sig = 0.1*diag(4)

    params = matrix(0, nburn + nmcmc, nparam)
    accept = double(nburn + nmcmc)
    params[1,] = chain_init

    cluster = dat[,1]
    time = dat[,2]
    snu = apply(matrix(dat[,3], ncol = 2), 1, sum)
    x = dat[,4:5]

    begin_time = as.numeric(Sys.time())
    for (i in 2:(nburn + nmcmc)){
        if (floor(i/display) == i/display && display > 0){
            curr_time = as.numeric(Sys.time()) - begin_time
            cat("\r   ", i, " / ", nburn+nmcmc, " -- ",
                trailing(100*i/(nburn+nmcmc), 2),"% -- Remaining: ",
                nice_time(curr_time*(nmcmc+nburn-i)/(i-1)), "            ", sep = "")
            }
        params[i,] = params[i-1,]
        # Update eta, beta0, beta1, alpha
        ind = c(1, 2, 3, 5)
        cand = params[i,]
        cand[ind] = mvrnorm(1, params[i-1, ind], cand_sig)
        if (all(cand[ind[c(1, 4)]] > 0)){
            curr_tval = calc_post(dat, params[i,])
            cand_tval = calc_post(dat, cand)
            if (log(runif(1)) <= cand_tval - curr_tval){
                params[i, ind] = cand[ind]
                accept[i] = 1
                }
            }
        if ((floor(i/window) == i/window) && (i <= nburn)){
            cand_sig = autotune(mean(accept[(i-window+1):i]),
                target = 0.234, k = k) * (cand_sig + window *
                var(params[(i-window+1):i, ind]) / i)
            }

        eta = params[i, 1]
        beta = params[i, 2:3]
        alpha = params[i, 5]
        w = params[i, 6:43]

        # Update gamma
        exb = exp(x %*% beta)
        params[i, 4] = rgamma(1, 0.001 + sum(snu),
            0.001 + sum(w[cluster] * (time^alpha) * exb))
        gamma = params[i, 4]

        # Update w
        params[i,6:43] = rgamma(38, eta + snu,
            eta + gamma * apply(matrix(time^alpha * exb, ncol = 2), 1, sum))

        if (i == (nburn + nmcmc) && display > 0){
            curr_time = as.numeric(Sys.time()) - begin_time
            cat("\r   ", i, " / ", nburn+nmcmc, " -- ",
                trailing(100, 2),"% -- Elapsed: ",
                nice_time(curr_time), "            \n", sep = "")
            }
        }

    # Discard the burn-in
    params = tail(params, nmcmc)
    accept = tail(accept, nmcmc)

    return (list(
        "params" = params,
        "accept" = accept,
        "cand_sig" = cand_sig))
    }
