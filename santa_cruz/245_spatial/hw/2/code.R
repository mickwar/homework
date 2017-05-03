### Karhunen-Loeve representation

L = 4*pi
J = 1000

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

par(mfrow = c(2,2), mar = c(4.1, 2.1, 2.1, 1.1))
for (j in c(5, 10, 20, 50)){
    set.seed(1)
    tmp = get.lambda.psi(J = j, L = 2*pi)
    lambda = tmp$lambda
    psi = tmp$psi
    ss = seq(-2*L, 2*L, length = 200)
    plot(ss, gen.Xs(ss, lambda, psi), type='l', ylim = c(-3, 3),
        bty = 'n', main = paste0("J = ",j), xlab = "", ylab = "")
    }
par(mfrow = c(1,1), mar = c(5.1, 4.1, 4.1, 2.1))
