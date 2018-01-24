library(glmnet)

make.sigma = function(p, rho){
    out = matrix(0, p, p)
    inds = expand.grid(1:p, 1:p)
    fill = rho^apply(inds, 1, function(x) abs(diff(x)))
    out[1:p^2] = fill
    return (out)
    }

for (i1 in 1:2){
    for (i2 in 1:2){
        for (
