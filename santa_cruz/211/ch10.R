a = function(t, v0)
    cbind(rep(0, length(t)), rep(-9.8, length(t)))
v = function(t, v0)
    cbind(0, a(t)[,2]*t)+matrix(rep(v0, length(t)), length(t), 2, byrow = TRUE)
r = function(t, v0)
    0.5*a(t)*t^2 + matrix(rep(v0, length(t)), length(t), 2, byrow = TRUE)*t


tseq = seq(0, 1, length = 1000)
x = r(tseq)
x = x[1:(which.max(x[,2] < 0)-1),]
xseq = tseq[1:NROW(x)]

plot(x, type = 'l')

fplot = function(v0, tseq, ...){
    if (missing(tseq))
        tseq = seq(0, 1, length = 1000)
    x = r(tseq, v0)
    if (any(x[,2] < 0))
        x = x[1:(which.max(x[,2] < 0)-1),]
    lines(x, ...)
    }

plot(0, type='n', xlim = c(0, 5), ylim = c(0, 3))
fplot(c(4, 7), tseq = seq(0, 5, length = 10000))

