set.seed(2)
n = 1000
x = rnorm(n)
sig = sqrt(0.01)
y = 0.3 + 0.4*x + 0.5*sin(2.7*x) + 1.1/(1+x^2) + rnorm(n, 0, sig)
points(x, y)

m = 100
z = seq(min(x) + 0.01, max(x) - 0.01, length = m)
