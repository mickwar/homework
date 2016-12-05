theta.s = c(2.037, 6.111, 30.555, 2.037)

# Pr(CR | T = S)
(theta.s[1] + theta.s[2]) / sum(theta.s)

# Pr(TOX | T = S)
(theta.s[2] + theta.s[4]) / sum(theta.s)

theta.e = c(2.037, 6.111, 30.555, 2.037)
theta.e = theta.e * 4 / sum(theta.e)

Nmax = 75
k = 4
n = 0
while (n <= Nmax){

    }
