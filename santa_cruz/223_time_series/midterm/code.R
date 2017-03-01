dat = read.csv("./googletrendsUCSC.csv")[,2]
times = read.csv("./googletrendsUCSC.csv")[,1]
as.numeric(times)

y = diff(dat)

plot(y, type='l', bty = 'n')
# axis(2)
# ind = round(seq(1, length(y), length = 6))
# axis(1, at = ind, labels = times[ind])
