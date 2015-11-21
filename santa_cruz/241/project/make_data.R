strain_rate = c(3.9e-3, 3e-3, 1e-9, 1e-7, 1.3e-3, 2.8e-3, 2.6e-3, 1e-9, 1.8e-3, 2.2e-3)
temperature = c(1073, 1273, 298, 298, 298, 473, 673, 77, 77, 873)

length(temperature)
length(strain_rate)

make = function(a, b, root, sd, seed){
    set.seed(seed)
    A = cbind(1, a, a^root)
    y = solve(A) %*% b
    x = seq(a[1], a[2], length = 100)
    X = cbind(1, x, x^root)
    return (list("strain_rate" = strain_rate[seed], "temperature" = temperature[seed],
        "x"=x, "y" = X %*% y + rnorm(100, 0, sd)))
    }

cols = rainbow(10)

dat = NULL
dat[[1]] = make(c(0.023, 0.24, 0.1), c(0.0014, 0.0035, 0.0028), 0.5, 0.00007, 1)
dat[[2]] = make(c(0.023, 0.19, 0.1), c(0.0017, 0.0030, 0.0025), 0.8, 0.00004, 2)
dat[[3]] = make(c(0.023, 0.18, 0.1), c(0.0020, 0.0035, 0.0030), 0.7, 0.00001, 3)
dat[[4]] = make(c(0.023, 0.18, 0.1), c(0.0030, 0.0039, 0.0036), 0.2, 0.00001, 4)
dat[[5]] = make(c(0.023, 0.13, 0.07), c(0.0057, 0.0070, 0.0065), 0.5, 0.00005, 5)
dat[[6]] = make(c(0.023, 0.26, 0.1), c(0.0036, 0.0059, 0.0047), 1.5, 0.00003, 6)
dat[[7]] = make(c(0.023, 0.26, 0.1), c(0.0023, 0.0046, 0.0033), 1.5, 0.00004, 7)
dat[[8]] = make(c(0.023, 0.17, 0.1), c(0.0088, 0.0101, 0.0093), 0.5, 0.00001, 8)
set.seed(8); dat[[8]]$y = -cos(dat[[8]]$x*19)/1500 + 0.0094 + rnorm(100, 0, 0.000005)
dat[[9]] = make(c(0.023, 0.16, 0.1), c(0.0093, 0.00962, 0.0094), 0.5, 0.00005, 9)
dat[[10]] = make(c(0.023, 0.22, 0.1), c(0.0020, 0.0037, 0.0030), 0.5, 0.00005, 10)

plot(0, type='n', xlim = c(0, 0.27), ylim = c(0, 0.012), xlab = "Plastic Strain",
    ylab = "Stress")
for (i in 1:length(dat))
    lines(dat[[i]]$x, dat[[i]]$y, type='l', col = cols[i], lwd = 2)

save("dat", file = "./chen_gray_fake.RData")

