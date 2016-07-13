# clear variables and close windows
rm(list = ls(all = TRUE))
graphics.off()

# install and load packages
install.packages("pracma")
library(pracma)
install.packages("timeSeries")
library(timeSeries)
install.packages("graphics")
library(graphics)
# functions
localHaltonSingleNumber = function(n, b) {
    n0 = n
    hn = 0
    f = 1/b
    while (n0 > 0) {
        n1 = floor(n0/b)
        r = n0 - n1 * b
        hn = hn + f * r
        f = f/b
        n0 = n1
    }
    hn
}
HaltonSequence = function(n, b) {
    hs = rep(0, n)
    for (idx in 1:n) {
        hs[idx] = localHaltonSingleNumber(idx, b)
    }
    hs
}
BoxMueller = function(hs1, hs2) {
    R = sqrt(-2 * log(hs1))
    Theta = 2 * pi * hs2
    P = R * sin(Theta)
    Q = R * cos(Theta)
    sequence = c(P, Q)
    sequence
}

AssetPathsHalton = function(S0, mu, sig, dt, steps, nsims) {
    pVec = primes(1e+05)
    a = 2 * steps
    bases = pVec[1:a]
    nu = mu - sig * sig/2
    epsilon = matrix(1, steps, nsims)
    for (idx in 1:steps) {
        epsH1 = HaltonSequence(nsims/2, bases[idx])
        epsH2 = HaltonSequence(nsims/2, bases[steps + idx])
        c = BoxMueller(epsH1, epsH2)
        epsilon[idx, ] = t(c)
    }
    b = exp(nu * dt + sig * sqrt(dt) * epsilon)
    f = matrix(1, steps, nsims)
    for (ids in 1:nsims) {
        e = b[, ids]
        f[, ids] = cumprod(e)
    }
    g = matrix(1, 1, nsims)
    h = rbind(g, f)
    S = S0 * h
    S
}
#### parameter setting
S0 = 50
# Price of underlying today
mu = 0.04
# expected return
sig = 0.1
# expected volatility
dt = 1/365
# time steps
steps = 50
# days to expiry
nsims = 1000
# Number of simulated paths
X = 50
# Strike at expiry
r = 0.03
# Risk free rate
T = dt * steps
# Plot the asset paths
S = AssetPathsHalton(S0, mu, sig, dt, steps, nsims)
matplot(S, lwd = 2, type = "l", lty = 1, ylim = c(min(colMins(S)), max(colMaxs(S))), col = 1:nsims, main = "Simulated Asset Paths", 
    xlab = "Days of holding the option", ylab = "Asset Price ")

# Price a standard Asian Put and Call option
i = X - colMeans(S)
PutPayoffT = rep(0, nsims)
CallPayoffT = rep(0, nsims)
for (idx in 1:nsims) {
    PutPayoffT[idx] = max(i[idx], 0)
    CallPayoffT[idx] = max(-i[idx], 0)
}
putPrice = mean(PutPayoffT) * exp(-r * T)
callPrice = mean(CallPayoffT) * exp(-r * T)
print(putPrice)
print(callPrice)







