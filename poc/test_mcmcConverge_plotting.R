## script to test rough convergence of model

setwd('~/Dropbox/hawaiiDimensions/metabarcoding/poc')
source('nimble_funs.R')

load('test_mcmc.RData')

effectiveSize(mod[[1]])
effectiveSize(mod[[2]])

plot(mod[[1]][-(1:100), 24], type = 'l')
points(mod[[2]][-(1:100), 24], type = 'l', col = hsv(1, 1, 1, alpha = 0.5))

plot(density(mod[[1]][-(1:100), 38]))
lines(density(mod[[2]][-(1:100), 38]), col = hsv(1, 1, 1, alpha = 0.5))
curve(dexp(x, 0.00001), add = TRUE, col = 'blue')
