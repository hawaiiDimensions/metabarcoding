## script to test rough convergence of model

setwd('~/Dropbox/hawaiiDimensions/metabarcoding/poc')
source('nimble_funs.R')

load('test_mcmc.RData')

effectiveSize(mod[[1]])
effectiveSize(mod[[2]])

plot(mod[[1]][, 1], type = 'l')
points(mod[[2]][, 1], type = 'l', col = hsv(1, 1, 1, alpha = 0.5))

plot(density(mod[[1]][, 1]))
abline(v = mean(mod[[1]][, 1]))
lines(density(mod[[2]][, 1]), col = hsv(1, 1, 1, alpha = 0.5))
abline(v = mean(mod[[2]][, 1]), col = 'red')
curve(dexp(x, 0.00001), add = TRUE, col = 'blue')


plot(colMeans(mod[[2]]), effectiveSize(mod[[2]]), log = 'xy')
