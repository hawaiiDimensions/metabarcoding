library(coda)
library(igraph)
library(nimble)

setwd('~/Dropbox/hawaiiDimensions/metabarcoding/poc')

diffMarkers <- read.csv('clean_diffMarkers.csv', as.is=TRUE)
pcrCycle <- read.csv('clean_pcrCycle.csv', as.is=TRUE)

## ================================
## define constants, data and inits
## ================================

## number of experiments
nexpr <- 23

## number of reads
set.seed(0)
metabConsts <- list(N = nexpr,
                   nread = rpois(nexpr, 10000))

a <- 4
b <- 1
set.seed(0)
pumpData <- list(x = rpois(ntime, pumpConsts$t * rgamma(ntime, a, b)))
pumpInits <- list(alpha = 1, beta = 1,
                  theta = rep(0.1, pumpConsts$N))


## set-up the model
metaBCode <- nimbleCode({
    for (i in 1:N){
        theta[i] ~ dgamma(alpha,beta)
        lambda[i] <- theta[i]*t[i]
        x[i] ~ dpois(lambda[i])
    }
    alpha ~ dexp(0.0001)
    beta ~ dgamma(1, 0.0001)
})



## combine model def and data, inits, etc into model
pump <- nimbleModel(code = pumpCode, name = 'pump', constants = pumpConsts,
                    data = pumpData, inits = pumpInits)

## compile model
Cpump <- compileNimble(pump)

## configure MCMC, add monitoring and compile
pumpConf <- configureMCMC(pump, print = TRUE)
pumpConf$addMonitors(c('alpha', 'beta', 'theta'))
pumpConf$addSampler(target = c('alpha', 'beta'), type = 'RW_block',
                    control = list(adaptInterval = 100))
pumpConf$setThin(10)

pumpMCMC <- buildMCMC(pumpConf)
CpumpMCMC <- compileNimble(pumpMCMC, project=pump, resetFunctions = TRUE)

## run mcmc
N <- 1000
niter <- N*pumpConf$thin
set.seed(0)
CpumpMCMC$run(niter)

## explore mcmc
samples <- as.matrix(CpumpMCMC$mvSamples)
acf(samples[, 'alpha'])
plot(samples[, 'alpha'], type='l', xlim=c(0, 100))
effectiveSize(samples[, 'alpha'])
plot(density(samples[, 'alpha']));abline(v=a)
plot(density(samples[, 'beta']));abline(v=b)


## do MLE instead
pump2 <- pump$newModel()
pumpMCEM <- buildMCEM(model = pump2, latentNodes = 'theta',
                      boxConstraints = list(list(c('alpha','beta'), c(0, Inf))))
pumpMLE <- pumpMCEM()