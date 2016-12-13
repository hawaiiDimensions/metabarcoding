library(nimble)

setwd('~/Dropbox/hawaiiDimensions/metabarcoding/poc')

## simulate data

## number of pools
Npool <- 100

## number of spp
Nspp <- 10

## coefficients
a <- seq(0.1, 10, length = Nspp)

## number of total reads per pool
Nreads <- round(runif(Npool, 5000, 30000))

## amount of DNA input per spp per pool (pool = rows; spp = columns)
amount_DNA <- matrix(runif(Nspp*Npool, 1, 100), nrow = Npool, ncol = Nspp)

## number of reads per spp per pool (pool = rows; spp = columns)
number_Reads <- t(sapply(1:nrow(amount_DNA), function(i) rmulti(1, Nreads[i], rdirch(1, amount_DNA[i, ]*a))))


## the model as NIMBLE code
modCode <- nimbleCode({
    for(i in 1:Npool) {
        alpha[i, 1:Nspp] <- x[i, 1:Nspp] * a[1:Nspp]
        p[i, 1:Nspp] ~ ddirch(alpha[i, 1:Nspp])
        y[i, 1:Nspp] ~ dmulti(p[i, 1:Nspp], Nreads[i])
    }
    
    ## priors
    for(j in 1:Nspp) {
        a[j] ~ dexp(0.00001)
    }

})


## model constants, data and inits
modConstants <- list(Nreads = Nreads, Npool = Npool, Nspp = Nspp, x = amount_DNA)
modData <- list(y = number_Reads)
modInits <- list(a = rep(1, Nspp), alpha = matrix(1, nrow = Npool, ncol = Nspp), 
                 p = matrix(rep(1/Nspp, Nspp), nrow = Npool, ncol = Nspp, byrow = TRUE))

## build model
mod <- nimbleModel(code = modCode, name = 'mod', constants = modConstants, 
                   data = modData, inits = modInits)

Cmod <- compileNimble(mod)

modConf <- configureMCMC(mod)
modConf$addMonitors('a')
modConf$setThin(20)
modConf$addSampler(target = sprintf('a[1:%s]', modConstants$Nspp), type = 'RW_block')

modMCMC <- buildMCMC(modConf)
CmodMCMC <- compileNimble(modMCMC, project = mod)

N <- 1000
burn <- 300
niter <- (N + burn) * modConf$thin
CmodMCMC$run(niter)

samp <- as.matrix(CmodMCMC$mvSamples)

plot(samp[-(1:burn), 1], type = 'l')
abline(h = a[1])
