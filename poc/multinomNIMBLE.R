library(nimble)

setwd('~/Dropbox/hawaiiDimensions/metabarcoding/poc')

## read data for calibrating simulation
diffMarkers <- read.csv('clean_diffMarkers.csv', as.is=TRUE)

## simulate data
Npool <- 1000
Nreads <- Npool * 1000
a <- 0.5
amount_DNA <- runif(Npool, 1, 100)
number_Reads <- as.vector(rmultinom(1, Nreads, rdirch(1, a*amount_DNA)))


## the model as NIMBLE code
modCode <- nimbleCode({
    alpha[1:Npool] <- x[1:Npool] * a
    p[1:Npool] ~ ddirch(alpha[1:Npool])
    y[1:Npool] ~ dmulti(p[1:Npool], Nreads)
    
    ## priors
    a ~ dexp(0.00001)
})


## model constants, data and inits
modConstants <- list(Nreads = Nreads, Npool = Npool, x = amount_DNA)
modData <- list(y = number_Reads)
modInits <- list(a = 1, alpha = rep(1, Npool), p = rep(1/Npool, Npool))

## build model
mod <- nimbleModel(code = modCode, name = 'mod', constants = modConstants, 
                   data = modData, inits = modInits)

Cmod <- compileNimble(mod)

modConf <- configureMCMC(mod)
modConf$addMonitors('a')
modConf$setThin(20)

modMCMC <- buildMCMC(modConf)
CmodMCMC <- compileNimble(modMCMC, project = mod)

N <- 1000
burn <- 300
niter <- (N + burn) * modConf$thin
CmodMCMC$run(niter)

samp <- as.matrix(CmodMCMC$mvSamples)

par(mfrow = c(1, 2), mar = c(3, 2, 0, 1) + 0.1)
plot(samp[-(1:burn), 'a'], type = 'l')
abline(h = a, col = 'red')
acf(samp[-(1:burn), 'a'])

