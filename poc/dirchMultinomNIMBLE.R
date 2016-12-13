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

## nimble functions for the dirichlet-multinomial distribution

ddirchmulti <- nimbleFunction(
    run = function(x = double(1), alpha = double(1), size = double(0),
                   log = integer(0, default = 0)) {
        returnType(double(0))
        alpha0 <- sum(alpha)
        
        ## new log prob that ignores 0's instead of throwing NaN/Inf
        lgammaSum <- numeric(length = length(x), value = 0, init = TRUE)
        for(i in 1:length(lgammaSum)) {
            if(x[i] > 0) {
                lgammaSum[i] <- log(x[i]) + lgamma(alpha[i]) +
                                        lgamma(x[i]) -
                                        lgamma(alpha[i] + x[i])
            }
        }
        logProb <- log(size) +
            lgamma(alpha0) + lgamma(size) - lgamma(alpha0 + size) -
            sum(lgammaSum)
        
        if(log) return(logProb)
        else return(exp(logProb))
    }
)

rdirchmulti <- nimbleFunction(
    run = function(n = double(0, default = 1), alpha = double(1), size = double(0)) {
        returnType(double(1))
        
        ## modified from MCMCpack to allow alpha_k = 0
        x <- numeric(length = length(alpha), value = 0, init = TRUE)
        for(i in 1:length(x)) x[i] <- rgamma(1, shape = alpha[i], rate = 1)
        p <- x/sum(x)
        
        return(rmulti(1, size = size, prob = p))
    }
)

## the model as NIMBLE code
modCode <- nimbleCode({
    for(i in 1:Npool) {
        alpha[i, 1:Nspp] <- x[i, 1:Nspp] * a[1:Nspp]
        y[i, 1:Nspp] ~ ddirchmulti(alpha[i, 1:Nspp], Nreads[i])
    }
    
    ## priors
    for(j in 1:Nspp) {
        a[j] ~ dexp(0.00001)
    }

})

## model constants, data and inits
modConstants <- list(Nreads = Nreads, Npool = Npool, Nspp = Nspp, x = amount_DNA)
modData <- list(y = number_Reads)
modInits <- list(a = rep(1, Nspp), alpha = matrix(1, nrow = Npool, ncol = Nspp))

## build model
mod <- nimbleModel(code = modCode, name = 'mod', constants = modConstants, 
                   data = modData, inits = modInits)

Cmod <- compileNimble(mod)

modConf <- configureMCMC(mod)
modConf$addMonitors('a')
modConf$setThin(50)
modConf$addSampler(target = sprintf('a[1:%s]', modConstants$Nspp), type = 'RW_block')

modMCMC <- buildMCMC(modConf)
CmodMCMC <- compileNimble(modMCMC, project = mod)

N <- 10000
burn <- 300
niter <- (N + burn) * modConf$thin
CmodMCMC$run(niter)

samp <- as.matrix(CmodMCMC$mvSamples)

i <- 1
plot(samp[-(1:burn), i], type = 'l')
abline(h = mean(samp[-(1:burn), i]), col = 'blue')
abline(h = a[i], col = 'red')
i <- i + 1