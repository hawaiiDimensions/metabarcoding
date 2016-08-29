library(nimble)

setwd('~/Dropbox/hawaiiDimensions/metabarcoding/poc')

## simulate data
N <- 1000
alpha <- 0.5
amount_DNA <- runif(N, 1, 100)
number_Reads <- rmultinom(1, N * 1000, rdirch(1, alpha*amount_DNA))


## the model as NIMBLE code
modCode <- nimbleCode({
    ## loop over all data
    for(i in 1:N) {
        ## u is the random intercept
        mu[i] <- b0 + b1*x1[i] + u[g[i]]
        y[i] ~ dnorm(mu[i], sd = 1/sqrt(tau))
    }
    
    ## random group intercepts
    for(j in 1:Ng) {
        u[j] ~ dnorm(0, sd = 1/sqrt(tau_g))
    }
    
    ## priors
    b0 ~ dnorm(0, sd = 1000)
    b1 ~ dnorm(0, sd = 1000)
    tau ~ dgamma(0.001, 0.001)
    tau_g ~ dgamma(0.001, 0.001)
})


## model constants, data and inits
modConstants <- list(N = n, Ng = ng, x1 = x1, g = g)
modData <- list(y = y)
modInits <- list(b0 = 0, b1 = 1, tau = 1, tau_g = 1, u = rep(0, ng))

## build model
mod <- nimbleModel(code = modCode, name = 'mod', constants = modConstants, 
                   data = modData, inits = modInits)

Cmod <- compileNimble(mod)

modConf <- configureMCMC(mod)
modConf$addMonitors(c('b0', 'b1', 'tau', 'tau_g'))
modConf$addSampler(target = c('tau', 'tau_g'), type = 'RW_block')
modConf$setThin(50)

modMCMC <- buildMCMC(modConf)
CmodMCMC <- compileNimble(modMCMC, project = mod)

N <- 1000
burn <- 50
niter <- (N + burn) * modConf$thin
CmodMCMC$run(niter)

samp <- as.matrix(CmodMCMC$mvSamples)


plot(samp[-(1:50), 'b0'], type = 'l')
abline(h = b0, col = 'red')
acf(samp[-(1:50), 'b0'])

plot(samp[-(1:50), 'b1'], type = 'l')
abline(h = b1, col = 'red')

plot(samp[-(1:50), 'tau_g'], type = 'l')
abline(h = 1/sig_g^2, col = 'red')
acf(samp[-(1:50), 'tau_g'])




diffMarkers <- read.csv('clean_diffMarkers.csv', as.is=TRUE)
pcrCycle <- read.csv('clean_pcrCycle.csv', as.is=TRUE)

## ================================
## define constants, data and inits
## ================================

## let x be the dataset to use
x <- diffMarkers
if(length(unique(x$PCR_cycles)) == 1) {
    treat <- 'Primer'
} else {
    treat <- 'PCR_cycles'
}

bla <- tapply(x$amount_DNA, list(as.factor(x$Pool), as.factor(x$Specimen), as.factor(x[, treat])), function(a) {
    max(a, na.rm=TRUE)
})


## aggregate on amount by different things then index by those things

## constants and explanitory variables
metabConsts <- list(npool = length(unique(x$Pool)),
                    nread = tapply(x$total_Reads, x$Pool, max),
                    ## rows of amountDNA correspond to pools, columns to specimens
                    amountDNA = tapply(x$amount_DNA, list(x$Pool, x$Specimen), function(a) {
                        if(all(is.na(a))) return(0)
                        else return(max(a, na.rm=TRUE))
                    }),
                    treatment = sort(unique(x[, treat])))





a <- 4
b <- 1
set.seed(0)
pumpData <- list(x = rpois(ntime, pumpConsts$t * rgamma(ntime, a, b)))
pumpInits <- list(alpha = 1, beta = 1,
                  theta = rep(0.1, pumpConsts$N))



code <- nimbleCode({
    y[1:K] ~ dmulti(p[1:K], n)
    p[1:K] ~ ddirch(alpha[1:K])
    log(alpha[1:K]) ~ dmnorm(alpha0[1:K], R[1:K, 1:K])
})

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