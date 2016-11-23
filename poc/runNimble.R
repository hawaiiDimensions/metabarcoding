## =============================================================
## functions to run nimble model and calculate model predictions
## and posterior R^2
## =============================================================

## function to run nimble model
## Nreads: number of total reads per pool
## amount_DNA: amount of DNA input per spp per pool (pool = rows; spp = columns)
## number_Reads: number of reads per spp per pool (pool = rows; spp = columns)

runNimble <- function(Nreads, amount_DNA, number_Reads, thin = 20, N = 1000, burn = 300) {
    ## number of pools
    Npool <- length(Nreads)
    
    ## number of spp
    Nspp <- ncol(amount_DNA)
    
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
    modConf$setThin(thin)
    modConf$addSampler(target = sprintf('a[1:%s]', modConstants$Nspp), 
                       type = 'RW_block', control = list(scale = 4))
    
    modMCMC <- buildMCMC(modConf)
    CmodMCMC <- compileNimble(modMCMC, project = mod)
    
    niter <- (N + burn) * modConf$thin
    CmodMCMC$run(niter)
    
    samp <- as.matrix(CmodMCMC$mvSamples)[-(1:burn), ]
    
    return(samp)
}

## function to calculate predicted values for multinomial-dirichlet model
## n is vector of number of trials (e.g. total number of reads)
## a is a vector of parameter estimates for dirichlet distrib

predictMultiDir <- function(n, a) {
    outer(n, a/sum(a))
}


## function to calculate Bayesian R^2 (not general, only for multi-dir)
## y is matrix of data (rows = diff values of n)
## n is vector of number of rials
## A is matrix of posterior parameter estimates for dirichlet dirstrib
## (one row per posterior sample)

bayesR2 <- function(y, n, A) {
    ## total variance of Y
    varY <- var(as.vector(y))
    
    varR <- sapply(1:nrow(A), function(i) {
        yhat <- predictMultiDir(n, A[i, ])
        var(as.vector(y - yhat))
    })
    
    return(1 - mean(varR) / varY)
}
