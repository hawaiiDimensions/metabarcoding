## =============================================================
## functions to run nimble model and calculate model predictions
## and posterior R^2
## =============================================================

library(nimble)
library(coda)
library(plyr)
library(reshape2)
library(parallel)

## function to run nimble model
## Nreads: number of total reads per pool
## amount_DNA: amount of DNA input per spp per pool (pool = rows; spp = columns)
## number_Reads: number of reads per spp per pool (pool = rows; spp = columns)

runNimble <- function(Nreads, amount_DNA, number_Reads, N = 10000, thin = 50, burn = 50) {
    ## trim data down to only those species captured at least once
    amount_DNA <- amount_DNA[, colSums(number_Reads) > 0]
    number_Reads <- number_Reads[, colSums(number_Reads) > 0]
    
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
    
    ## make block sampler have bigger step for variables with big estimated slopes
    
    m <- sapply(1:ncol(modConstants$x), function(i) {
        lm(c(0, modData$y[, i]) ~ c(0, modConstants$x[, i]))$coefficients[2]
    })
    names(m) <- NULL
    m[m <= 0] <- min(m[m > 0])
    m <- 1 + 5 * (log(m) - log(min(m)))
    
    modConf$addSampler(target = sprintf('a[1:%s]', modConstants$Nspp),
                       type = 'RW_block', 
                       control = list(propCov = diag(m, modConstants$Nspp)))
    
    ## build mcmc
    
    modMCMC <- buildMCMC(modConf)
    CmodMCMC <- compileNimble(modMCMC, project = mod)
    
    niter <- (N + burn) * modConf$thin
    CmodMCMC$run(niter)
    
    samp <- as.matrix(CmodMCMC$mvSamples)[-(1:burn), ]
    
    return(samp)
}


## function to take data.frame, parse it by explanitory variables and run
## nimble model
## dat: the data frame
## x: the explanitory variable that differs across experiments (e.g. marker, num cycle, etc)

buildNRunMod <- function(dat, x, N = 10000, thin = 50, burn = 50) {
    ## number or reads for each primer and pool combination (rows:x, cols:pools)
    totReads <- acast(melt(dat, c(x, 'Pool'), 'total_Reads'), 
                      as.formula(paste(x, 'Pool', sep = ' ~ ')), 
                      value.var =  'value', max, na.rm = TRUE, fill = 0)
    
    ## amount of DNA from each species in each pool (rows:pools, cols:species)
    amountDNA <- acast(melt(dat, c('Pool', 'Specimen'), 'amount_DNA'), Pool ~ Specimen, 
                       value.var =  'value', max, na.rm = TRUE, fill = 0)
    
    ## number of reads per primer, pool, species combo (dim1:x, dim2:pool, dim3:species)
    numReads <- acast(melt(dat, c(x, 'Pool', 'Specimen'), 'number_Reads'), 
                      as.formula(paste(x, 'Pool', 'Specimen', sep = ' ~ ')),
                      value.var =  'value', max, na.rm = TRUE, fill = 0)
    
    ## loop over explanitory var, fitting model to each and calculating: 
    ## R2 
    ## effective sample size
    ## Geweke's convergence test
    
    # out <- mclapply(1:nrow(totReads), mc.cores = 4, FUN = function(i) {
    out <- lapply(1:nrow(totReads), function(i) {
        ## model parameters
        modPar <- try(runNimble(totReads[i, ], amountDNA, numReads[i, , ], 
                            N = N, thin = thin, burn = burn))
        if(class(modPar) == 'try-error') browser()
        
        ## return R2 and (across all a's) min effective size and Geweke's test
        out <- list(par = modPar, 
                    summ = c(R2 = bayesR2(numReads[i, , ], totReads[i, ], modPar), 
                             minESS = min(effectiveSize(modPar)), 
                             nGewekeFail = sum(abs(geweke.diag(modPar)$z) > 1.96)))
        
        if(class(out) == 'try-error') browser()
        return(out)
    })
    
    
    ## extract posterior parameter samples and summary from output
    outPar <- lapply(out, function(x) x$par)
    names(outPar) <- rownames(totReads)
    
    outSumm <- data.frame(rownames(totReads), t(sapply(out, function(x) x$summ)))
    names(outSumm)[1] <- x
    
    return(list(par = outPar, summ = outSumm))
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
    
    out <- c(1 - mean(varR) / varY, quantile(1 - varR / varY, c(0.025, 0.975)))
    names(out) <- c('mean', 'ciLo', 'ciHi')
    
    return(out)
}
