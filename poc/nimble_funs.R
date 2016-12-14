## =============================================================
## functions to run nimble model and calculate model predictions
## and posterior R^2
## =============================================================

library(nimble)
library(coda)
library(plyr)
library(reshape2)
library(parallel)


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


## function to run nimble model
## Nreads: number of total reads per pool
## amount_DNA: amount of DNA input per spp per pool (pool = rows; spp = columns)
## number_Reads: number of reads per spp per pool (pool = rows; spp = columns)

runNimble <- function(Nreads, amount_DNA, number_Reads, N = 10000, thin = 50, burn = 50) {
    ## number of pools
    Npool <- length(Nreads)
    
    ## number of spp
    Nspp <- ncol(amount_DNA)
    
    ## the model as NIMBLE code
    modCode <- nimbleCode({
        for(i in 1:Npool) {
            alpha[i, 1:Nspp] <- x[i, 1:Nspp] * a[1:Nspp]
            # p[i, 1:Nspp] ~ ddirch(alpha[i, 1:Nspp])
            # y[i, 1:Nspp] ~ dmulti(p[i, 1:Nspp], Nreads[i])
            y[i, 1:Nspp] ~ ddirchmulti(alpha[i, 1:Nspp], Nreads[i])
        }
        
        ## priors
        for(j in 1:Nspp) {
            a[j] ~ dexp(0.00001)
        }
    })
    
    ## dirichlet can't have alpha_i = 0, so add a little to any case where there are 0's
    # if(any(amount_DNA == 0)) amount_DNA <- amount_DNA + .Machine$double.eps
    
    ## model constants, data and inits
    modConstants <- list(Nreads = Nreads, Npool = Npool, Nspp = Nspp, x = amount_DNA)
    modData <- list(y = number_Reads)
    modInits <- list(a = rep(1, Nspp), alpha = matrix(1, nrow = Npool, ncol = Nspp)#, 
                     # p = matrix(rep(1/Nspp, Nspp), nrow = Npool, ncol = Nspp, byrow = TRUE)
                     )
    
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
                       control = list(scale = 4))
    
    ## build mcmc
    
    modMCMC <- buildMCMC(modConf)
    CmodMCMC <- compileNimble(modMCMC, project = mod)
    
    niter <- (N + burn) * modConf$thin
    CmodMCMC$run(niter)
    
    samp <- as.matrix(CmodMCMC$mvSamples)[-(1:burn), ]
    samp
    
    return(samp)
}


## function to take data.frame, parse it by explanitory variables and run
## nimble model
## dat: the data frame
## x: the explanitory variable that differs across experiments (e.g. marker, num cycle, etc)

buildNRunMod <- function(dat, x, N = 10000, thin = 50, burn = 50) {
    ## format data
    fdat <- .formatData(dat, x)
    
    # loop over explanitory var, fitting model to each and calculating:
    # R2
    # effective sample size
    # Geweke's convergence test
    # out <- mclapply(1:nrow(totReads), mc.cores = 4, FUN = function(i) {
    out <- lapply(1:nrow(fdat[[1]]), function(i) {
        .internalRunNimble(fdat, i, N, thin, burn)
    })
    
    ## extract posterior parameter samples and summary from output
    outPar <- lapply(out, function(x) x$par)
    names(outPar) <- rownames(totReads)
    
    outSumm <- data.frame(rownames(totReads), t(sapply(out, function(x) x$summ)))
    names(outSumm)[1] <- x
    
    return(list(par = outPar, summ = outSumm))
}

## helper function to get data into right format
.formatData <- function(dat, x) {
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
    
    return(list(totReads, amountDNA, numReads))
}


## helper function to prep already formatted data (fdat) for nimble model
.prepData <- function(fdat, i) {
    ## extract needed data
    thisNreads <- fdat[[1]][i, ]
    thisAmount_DNA <- fdat[[2]]
    thisNumber_Reads <- fdat[[3]][i, , ]
    
    ## trim data down to only those pools with any reads at all
    thisAmount_DNA <- thisAmount_DNA[thisNreads > 0, ]
    thisNumber_Reads <- thisNumber_Reads[thisNreads > 0, ]
    thisNreads <- thisNreads[thisNreads > 0]
    
    ## trim data down to only those species captured at least once
    thisAmount_DNA <- thisAmount_DNA[, colSums(thisNumber_Reads) > 0]
    thisNumber_Reads <- thisNumber_Reads[, colSums(thisNumber_Reads) > 0]
    
    return(list(thisNreads, thisAmount_DNA, thisNumber_Reads))
}

## helper function to do all the nitty gritty of running and summarizing 
## the nimble model
.internalRunNimble <- function(fdat, i, N, thin, burn) {
    ## extract needed data
    pdat <- .prepData(fdat, i)
    thisNreads <- pdat[[1]]
    thisAmount_DNA <- pdat[[2]]
    thisNumber_Reads <- pdat[[3]]
    
    ## model parameters
    modPar <- try(runNimble(thisNreads, thisAmount_DNA, thisNumber_Reads,
                            N = N, thin = thin, burn = burn))
    
    ## error catching
    if(class(modPar) == 'try-error') {
        return(list(par = NA, summ = NA))
    } else {
        ## return R2 and (across all a's) min effective size and Geweke's test
        out <- try(list(par = modPar,
                        summ = c(R2 = bayesR2(thisNumber_Reads, thisNreads, 
                                              thisAmount_DNA, modPar),
                                 minESS = min(effectiveSize(modPar)),
                                 nGewekeFail = sum(abs(geweke.diag(modPar)$z) > 1.96))))
        if(class(out) == 'try-error') {
            return(list(par = NA, summ = NA))
        } else {
            return(out)
        }
    }
}


## function to calculate predicted values for multinomial-dirichlet model
## n is vector of number of trials (e.g. total number of reads)
## x is explanitory var (matrix, nrow = length(n))
## a is a vector of parameter estimates for dirichlet distrib (length = ncol(x))

predictMultiDir <- function(n, x, a) {
    alpha <- t(t(x) * a)
    alpha <- alpha / rowSums(alpha)
    
    as.vector(alpha) * rep(n, ncol(x))
}


## function to calculate Bayesian R^2 (not general, only for multi-dir)
## y is matrix of data (rows = diff values of n)
## n is vector of number of trials
## x is matrix of explanitory variable (same dim as y)
## A is matrix of posterior parameter estimates for dirichlet dirstrib
## (one row per posterior sample)

bayesR2 <- function(y, n, x, A) {
    ## total variance of Y
    varY <- var(as.vector(y))
    
    varR <- sapply(1:nrow(A), function(i) {
        yhat <- predictMultiDir(n, x, A[i, ])
        var(as.vector(y - yhat))
    })
    
    out <- c(1 - mean(varR) / varY, quantile(1 - varR / varY, c(0.025, 0.975)))
    names(out) <- c('mean', 'ciLo', 'ciHi')
    
    return(out)
}


## function to calculate model-predicted amount of DNA given fitted parameters and data
## y is matrix of data (number of reads; rows = diff values of n)
## n is vector of number of trials (total reads)
## x is matrix of explanitory variable (amount DNA; same dim as y)
## A is matrix of posterior parameter estimates for dirichlet dirstrib
## (one row per posterior sample)

predictDNA <- function(y, n, x, A) {
    out <- sapply(1:nrow(A), function(j) {
        a <- A[j, ]
        
        xcor <- lapply(1:nrow(y), function(i) {
            ## eigen vector solution to get predicted x
            B <- matrix(rep(y[i, ]/(n[i] * a), length(a)), ncol = length(a)) %*% diag(a)
            v1 <- Re(eigen(B)$vectors[, 1])
            
            return(v1 * sum(x[i, ]) / sum(v1))
        })
        xcor <- unlist(xcor)
        
        return(xcor)
    })
    
    out <- aaply(out, 1, function(X) c(mean(X), quantile(X, prob = c(0.05, 0.975))))
    dimnames(out) <- NULL
    out <- cbind(as.vector(t(x)), out)
    colnames(out) <- c('x', 'pred_mean', 'pred_ciLo', 'pred_ciHi')
    
    return(out)
}
