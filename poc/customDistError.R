library(nimble)

sessionInfo()
# R version 3.3.0 (2016-05-03)
# Platform: x86_64-apple-darwin13.4.0 (64-bit)
# Running under: OS X 10.11.4 (El Capitan)
# 
# locale:
# [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
# 
# attached base packages:
# [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
# [1] nimble_0.5-1
# 
# loaded via a namespace (and not attached):
# [1] magrittr_1.5     igraph_1.0.1     coda_0.18-1      codetools_0.2-14 grid_3.3.0       
# [6] lattice_0.20-33 


dmyexp <- nimbleFunction(
    run = function(x = double(0), rate = double(0, default = 1),
                   log = integer(0, default = 0)) {
        returnType(double(0))
        logProb <- log(rate) - x*rate
        if(log) return(logProb)
        else return(exp(logProb))
    })
    
rmyexp <- nimbleFunction(
    run = function(n = integer(0), rate = double(0, default = 1)) {
        returnType(double(0))
        if(n != 1) print("rmyexp only allows n = 1; using n = 1.")
        dev <- runif(1, 0, 1)
        return(-log(1-dev) / rate)
    })
    
pmyexp <- nimbleFunction(run = function(q = double(0), rate = double(0, default = 1),
                                        lower.tail = integer(0, default = 1),
                                        log.p = integer(0, default = 0)) {
    returnType(double(0))
    if(!lower.tail) {
        logp <- -rate * q
        if(log.p) return(logp)
        else return(exp(logp))
    } else {
        p <- 1 - exp(-rate * q)
        if(!log.p) return(p)
        else return(log(p))
    }
})

qmyexp <- nimbleFunction(
    run = function(p = double(0), rate = double(0, default = 1),
                   lower.tail = integer(0, default = 1),
                   log.p = integer(0, default = 0)) {
        returnType(double(0))
        if(log.p) p <- exp(p)
        if(!lower.tail) p <- 1 - p
        return(-log(1 - p) / rate)
    })
    
ddirchmulti <- nimbleFunction(
    run = function(x = double(1), alpha = double(1), size = double(0),
                   log = integer(0, default = 0)) {
        returnType(double(0))
        logProb <- lgamma(size) - sum(lgamma(x)) + lgamma(sum(alpha)) -
            sum(lgamma(alpha)) + sum(lgamma(alpha + x)) - lgamma(sum(alpha) +
                                                                     size)
        if(log) return(logProb)
        else return(exp(logProb))
    })
    
rdirchmulti <- nimbleFunction(
    run = function(n = integer(0), alpha = double(1), size = double(0)) {
        returnType(double(1))
        if(n != 1) print("rdirchmulti only allows n = 1; using n = 1.")
        p <- rdirch(1, alpha)
        return(rmulti(1, size = size, prob = p))
    })

code <- nimbleCode({
    y[1:K] ~ ddirchmulti(alpha[1:K], n)
    for(i in 1:K) {
        alpha[i] ~ dmyexp(1/3)
    }
})

# model <- nimbleModel(code, constants = list(K = 5))
# defining model...
# Error in BUGSdeclClassObject$setup(code[[i]], contextID, lineNumber) : 
#   Improper syntax for stochastic declaration: y[1:K] ~ ddirchmulti(alpha[1:K], n)


## trying to see if registering custom functions fixes anything
registerDistributions(list(
    ddirchmulti = list(
        BUGSdist = 'ddirchmulti(alpha, size)',
        discrete = TRUE,
        types = c('value = integer(1)', 'alpha = double(1)', 'size = integer(0)')
        ), 
    dmyexp = list(
        BUGSdist = "dmyexp(rate, scale)",
        BUGSdist = "dmyexp(rate = 1/scale)",
        Rdist = "dmyexp(rate = 1/scale)",
        altParams = c("scale = 1/rate", "mean = 1/rate"),
        pqAvail = TRUE,
        range = c(0, Inf)
        )
    ))
    
model <- nimbleModel(code, constants = list(K = 5, n = 25), data = list(y = rep(5, 5)),
                    inits = list(alpha = rep(1, 5)))
model$calculate('y')
model$resetData()
model$simulate('y')
model$y
model$calculate('y')

# defining model...
# building model...
# checking model...   (use nimbleModel(..., check = FALSE) to skip model check)
# Error in model$check() : 
#   Dimension of distribution argument(s) 'value,alpha' does not match required dimension(s) 
#   for the distribution 'ddirchmulti'. Necessary dimension(s) are 0,0.
