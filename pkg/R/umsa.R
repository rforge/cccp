##
## Adjusting LNL-variable by maximum step-size
setMethod("umsa", signature = "LNLV", function(u, alpha, init = TRUE){
    u@u <- .smsa_nl(u@u, alpha, init)
    u
})
##
## Adjusting SOC-variable by maximum step-size
setMethod("umsa", signature = "SOCV", function(u, alpha, init = TRUE){
    u@u <- .smsa_s(u@u, alpha, init)
    u
})
##
## Adjusting PSD-variable by maximum step-size
setMethod("umsa", signature = "PSDV", function(u, alpha, init = TRUE, sigma = NULL, lambda = NULL){
    if(init){
        u@u <- .smsa1_p(u@u, alpha, u@dims)
    } else {
        u@u <- .smsa2_p(u@u, alpha, sigma, lambda@u, u@dims)
    }
    u
})
