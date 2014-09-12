##
## Maximum step-size for NNO-variables
setMethod("umss", signature = "LNLV", function(u){
    list(ms = -min(u@u), evd = NULL)
})
##
## Maximum step-size for SOC-variables
setMethod("umss", signature = "SOCV", function(u){
    list(ms = drop(sqrt(crossprod(u@u[-1]))) - u@u[1], evd = NULL)
})
##
## Maximum step-size for PSD-variables
setMethod("umss", signature = "PSDV", function(u){
    dim(u@u) <- c(u@dims, u@dims)
    evd <- eigen(u@u)
    ev <- evd$values
    list(ms = -ev[u@dims], evd = evd)
})
