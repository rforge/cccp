##
## Maximum step-size for NNO-variables
setMethod("umss", signature = "LNLV", function(u){
    .smss_nl(u@u)
})
##
## Maximum step-size for SOC-variables
setMethod("umss", signature = "SOCV", function(u){
    .smss_s(u@u)
})
##
## Maximum step-size for PSD-variables
setMethod("umss", signature = "PSDV", function(u){
    .smss_p(u@u, u@dims)
})
