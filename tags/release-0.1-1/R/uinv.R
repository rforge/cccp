##
## inverse of product of points in LNL cones
setMethod("uinv", signature = c("LNLV", "LNLV"), function(u, v){
    u@u <- .sinv_nl(u@u, v@u)
    u
})
##
## inverse of product of points in SOC cones
setMethod("uinv", signature = c("SOCV", "SOCV"), function(u, v){
    u@u <- .sinv_s(u@u, v@u)
    u
})
##
## inverse of product of points in PSD cones
setMethod("uinv", signature = c("PSDV", "PSDV"), function(u, v){
    u@u <- matrix(.sinv_p(u@u, v@u, u@dims))
    u
})
