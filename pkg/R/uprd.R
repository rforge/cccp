##
## product of points in LNL cones
setMethod("uprd", signature = c("LNLV", "LNLV"), function(u, v){
    u@u <- .sprd_nl(u@u, v@u)
    u
})
setMethod("uprd", signature = c("LNLV", "missing"), function(u, v){
    u@u <- .sprd_nl(u@u, u@u)
    u
})
##
## product of points in SOC cones
setMethod("uprd", signature = c("SOCV", "SOCV"), function(u, v){
    u@u <- .sprd_s(u@u, v@u)
    u
})
setMethod("uprd", signature = c("SOCV", "missing"), function(u, v){
    u@u <- .sprd_s(u@u, u@u)
    u
})
##
## product of points in PSD cones
setMethod("uprd", signature = c("PSDV", "PSDV"), function(u, v){
    u@u <- .sprd_p(u@u, v@u, u@dims)
    u
})
setMethod("uprd", signature = c("PSDV", "missing"), function(u, v){
    u@u <- .sprd_p(u@u, u@u, u@dims)
    u
})
