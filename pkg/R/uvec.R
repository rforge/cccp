##
## uvec-method for objects of class union LNLV
setMethod("uvec", signature = "LNLV", function(u){
    u
})
##
## uvec-method for SOCV
setMethod("uvec", signature = "SOCV", function(u){
    u
})
##
## uvec-method for PSDV
setMethod("uvec", signature = "PSDV", function(u){
    dim(u@u) <- c(u@dims^2, 1)
    u
})
