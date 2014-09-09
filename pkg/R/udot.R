##
## inner-product between two vectors
setMethod("udot", signature = c("numeric", "numeric"), function(u, v){
    drop(crossprod(u, v))
})
setMethod("udot", signature = c("numeric", "missing"), function(u, v){
    drop(crossprod(u))
})
##
## inner-product between points in NNO cones
setMethod("udot", signature = c("NNOV", "NNOV"), function(u, v){
    drop(crossprod(u@u, v@u))
})
setMethod("udot", signature = c("NNOV", "missing"), function(u, v){
    drop(crossprod(u@u))
})
##
## inner-product between points in SOC cones
setMethod("udot", signature = c("SOCV", "SOCV"), function(u, v){
    drop(crossprod(u@u, v@u))
})
setMethod("udot", signature = c("SOCV", "missing"), function(u, v){
    drop(crossprod(u@u))
})
##
## inner-product between points in PSD cones
setMethod("udot", signature = c("PSDV", "PSDV"), function(u, v){
    drop(sum(u@u * v@u))
})
setMethod("udot", signature = c("PSDV", "missing"), function(u, v){
    drop(sum(u@u^2))
})
