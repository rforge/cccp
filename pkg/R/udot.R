##
## inner-product between two vectors
setMethod("udot", signature = c("numeric", "numeric"), function(u, v){
    .sdot_nls(u, v)
})
setMethod("udot", signature = c("numeric", "missing"), function(u, v){
    .sdot_nls(u, u)
})
##
## inner-product between points in LNL cones
setMethod("udot", signature = c("LNLV", "LNLV"), function(u, v){
    .sdot_nls(u@u, v@u)
})
setMethod("udot", signature = c("LNLV", "missing"), function(u, v){
    .sdot_nls(u@u, u@u)
})
##
## inner-product between points in SOC cones
setMethod("udot", signature = c("SOCV", "SOCV"), function(u, v){
    .sdot_nls(u@u, v@u)
})
setMethod("udot", signature = c("SOCV", "missing"), function(u, v){
    .sdot_nls(u@u, u@u)
})
##
## inner-product between points in PSD cones
setMethod("udot", signature = c("PSDV", "PSDV"), function(u, v){
    u <- matrix(u@u, u@dims, u@dims)
    v <- matrix(v@u, v@dims, v@dims)
    a1 <- crossprod(diag(u), diag(v))
    a2 <- 2 * crossprod(u[lower.tri(u)], v[lower.tri(v)])
    drop(a1 + a2)
})
setMethod("udot", signature = c("PSDV", "missing"), function(u, v){
    u <- matrix(u@u, u@dims, u@dims)
    a1 <- crossprod(diag(u))
    a2 <- 2 * crossprod(u[lower.tri(u)])
    drop(a1 + a2)
})
