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
    .sdot_p(u@u, v@u, u@dims)
})
setMethod("udot", signature = c("PSDV", "missing"), function(u, v){
    .sdot_p(u@u, u@u, u@dims)
})
