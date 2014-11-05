##
## scaling constant for points in SOC cones
setMethod("jdot", signature = c("numeric", "numeric"), function(u, v){
    .jdot(matrix(u), matrix(v))
})
setMethod("jdot", signature = c("SOCV", "SOCV"), function(u, v){
    .jdot(u@u, v@u)
})
setMethod("jdot", signature = c("SOCV", "missing"), function(u, v){
    .jdot(u@u, u@u)
})
