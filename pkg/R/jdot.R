##
## scaling constant for points in SOC cones
setMethod("jdot", signature = c("numeric", "numeric"), function(u, v){
    u[1] * v[1] - drop(crossprod(u[-1], v[-1]))
})
setMethod("jdot", signature = c("SOCV", "SOCV"), function(u, v){
    u@u[1, 1] * v@u[1, 1] - drop(crossprod(u@u[-1, 1], v@u[-1, 1]))
})
setMethod("jdot", signature = c("SOCV", "missing"), function(u, v){
    u@u[1, 1]^2 - drop(crossprod(u@u[-1, 1]))
})
