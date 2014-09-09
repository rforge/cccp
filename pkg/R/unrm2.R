##
## Norm of a vector
setMethod("unrm2", signature = "numeric", function(u){
    sqrt(udot(u))
})
##
## Norm of a point in NNO cones
setMethod("unrm2", signature = "NNOV", function(u){
    sqrt(udot(u))
})
##
## Norm of a point in SOC cones
setMethod("unrm2", signature = "SOCV", function(u){
    sqrt(udot(u))
})
##
## Norm of a point in SOC cones
setMethod("unrm2", signature = "PSDV", function(u){
    sqrt(udot(u))
})

