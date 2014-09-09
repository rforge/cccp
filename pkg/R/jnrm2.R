##
## Norm of a scaling point in SOC cones
setMethod("jnrm2", signature = "SOCV", function(u){
    a <- drop(sqrt(crossprod(u@u[-1])))
    sqrt(u@u[1, 1] - a) * sqrt(u@u[1, 1] + a)
})

