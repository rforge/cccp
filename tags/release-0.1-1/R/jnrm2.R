##
## Norm of a scaling point in SOC cones
setMethod("jnrm2", signature = "SOCV", function(u){
    .jnrm2(u@u)
})

