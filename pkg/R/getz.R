##
## Value of variable 'z'
setMethod("getz", signature = "PDV", function(object){
    lapply(object@z, function(z) umat(z)@u)
})
setMethod("getz", signature = "CPS", function(object){
    lapply(object@z, function(z) umat(z)@u)
})
