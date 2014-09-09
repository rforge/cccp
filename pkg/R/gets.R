##
## Value of variable 's'
setMethod("gets", signature = "PDV", function(object){
    lapply(object@s, function(s) umat(s)@u)
})
setMethod("gets", signature = "CPS", function(object){
    lapply(object@s, function(s) umat(s)@u)
})
