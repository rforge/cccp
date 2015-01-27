##
## Value of variable 'x'
setMethod("getx", signature = "PDV", function(object){
    object@x
})
setMethod("getx", signature = "CPS", function(object){
    object@x
})
