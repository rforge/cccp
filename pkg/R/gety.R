##
## Value of variable 'y'
setMethod("gety", signature = "PDV", function(object){
    object@y
})
setMethod("gety", signature = "CPS", function(object){
    object@y
})
