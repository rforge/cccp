##
## umat-method for objects of class NNOV
setMethod(f = "umat", "NNOV", definition = function(u){
    u
})
##
## umat-method for objects of class SOCV
setMethod(f = "umat", "SOCV", definition = function(u){
    u
})
##
## umat-method for objects of class PSDV
setMethod(f = "umat", "PSDV", definition = function(u){
    dim(u@u) <- c(u@dims, u@dims)
    u
})
