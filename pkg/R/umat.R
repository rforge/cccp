##
## umat-method for objects of class union LNLV
setMethod(f = "umat", "LNLV", definition = function(u){
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
