##
## Initial point for NNO-slack variables
setMethod("initp", signature = "NNOC", function(object){
    new("NNOV", u = Matrix(1, nrow = object@dims, ncol = 1), dims = object@dims)
})
##
## Initial point for NLF-slack variables
setMethod("initp", signature = "NLFC", function(object){
    new("NLFV", u = Matrix(1, nrow = object@dims, ncol = 1), dims = object@dims)
})
##
## Initial point for SOC-slack variables
setMethod("initp", signature = "SOCC", function(object){
    p <- Matrix(0, nrow = object@dims, ncol = 1)
    p[1, 1] <- 1
    new("SOCV", u = p, dims = object@dims)
})
##
## Initial point for PSD-slack variables
setMethod("initp", signature = "PSDC", function(object){
    p <- Diagonal(object@dims)
    dim(p) <- c(object@dims^2, 1)
    new("SOCV", u = p, dims = object@dims)
})

