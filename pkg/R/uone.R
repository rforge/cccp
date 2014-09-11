##
## One-element for points in LNL cones
setMethod("uone", signature = "LNLV", function(u){
    new(class(u), u = Matrix(1, nrow = u@dims, ncol = 1), dims = u@dims)
})
##
## One-element for points in SOC cones
setMethod("uone", signature = "SOCV", function(u){
    ans <- Matrix(0, nrow = u@dims, ncol = 1)
    ans[1, 1] <- 1.0
    new("SOCV", u = ans, dims = u@dims)
})
##
## One-element for points in PSD cones
setMethod("uone", signature = "PSDV", function(u){
    ans <- Diagonal(u@dims)
    dim(ans) <- c(u@dims^2, 1)
    new("PSDV", u = ans, dims = u@dims)
})
