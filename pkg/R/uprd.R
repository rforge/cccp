##
## product of points in LNL cones
setMethod("uprd", signature = c("LNLV", "LNLV"), function(u, v){
    new(class(u), u = u@u * v@u, dims = u@dims)
})
setMethod("uprd", signature = c("LNLV", "missing"), function(u, v){
    ## similar to misc.ssqr
    new(class(u), u = u@u^2, dims = u@dims)
})
##
## product of points in SOC cones
setMethod("uprd", signature = c("SOCV", "SOCV"), function(u, v){
    ans <- u@u
    ans[1, 1] <- crossprod(u@u, v@u)
    ans[-1, 1] <- u@u[1, 1] * v@u[-1, 1] + v@u[1, 1] * u@u[-1, 1]
    new("SOCV", u = ans, dims = u@dims)
})
setMethod("uprd", signature = c("SOCV", "missing"), function(u, v){
    v <- u
    v@u[1, 1] <- crossprod(u@u)
    v@u[-1, 1] <- 2 * u@u[1, 1] * v@u[-1, 1]
    v
})
##
## product of points in PSD cones
setMethod("uprd", signature = c("PSDV", "PSDV"), function(u, v){
    dim(u@u) <- c(u@dims, u@dims)
    dim(v@u) <- c(v@dims, v@dims)
    ans <- 0.5 * (u@u %*% v@u + v@u %*% u@u)
    dim(ans) <- c(u@dims^2, 1)
    new("PSDV", u = ans, dims = u@dims)
})
setMethod("uprd", signature = c("PSDV", "missing"), function(u, v){
    dim(u@u) <- c(u@dims, u@dims)
    u@u <- crossprod(u@u)
    dim(u@u) <- c(u@dims^2, 1)
    u
})
