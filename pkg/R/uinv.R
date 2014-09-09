##
## inverse of product of points in NNO cones
setMethod("uinv", signature = c("NNOV", "NNOV"), function(u, v){
    new("NNOV", u = u@u / v@u, dims = u@dims)
})
##
## inverse of product of points in SOC cones
setMethod("uinv", signature = c("SOCV", "SOCV"), function(u, v){
    aa <- jnrm2(v)^2
    cc <- u@u[1]
    dd <- udot(u@u[-1], v@u[-1])
    u@u[1] <- cc * v@u[1] - dd
    u@u[-1] <- aa / v@u[1] * u@u[-1]
    u@u[-1] <- (dd / v@u[1] - cc) * v@u[-1] + u@u[-1]
    ans <- u@u / aa
    new("SOCV", u = ans, dims = u@dims)
})
##
## inverse of product of points in PSD cones
setMethod("uinv", signature = c("PSDV", "PSDV"), function(u, v){
    Gamma <- matrix(NA, nrow = u@dims, ncol = u@dims)
    dim(u@u) <- c(u@dims, u@dims)
    dim(v@u) <- c(v@dims, v@dims)
    for(i in 1:u@dims){
        for(j in 1:u@dims){
            Gamma[i, j] <- 2 / (v@u[i, i] + v@u[j, j])
        }
    }
    ans <- u@u * Gamma 
    dim(ans) <- c(u@dims^2, 1)
    new("PSDV", u = ans, dims = u@dims)
})
