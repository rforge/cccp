##
## Scaling of NNO-points by Hessian of Log-Barrier function  
setMethod("uslb", signature = c("NNOV", "NNOV"), function(u, lambda, inv = FALSE){
    if(inv){
        u@u <- u@u * lambda@u
    } else {
        u@u <- u@u / lambda@u
    }
    u
})
##
## Scaling of SOC-points by Hessian of Log-Barrier function  
setMethod("uslb", signature = c("SOCV", "SOCV"), function(u, lambda, inv = FALSE){
    a <- jnrm2(lambda)
    if(inv){
        lx <- crossprod(lambda@u, u@u) / a
    } else {
        lx <- jdot(lambda, u) / a
    }
    u0 <- u@u[1, 1]
    u@u[1, 1] <- lx
    cc <- drop((lx + u0) / (lambda@u[1, 1] / a + 1) / a)
    if(!inv){
        cc <- -cc
        a <- 1 / a
    }
    u@u[-1, 1] <- cc * lambda@u[-1, 1] + u@u[-1, 1]
    u@u <- a * u@u
    u
})
##
## Scaling of PSD-points by Hessian of Log-Barrier function  
setMethod("uslb", signature = c("PSDV", "PSDV"), function(u, lambda, inv = FALSE){
    m <- u@dims
    dim(lambda@u) <- c(m, m)
    dim(u@u) <- c(u@dims, u@dims)
    l <- diag(x = lambda@u)
    for(j in 1:m){
        lscaled <- sqrt(l[j]) * sqrt(l)
        if(inv){
            u@u[, j] <- u@u[, j] * lscaled
        } else {
            u@u[, j] <- u@u[, j] / lscaled
        }
    }
    dim(u@u) <- c(u@dims^2, 1)
    u
})
