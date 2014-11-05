##
## Compute Nesterov-Todd scalings for variables in NNO cones
setMethod("ntsc", signature = c("NNOV", "NNOV"), function(s, z){
    d <- sqrt(s@u / z@u)
    di <- sqrt(z@u / s@u)
    lambda <- new("NNOV", u = sqrt(s@u * z@u), dims = s@dims)
    new("NNOS", W = list(d = d, di = di, lambda = lambda))
})
## Initial scaling
setMethod("ntsc", signature = c("NNOC", "missing"), function(s, z){
    di <- d <- matrix(1, nrow = s@dims, ncol = 1)
    lambda <- new("NNOV", u = matrix(0, nrow = s@dims, ncol = 1), dims = s@dims)
    new("NNOS", W = list(d = d, di = di, lambda = lambda))
})
##
## Compute Nesterov-Todd scalings for variables in SOC cones
setMethod("ntsc", signature = c("SOCV", "SOCV"), function(s, z){
    aa <- jnrm2(s)
    bb <- jnrm2(z)
    beta <- sqrt(aa / bb)
    cc <- sqrt((udot(s, z) / aa / bb + 1.0) / 2.0)
    v <- -z@u / bb
    v[1] <- -v[1]
    v <- s@u / aa + v
    v <- 1.0 / 2.0 / cc * v
    v[1] <- v[1] + 1.0
    v <- 1.0 / sqrt(2.0 * v[1]) * v
    lambda <- matrix(0, nrow = s@dims, ncol = 1)
    lambda[1, 1] <- cc
    dd <- 2 * cc + s@u[1, 1] / aa + z@u[1, 1] / bb
    lambda[-1, 1] <- s@u[-1, 1]
    lambda[-1, 1] <- (cc + z@u[1, 1] / bb) / dd / aa * lambda[-1, 1]
    lambda[-1, 1] <- (cc + s@u[1, 1] / aa) / dd / bb * z@u[-1, 1] + lambda[-1, 1]
    lambda <- new("SOCV", u = sqrt(aa * bb) * lambda, dims = s@dims)
    new("SOCS", W = list(v = v, beta = beta, lambda = lambda))
})
## Initial scaling
setMethod("ntsc", signature = c("SOCC", "missing"), function(s, z){
    beta <- 1.0
    v <- matrix(0, nrow = s@dims, ncol = 1)
    v[1, 1] <- 1.0
    lambda <- new("SOCV", u = matrix(0, nrow = s@dims, ncol = 1), dims = s@dims)
    new("SOCS", W = list(v = v, beta = beta, lambda = lambda))
})
##
## Compute Nesterov-Todd scalings for variables in PSD cones
setMethod("ntsc", signature = c("PSDV", "PSDV"), function(s, z){
    dim(s@u) <- c(s@dims, s@dims)
    dim(z@u) <- c(z@dims, z@dims)
    Us <- chol(s@u)
    Us[lower.tri(Us)] <- 0
    Uz <- chol(z@u)
    UzUs <- tcrossprod(Uz, Us)
    SVD <- svd(UzUs)
    lsqrt <- sqrt(SVD$d)
    r <- solve(Uz) %*% SVD$u %*% diag(x = lsqrt)
    rti <- crossprod(Uz,SVD$u) %*% diag(x = 1 / lsqrt)
    l <- diag(x = SVD$d)
    dim(l) <- c(s@dims^2, 1)
    lambda <- new("PSDV", u = l, dims = s@dims)
    new("PSDS", W = list(r = r, rti = rti, lambda = lambda))
})
## Initial scaling (all variables (s, z, lambda) in column-stacked order)
setMethod("ntsc", signature = c("PSDC", "missing"), function(s, z){
    r <- rti <- diag(s@dims)
    lambda <- new("PSDV", u = matrix(0, nrow = s@dims^2, ncol = 1), dims = s@dims)
    new("PSDS", W = list(r = r, rti = rti, lambda = lambda))
})
##
## Compute Nesterov-Todd scalings for variables in nonlinear constraints
setMethod("ntsc", signature = c("NLFV", "NLFV"), function(s, z){
    dnl <- sqrt(s@u / z@u)
    dnli <- sqrt(z@u / s@u)
    lambda <- new("NLFV", u = sqrt(s@u * z@u), dims = s@dims)
    new("NLFS", W = list(dnl = dnl, dnli = dnli, lambda = lambda))
})
## Initial scaling
setMethod("ntsc", signature = c("NLFC", "missing"), function(s, z){
    dnli <- dnl <- matrix(1, nrow = s@dims, ncol = 1)
    lambda <- new("NLFV", u = matrix(0, nrow = s@dims, ncol = 1), dims = s@dims)
    new("NLFS", W = list(dnl = dnl, dnli = dnli, lambda = lambda))
})
