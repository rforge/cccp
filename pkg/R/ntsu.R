##
## Update Nesterov-Todd scalings and lambda for variables in NNO cones
setMethod("ntsu", signature = c("NNOS", "NNOV", "NNOV"), function(W, s, z){
    s <- sqrt(s@u)
    z <- sqrt(z@u)
    d <- W@W[["d"]] * s / z
    di <- 1 / d
    l <- W@W[["lambda"]]
    l@u[, 1] <- s * z
    new("NNOS", W = list(d = d, di = di, lambda = l))
})
##
## Update Nesterov-Todd scalings and lambda for variables in SOC cones
setMethod("ntsu", signature = c("SOCS", "SOCV", "SOCV"), function(W, s, z){
    beta <- W@W[["beta"]]
    v <- W@W[["v"]]
    lambda <- W@W[["lambda"]]@u
    ln <- jnrm2(W@W[["lambda"]])
    aa <- jnrm2(s)
    s@u <- s@u / aa
    bb <- jnrm2(z)
    z@u <- z@u / bb
    cc <- sqrt((1 + udot(s, z)) / 2.0)
    vs <- udot(v[, 1], s@u[, 1])
    vz <- jdot(v[, 1], z@u[, 1])
    vq <- (vs + vz) / 2.0 / cc
    vu <- vs - vz
    lambda[1] <- cc
    wk0 <- 2 * v[1, 1] * vq - (s@u[1, 1] + z@u[1, 1]) / 2.0 / cc
    dd <- (v[1, 1] * vu - s@u[1, 1] / 2.0 + z@u[1, 1] / 2.0) / (wk0 + 1.0)
    lambda[-1, 1] <- v[-1, 1]
    lambda[-1, 1] <- 2.0 * (-dd * vq + 0.5 * vu) * lambda[-1, 1]
    lambda[-1, 1] <- 0.5 * (1.0 - dd / cc) * s@u[-1, 1] + lambda[-1, 1]
    lambda[-1, 1] <- 0.5 * (1.0 + dd / cc) * z@u[-1, 1] + lambda[-1, 1]
    lambda <- sqrt(aa * bb) * lambda
    l <- new("SOCV", u = lambda, dims = s@dims)
    v <- 2.0 * vq * v
    v[1, 1] <- v[1, 1] - s@u[1, 1] / 2.0 / cc
    v[-1, 1] <- 0.5 / cc * s@u[-1, 1] + v[-1, 1]
    v <- -0.5 / cc * z@u + v
    v[1, 1] <- v[1, 1] + 1.0
    v <- 1.0 / sqrt(2.0 * v[1, 1]) * v
    beta <- beta * sqrt(aa / bb)
    new("SOCS", W = list(v = v, beta = beta, lambda = l))
})
##
## Update Nesterov-Todd scalings and lambda for variables in PSD cones
setMethod("ntsu", signature = c("PSDS", "PSDV", "PSDV"), function(W, s, z){
    dim(s@u) <- c(s@dims, s@dims)
    dim(z@u) <- c(z@dims, z@dims)
    ZS <- crossprod(z@u, s@u)
    SVD <- svd(ZS)
    l <- diag(x = SVD$d)
    lisqrt <- diag(x = 1 / sqrt(SVD$d))
    r <- W@W[["r"]] %*% s@u %*% SVD$v %*% lisqrt
    rti <- W@W[["rti"]] %*% z@u %*% SVD$u %*% lisqrt
    W@W[["r"]] <- r 
    W@W[["rti"]] <- rti
    dim(l) <- c(s@dims^2, 1)
    W@W[["lambda"]]@u <- l
    W
})
##
## Update Nesterov-Todd scalings and lambda for variables pertinent to nonlinear constraints
setMethod("ntsu", signature = c("NLFS", "NLFV", "NLFV"), function(W, s, z){
    s <- sqrt(s@u)
    z <- sqrt(z@u)
    dnl <- W@W[["dnl"]] * s / z
    dnli <- 1 / dnl
    l <- W@W[["lambda"]]
    l@u[, 1] <- s * z
    new("NLFS", W = list(dnl = dnl, dnli = dnli, lambda = l))
})
