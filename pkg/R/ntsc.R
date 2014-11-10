##
## Compute Nesterov-Todd scalings for variables in nonlinear constraints
setMethod("ntsc", signature = c("NLFV", "NLFV"), function(s, z){
    .ntsc_n(s@u, z@u)
})
## Initial scaling
setMethod("ntsc", signature = c("NLFC", "missing"), function(s, z){
    dnli <- dnl <- matrix(1, nrow = s@dims, ncol = 1)
    lambda <- new("NLFV", u = matrix(0, nrow = s@dims, ncol = 1), dims = s@dims)
    new("NLFS", W = list(dnl = dnl, dnli = dnli, lambda = lambda))
})
##
## Compute Nesterov-Todd scalings for variables in NNO cones
setMethod("ntsc", signature = c("NNOV", "NNOV"), function(s, z){
    .ntsc_l(s@u, z@u)
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
    .ntsc_s(s@u, z@u)
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
    .ntsc_p(s@u, z@u, s@dims)    
})
## Initial scaling (all variables (s, z, lambda) in column-stacked order)
setMethod("ntsc", signature = c("PSDC", "missing"), function(s, z){
    r <- rti <- diag(s@dims)
    lambda <- new("PSDV", u = matrix(0, nrow = s@dims^2, ncol = 1), dims = s@dims)
    new("PSDS", W = list(r = r, rti = rti, lambda = lambda))
})
