##
## Update Nesterov-Todd scalings and lambda for variables pertinent to nonlinear constraints
setMethod("ntsu", signature = c("NLFS", "NLFV", "NLFV"), function(W, s, z){
    .ntsu_n(W@W, s@u, z@u)
})
##
## Update Nesterov-Todd scalings and lambda for variables in NNO cones
setMethod("ntsu", signature = c("NNOS", "NNOV", "NNOV"), function(W, s, z){
    .ntsu_l(W@W, s@u, z@u)
})
##
## Update Nesterov-Todd scalings and lambda for variables in SOC cones
setMethod("ntsu", signature = c("SOCS", "SOCV", "SOCV"), function(W, s, z){
    .ntsu_s(W@W, s@u, z@u)
})
##
## Update Nesterov-Todd scalings and lambda for variables in PSD cones
setMethod("ntsu", signature = c("PSDS", "PSDV", "PSDV"), function(W, s, z){
    .ntsu_p(W@W, s@u, z@u, s@dims)
})
