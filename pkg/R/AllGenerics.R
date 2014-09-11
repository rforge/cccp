##
## Generic for solving convex programs
setGeneric("cps", function(cpd, ...) standardGeneric("cps"))
##
## Generic for value of primal objective
setGeneric("pobj", function(pdv, cpd) standardGeneric("pobj"))
##
## Generic for value of dual objective
setGeneric("dobj", function(pdv, cpd) standardGeneric("dobj"))
##
## Generic for duality gap
setGeneric("dgap", function(pdv) standardGeneric("dgap"))
##
## Generic for centrality residuals
setGeneric("rcent", function(pdv, cpd) standardGeneric("rcent"))
##
## Generic for dual residuals
setGeneric("rdual", function(pdv, cpd) standardGeneric("rdual"))
##
## Generic for primal residuals
setGeneric("rprim", function(pdv, cpd) standardGeneric("rprim"))
##
## Generic for Nesterov-Todd scalings
setGeneric("ntsc", function(s, z) standardGeneric("ntsc"))
##
## Generic for updating Nesterov-Todd scalings and lambda
setGeneric("ntsu", function(W, s, z) standardGeneric("ntsu"))
##
## Generic for scaling constant of SOC-variables
setGeneric("jdot", function(u, v) standardGeneric("jdot"))
##
## Generic for the norm of cone-variable
setGeneric("jnrm2", function(u) standardGeneric("jnrm2"))
##
## Generic for inner-product between cone-variables
setGeneric("udot", function(u, v) standardGeneric("udot"))
##
## Generic for inverse of product between cone-variable
setGeneric("uinv", function(u, v) standardGeneric("uinv"))
##
## Generic for unpacking object in packed-storage form
setGeneric("umat", function(u) standardGeneric("umat"))
##
## Generic for adjusting cone variables by maximum step-size
setGeneric("umsa", function(u, ...) standardGeneric("umsa"))
##
## Generic for maximum step-size of cone-variables
setGeneric("umss", function(u) standardGeneric("umss"))
##
## Generic for the norm of cone-variable
setGeneric("unrm2", function(u) standardGeneric("unrm2"))
##
## Generic for one-element of cone-variable
setGeneric("uone", function(u) standardGeneric("uone"))
##
## Generic for product of cone-variable(s)
setGeneric("uprd", function(u, v) standardGeneric("uprd"))
##
## Generic for NT-scaling of cone constraints and variables
setGeneric("usnt", function(u, W, ...) standardGeneric("usnt"))
##
## Generic for Log-Barrier-scaling of cone constraints and variables
setGeneric("uslb", function(u, lambda, ...) standardGeneric("uslb"))
##
## Generic for packed-storage of symmetric matrices
setGeneric("uvec", function(u) standardGeneric("uvec"))
##
## Generic for extractor of x-variables
setGeneric("getx", function(object) standardGeneric("getx"))
##
## Generic for extractor of y-variables
setGeneric("gety", function(object) standardGeneric("gety"))
##
## Generic for extractor of s-variables
setGeneric("gets", function(object) standardGeneric("gets"))
##
## Generic for extractor of z-variables
setGeneric("getz", function(object) standardGeneric("getz"))
##
## Generic for computing a solution of a KKT-system
setGeneric("kktSOL", function(cpd, ...) standardGeneric("kktSOL"))
##
## Generic for computing residuals in KKT-system
setGeneric("kktRES", function(cpd, u, v, ...) standardGeneric("kktRES"))
##
## Generic for initial point of primal and dual slack-variables
setGeneric("initp", function(object) standardGeneric("initp"))
