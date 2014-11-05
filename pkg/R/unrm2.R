##
## Norm of a vector
setMethod("unrm2", signature = "numeric", function(u){
    .snrm2_nls(u)
})
##
## Norm of a point in LNL cones
setMethod("unrm2", signature = "LNLV", function(u){
    .snrm2_nls(drop(u@u))
})
##
## Norm of a point in SOC cones
setMethod("unrm2", signature = "SOCV", function(u){
    .snrm2_nls(drop(u@u))
})
##
## Norm of a point in SOC cones
setMethod("unrm2", signature = "PSDV", function(u){
    .snrm2_p(drop(u@u), u@dims)
})

