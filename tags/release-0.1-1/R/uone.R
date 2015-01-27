##
## One-element for points in LNL cones
setMethod("uone", signature = "LNLV", function(u){
    u@u <- matrix(.sone_nl(u@u), ncol = 1)
    u
})
##
## One-element for points in SOC cones
setMethod("uone", signature = "SOCV", function(u){
    u@u <- matrix(.sone_s(u@u), ncol = 1)
    u
})
##
## One-element for points in PSD cones
setMethod("uone", signature = "PSDV", function(u){
    u@u <- matrix(.sone_p(u@dims), ncol = 1)
    u
})
