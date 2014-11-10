##
## Scaling of LNL-points by Hessian of Log-Barrier function  
setMethod("uslb", signature = c("LNLV", "LNLV"), function(u, lambda, inv = FALSE){
    u@u <- .sslb_nl(u@u, lambda@u, inv)
    u
})
##
## Scaling of SOC-points by Hessian of Log-Barrier function  
setMethod("uslb", signature = c("SOCV", "SOCV"), function(u, lambda, inv = FALSE){
    u@u <- .sslb_s(u@u, lambda@u, inv)
    u
})
##
## Scaling of PSD-points by Hessian of Log-Barrier function  
setMethod("uslb", signature = c("PSDV", "PSDV"), function(u, lambda, inv = FALSE){
    u@u = .sslb_p(u@u, lambda@u, inv, u@dims)
    u
})
