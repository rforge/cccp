##
## Nesterov-Todd scaling for NNO-points
## scaling with Matrix
setMethod("usnt", signature = c("matrix", "NNOS"), function(u, W, trans = FALSE, inv = FALSE){
    u <- .ssnt_l(u, W@W, inv)
    new("NNOV", u = u, dims = nrow(u))
})
##
## Nesterov-Todd scaling for SOC-points
## scaling with Matrix
setMethod("usnt", signature = c("matrix", "SOCS"), function(u, W, trans = FALSE, inv = FALSE){
    u <- .ssnt_s(u, W@W, inv)
    new("SOCV", u = u, dims = nrow(u))
})
##
## Nesterov-Todd scaling for PSD-points
## scaling with Matrix (columns in stacked format)
setMethod("usnt", signature = c("matrix", "PSDS"), function(u, W, trans = FALSE, inv = FALSE){
    u <- .ssnt_p(u, W@W, trans, inv)
    new("PSDV", u = u, dims = ncol(W@W[["r"]]))
})
##
## Nesterov-Todd scaling for NLC-points
## scaling with Matrix
setMethod("usnt", signature = c("matrix", "NLFS"), function(u, W, trans = FALSE, inv = FALSE){
    u <- .ssnt_n(u, W@W, inv)
    new("NLFV", u = u, dims = nrow(u))
})
