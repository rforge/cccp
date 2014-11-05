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
    if(inv){
        w <- W@W[["rti"]]
        tt <- trans
    } else {
        w <- W@W[["r"]]
        if(trans == FALSE){
            tt <- TRUE
        } else {
            tt <- FALSE
        }
    }
    m <- ncol(w)
    for(j in 1:ncol(u)){
        U <- matrix(u[, j], nrow = m, ncol = m)
        diag(U) <- 0.5 * diag(U)
        U[upper.tri(U)] <- 0
        if(tt){
            a <- U %*% w
            ans <- crossprod(w, a) + crossprod(a, w)
        } else {
            a <- w %*% U
            ans <- tcrossprod(w, a) + tcrossprod(a, w)
        }
        dim(ans) <- c(m^2, 1)
        u[, j] <- ans
    }
    new("PSDV", u = u, dims = m)
})
##
## Nesterov-Todd scaling for NLC-points
## scaling with Matrix
setMethod("usnt", signature = c("matrix", "NLFS"), function(u, W, trans = FALSE, inv = FALSE){
    if(inv){
        w <- W@W[["dnli"]]
    } else {
        w <- W@W[["dnl"]]
    }
    for(j in 1:ncol(u)) u[, j] <- u[, j] * w
    new("NLFV", u = u, dims = nrow(u))
})
