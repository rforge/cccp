##
## Nesterov-Todd scaling for NNO-points
## scaling with Matrix
setMethod("usnt", signature = c("Matrix", "NNOS"), function(u, W, trans = FALSE, inv = FALSE){
    if(inv){
        w <- W@W[["di"]]
    } else {
        w <- W@W[["d"]]
    }
    for(j in 1:ncol(u)) u[, j] <- u[, j] * w
    new("NNOV", u = u, dims = nrow(u))
})
##
## Nesterov-Todd scaling for SOC-points
## scaling with Matrix
setMethod("usnt", signature = c("Matrix", "SOCS"), function(u, W, trans = FALSE, inv = FALSE){
    if(inv){
        u[1, ] <- -u[1, ]
    }
    w <- crossprod(u, W@W[["v"]])
    u[1,] <- -u[1,]
    u <- 2 * tcrossprod(W@W[["v"]], w) + u
    if(inv){
        u[1,] <- -u[1,]
        a <- 1 / W@W[["beta"]]
    } else {
        a <- W@W[["beta"]]
    }
    u <- a * u
    new("SOCV", u = u, dims = nrow(u))
})
##
## Nesterov-Todd scaling for PSD-points
## scaling with Matrix (columns in stacked format)
setMethod("usnt", signature = c("Matrix", "PSDS"), function(u, W, trans = FALSE, inv = FALSE){
    if(inv){
        w <- W@W[["rti"]]
        trans <- !trans
    } else {
        w <- W@W[["r"]]
    }
    m <- ncol(w)
    for(j in 1:ncol(u)){
        U <- Matrix(u[, j], nrow = m, ncol = m)
        if(trans){
            ans <- w %*% tcrossprod(U, w)
        } else {
            ans <- crossprod(w, U) %*% w
        }
        dim(ans) <- c(m^2, 1)
        u[, j] <- ans
    }
    new("PSDV", u = u, dims = m)
})
##
## Nesterov-Todd scaling for NLC-points
## scaling with Matrix
setMethod("usnt", signature = c("Matrix", "NLFS"), function(u, W, trans = FALSE, inv = FALSE){
    if(inv){
        w <- W@W[["dnli"]]
    } else {
        w <- W@W[["dnl"]]
    }
    for(j in 1:ncol(u)) u[, j] <- u[, j] * w
    new("NLFV", u = u, dims = nrow(u))
})
