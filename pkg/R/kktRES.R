##
## Method for determining residuals of a KKT-system for LP programs
setMethod("kktRES", signature = c("DEFLP", "PDV", "PDV"), function(cpd, u, v, W, sdv){
    ResPdv <- new("PDV")
    idx <- 1:cpd@k
    ## dual residual
    ## v@x - A'*u@y - G' * w^{-1} * u@z - c * u@tau / sdv["dg"]
    Ay <- drop(crossprod(cpd@A, u@y))
    Winvz <- lapply(idx, function(j) usnt(u@z[[j]]@u, W[[j]], inv = TRUE, trans = FALSE))
    GWinvz <- Reduce("+", lapply(idx, function(j) drop(crossprod(cpd@conecon[[j]]@G, Winvz[[j]]@u))))
    qtaudg <- cpd@q * u@tau / sdv["dg"] 
    ResPdv@x <- v@x - Ay - GWinvz - qtaudg
    ## primal residual
    ## v@y + A * u@x - b * u@tau / sdv["dg"]
    ResPdv@y <- v@y + drop(cpd@A %*% u@x) - drop(cpd@b * u@tau / sdv["dg"])
    ## centrality residual (dual)
    ## v@z + G * u@x - h * u@tau / sdv["dg"] + W'u@s
    Wts <- lapply(idx, function(j) usnt(u@s[[j]]@u, W[[j]], inv = FALSE, trans = TRUE))
    ResPdv@z <- lapply(idx, function(j){
        vz <- v@z[[j]]
        vz@u <- vz@u + drop(cpd@conecon[[j]]@G %*% u@x) - cpd@conecon[[j]]@h@u * u@tau / sdv["dg"] + Wts[[j]]@u
        vz
    })
    ## centrality residual (primal)
    ## v@s + lambda o (u@z + u@s)
    ResPdv@s <- lapply(idx, function(j){
        zs <- u@s[[j]]
        zs@u <- zs@u + u@z[[j]]@u
        lzs <- uprd(W[[j]]@W[["lambda"]], zs)
        lzs@u <- v@z[[j]]@u + lzs@u
        lzs
    })                   
    ## self-dual residuals (tau, kappa)
    ## v@tau + c' * u@x + b' * u@y + h' W^{-1} * u@z + sdv["dg"] * u@kappa
    hWz <- sum(unlist(lapply(idx, function(j) udot(cpd@conecon[[j]]@h, Winvz[[j]]))))
    ResPdv@tau <- v@tau + udot(cpd@q, u@x) + drop(crossprod(cpd@b, u@y)) + hWz + sdv["dg"] * u@kappa
    ## v@kappa + sdv["lg"] * (u@tau + u@kappa)
    ResPdv@kappa <- v@kappa + sdv["lg"] * (u@tau + u@kappa)
    return(ResPdv)
})
##
## Method for determining residuals of a KKT-system for QP programs
setMethod("kktRES", signature = c("DEFQP", "PDV", "PDV"), function(cpd, u, v, W){
    ResPdv <- new("PDV")
    idx <- 1:cpd@k
    ## dual residual
    ## v@x - P * ux - A'*u@y - G' * w^{-1} * u@z
    Px <- drop(crossprod(cpd@P, u@x))
    Ay <- drop(crossprod(cpd@A, u@y))
    Winvz <- lapply(idx, function(j) usnt(u@z[[j]]@u, W[[j]], inv = TRUE, trans = FALSE))
    GWinvz <- Reduce("+", lapply(idx, function(j) drop(crossprod(cpd@conecon[[j]]@G, Winvz[[j]]@u))))
    ResPdv@x <- v@x - Px - Ay - GWinvz
    ## primal residual
    ## v@y - A * u@x 
    ResPdv@y <- v@y - drop(cpd@A %*% u@x)
    ## centrality residual (dual)
    ## v@z - G * u@x - W'u@s
    Wts <- lapply(idx, function(j) usnt(u@s[[j]]@u, W[[j]], inv = FALSE, trans = TRUE))
    ResPdv@z <- lapply(idx, function(j){
        vz <- v@z[[j]]
        vz@u <- vz@u - drop(cpd@conecon[[j]]@G %*% u@x) - Wts[[j]]@u
        vz
    })
    ## centrality residual (primal)
    ## v@s - lambda o (u@z + u@s)
    ResPdv@s <- lapply(idx, function(j){
        zs <- u@s[[j]]
        zs@u <- zs@u + u@z[[j]]@u
        lzs <- uprd(W[[j]]@W[["lambda"]], zs)
        lzs@u <- v@z[[j]]@u - lzs@u
        lzs
    })                   
    return(ResPdv)
})

