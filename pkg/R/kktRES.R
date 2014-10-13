##
## Method for determining residuals of a KKT-system for LP programs
setMethod("kktRES", signature = c("DEFLP", "PDV", "PDV"), function(cpd, u, v, W, sdv){
    idx <- 1:cpd@k
    ## dual residual
    ## v@x := v@x - A'*u@y - G' * W^{-1} * u@z - c * u@tau / sdv["dg"]
    Ay <- drop(crossprod(cpd@A, u@y))
    Winvz <- lapply(idx, function(j) usnt(u@z[[j]]@u, W[[j]], inv = TRUE, trans = FALSE))
    GWinvz <- Reduce("+", lapply(idx, function(j) drop(crossprod(cpd@cList[[j]]@G, Winvz[[j]]@u))))
    qtaudg <- cpd@q * u@tau / sdv["dg"] 
    v@x <- v@x - Ay - GWinvz - qtaudg
    ## primal residual
    ## v@y := v@y + A * u@x - b * u@tau / sdv["dg"]
    v@y <- v@y + drop(cpd@A %*% u@x) - drop(cpd@b * u@tau / sdv["dg"])
    ## centrality residual (dual)
    ## v@z := v@z + G * u@x - h * u@tau / sdv["dg"] + W'u@s
    Wts <- lapply(idx, function(j) usnt(u@s[[j]]@u, W[[j]], inv = FALSE, trans = TRUE))
    v@z <- lapply(idx, function(j){
        vz <- v@z[[j]]
        vz@u <- vz@u + drop(cpd@cList[[j]]@G %*% u@x) - cpd@cList[[j]]@h@u * u@tau / sdv["dg"] + Wts[[j]]@u
        vz
    })
    ## centrality residual (primal)
    ## v@s := v@s + lambda o (u@z + u@s)
    v@s <- lapply(idx, function(j){
        zs <- u@s[[j]]
        zs@u <- zs@u + u@z[[j]]@u
        lzs <- uprd(W[[j]]@W[["lambda"]], zs)
        lzs@u <- v@s[[j]]@u + lzs@u
        lzs
    })                   
    ## self-dual residuals (tau, kappa)
    ## v@tau := v@tau + c' * u@x + b' * u@y + h' W^{-1} * u@z + sdv["dg"] * u@kappa
    hWz <- sum(unlist(lapply(idx, function(j) udot(cpd@cList[[j]]@h, Winvz[[j]]))))
    v@tau <- v@tau + udot(cpd@q, u@x) + drop(crossprod(cpd@b, u@y)) + hWz + sdv["dg"] * u@kappa
    ## v@kappa := v@kappa + sdv["lg"] * (u@tau + u@kappa)
    v@kappa <- v@kappa + sdv["lg"] * (u@tau + u@kappa)
    return(v)
})
##
## Method for determining residuals of a KKT-system for QP programs
setMethod("kktRES", signature = c("DEFQP", "PDV", "PDV"), function(cpd, u, v, W){
    idx <- 1:cpd@k
    ## dual residual
    ## v@x := v@x - P * ux - A'*u@y - G' * W^{-1} * u@z
    Px <- drop(crossprod(cpd@P, u@x))
    Ay <- drop(crossprod(cpd@A, u@y))
    Winvz <- lapply(idx, function(j) usnt(u@z[[j]]@u, W[[j]], inv = TRUE, trans = FALSE))
    GWinvz <- Reduce("+", lapply(idx, function(j) drop(crossprod(cpd@cList[[j]]@G, Winvz[[j]]@u))))
    v@x <- v@x - Px - Ay - GWinvz
    ## primal residual
    ## v@y := v@y - A * u@x 
    v@y <- v@y - drop(cpd@A %*% u@x)
    ## centrality residual (dual)
    ## v@z := v@z - G * u@x - W'u@s
    Wts <- lapply(idx, function(j) usnt(u@s[[j]]@u, W[[j]], inv = FALSE, trans = TRUE))
    v@z <- lapply(idx, function(j){
        vz <- v@z[[j]]
        vz@u <- vz@u - drop(cpd@cList[[j]]@G %*% u@x) - Wts[[j]]@u
        vz
    })
    ## centrality residual (primal)
    ## v@s := v@s - lambda o (u@z + u@s)
    v@s <- lapply(idx, function(j){
        zs <- u@s[[j]]
        zs@u <- zs@u + u@z[[j]]@u
        lzs <- uprd(W[[j]]@W[["lambda"]], zs)
        lzs@u <- v@s[[j]]@u - lzs@u
        lzs
    })                   
    return(v)
})
##
## Method for determining residuals of a KKT-system for LP programs with nonlinear constraints
setMethod("kktRES", signature = c("DEFNL", "PDV", "PDV"), function(cpd, u, v, W, sdv){
    idx <- 1:cpd@k
    ## dual residual
    ## v@x := v@x - H * ux - A'*u@y - G' * W^{-1} * u@z
    Hx <- drop(crossprod(cpd@H, u@x))
    Ay <- drop(crossprod(cpd@A, u@y))
    Winvz <- lapply(idx, function(j) usnt(u@z[[j]]@u, W[[j]], inv = TRUE, trans = FALSE))
    GWinvz <- Reduce("+", lapply(idx, function(j) drop(crossprod(cpd@cList[[j]]@G, Winvz[[j]]@u))))
    v@x <- v@x - Hx - Ay - GWinvz
    ## primal residual
    ## v@y := v@y - A * u@x 
    v@y <- v@y - drop(cpd@A %*% u@x)
    ## centrality residual (dual)
    ## v@z := v@z - G * u@x - W'u@s
    Wts <- lapply(idx, function(j) usnt(u@s[[j]]@u, W[[j]], inv = FALSE, trans = TRUE))
    v@z <- lapply(idx, function(j){
        vz <- v@z[[j]]
        vz@u <- vz@u - drop(cpd@cList[[j]]@G %*% u@x) - Wts[[j]]@u
        vz
    })
    ## centrality residual (primal)
    ## v@s := v@s - lambda o (u@z + u@s)
    v@s <- lapply(idx, function(j){
        zs <- u@s[[j]]
        zs@u <- zs@u + u@z[[j]]@u
        lzs <- uprd(W[[j]]@W[["lambda"]], zs)
        lzs@u <- v@s[[j]]@u - lzs@u
        lzs
    })                   
    return(v)  
})
##
## Method for determining residuals of a KKT-system for convex programs with nonlinear constraints
setMethod("kktRES", signature = c("DEFCP", "PDV", "PDV"), function(cpd, u, v, W, sdv){
    idx <- 1:cpd@k
    ## dual residual
    ## v@x := v@x - H * ux - A'*u@y - G' * W^{-1} * u@z
    Hx <- drop(crossprod(cpd@H, u@x))
    Ay <- drop(crossprod(cpd@A, u@y))
    Winvz <- lapply(idx, function(j) usnt(u@z[[j]]@u, W[[j]], inv = TRUE, trans = FALSE))
    GWinvz <- Reduce("+", lapply(idx, function(j) drop(crossprod(cpd@cList[[j]]@G, Winvz[[j]]@u))))
    v@x <- v@x - Hx - Ay - GWinvz
    ## primal residual
    ## v@y := v@y - A * u@x 
    v@y <- v@y - drop(cpd@A %*% u@x)
    ## centrality residual (dual)
    ## v@z := v@z - G * u@x - W'u@s
    Wts <- lapply(idx, function(j) usnt(u@s[[j]]@u, W[[j]], inv = FALSE, trans = TRUE))
    v@z <- lapply(idx, function(j){
        vz <- v@z[[j]]
        vz@u <- vz@u - drop(cpd@cList[[j]]@G %*% u@x) - Wts[[j]]@u
        vz
    })
    ## centrality residual (primal)
    ## v@s := v@s - lambda o (u@z + u@s)
    v@s <- lapply(idx, function(j){
        zs <- u@s[[j]]
        zs@u <- zs@u + u@z[[j]]@u
        lzs <- uprd(W[[j]]@W[["lambda"]], zs)
        lzs@u <- v@s[[j]]@u - lzs@u
        lzs
    })                   
    return(v)  
})
