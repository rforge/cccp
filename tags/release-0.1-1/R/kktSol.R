##
## Method for determining a solution of a KKT-system for LP-programs
setMethod("kktSOL", signature = "DEFLP", function(cpd, SolKkt2, SolKkt1, W, WhL, sdv, kktslv, refine = FALSE){
    ## If iterative refinement
    if(refine){
        WVar <- SolKkt2
        SolKkt2 <- kktSOL(cpd, SolKkt2, SolKkt1, W, WhL, sdv, kktslv, refine = FALSE)
        for(i in 0:1){
            WVar2 <- WVar
            WVar2 <- kktRES(cpd, SolKkt2, WVar2, W, sdv)
            SolKkt2 <- kktSOL(cpd, WVar2, SolKkt1, W, WhL, sdv, kktslv, refine = FALSE)
        }
        return(SolKkt2)
    }
    ## Solving second KKT-sytem
    SolKkt2@s <- lapply(1:cpd@k, function(j) uinv(SolKkt2@s[[j]], W[[j]]@W[["lambda"]]))
    SolKkt2@s <- lapply(SolKkt2@s, function(s){
        s@u <- -s@u
        s})
    Ws3 <- lapply(1:cpd@k, function(j) usnt(SolKkt2@s[[j]]@u, W[[j]], inv = FALSE, trans = TRUE))
    SolKkt2@z <- lapply(1:cpd@k, function(j){
        SolKkt2@z[[j]]@u <- SolKkt2@z[[j]]@u + Ws3[[j]]@u
        SolKkt2@z[[j]]@u <- -SolKkt2@z[[j]]@u
        SolKkt2@z[[j]]
    })
    ans <- with(kktslv@items, kktslv@f(x = SolKkt2@x, y = -SolKkt2@y, z = SolKkt2@z))
    ## Combining solutions: SolKkt1 and SolKkt2
    SolKkt2@kappa <- -SolKkt2@kappa / sdv["lg"]
    SolKkt2@tau <- SolKkt2@tau + SolKkt2@kappa / sdv["dgi"]
    cx <- udot(cpd@q, ans$x)
    by <- drop(crossprod(cpd@b, ans$y))
    Whz <- sum(unlist(lapply(1:cpd@k, function(j) udot(WhL[[j]], ans$z[[j]]))))
    nomin <- sdv["dgi"] * (SolKkt2@tau + cx + by + Whz)
    denom <- 1.0 + sum(unlist(lapply(1:cpd@k, function(j) udot(SolKkt1@z[[j]]))))
    SolKkt2@tau <- nomin / denom
    SolKkt2@x <- ans$x + SolKkt2@tau * SolKkt1@x 
    SolKkt2@y <- ans$y + SolKkt2@tau * SolKkt1@y
    SolKkt2@z <- sapply(1:cpd@k, function(j){
        ans$z[[j]]@u <- ans$z[[j]]@u + SolKkt2@tau * SolKkt1@z[[j]]@u
        ans$z[[j]]
    })
    SolKkt2@s <- lapply(1:cpd@k, function(j){
        SolKkt2@s[[j]]@u <- SolKkt2@s[[j]]@u - SolKkt2@z[[j]]@u
        SolKkt2@s[[j]]
    })
    SolKkt2@kappa <- SolKkt2@kappa - SolKkt2@tau
    return(SolKkt2)
})
##
## Method for determining a solution of a KKT-system for QP-programs
setMethod("kktSOL", signature = "DEFQP", function(cpd, SolKkt, W, kktslv, refine = FALSE){
    ## If iterative refinement
    if(refine){
        WVar <- SolKkt
        SolKkt <- kktSOL(cpd, SolKkt, W, kktslv, refine = FALSE)
        for(i in 0:1){
            WVar2 <- WVar
            WVar2 <- kktRES(cpd, SolKkt, WVar2, W)
            SolKkt <- kktSOL(cpd, WVar2, W, kktslv, refine = FALSE)
        }
        return(SolKkt)
    }
    ## Solving KKT-sytem
    SolKkt@s <- lapply(1:cpd@k, function(j) uinv(SolKkt@s[[j]], W[[j]]@W[["lambda"]]))
    Ws3 <- lapply(1:cpd@k, function(j) usnt(SolKkt@s[[j]]@u, W[[j]], inv = FALSE, trans = TRUE))
    SolKkt@z <- lapply(1:cpd@k, function(j){
        SolKkt@z[[j]]@u <- SolKkt@z[[j]]@u - Ws3[[j]]@u
        SolKkt@z[[j]]
    })
    ans <- with(kktslv@items, kktslv@f(x = SolKkt@x, y = SolKkt@y, z = SolKkt@z))
    SolKkt@x <- ans$x 
    SolKkt@y <- ans$y
    SolKkt@z <- ans$z
    SolKkt@s <- lapply(1:cpd@k, function(j){
        SolKkt@s[[j]]@u <- SolKkt@s[[j]]@u - SolKkt@z[[j]]@u
        SolKkt@s[[j]]
    })
    return(SolKkt)
})
##
## Method for determining a solution of a KKT-system for LP-programs with nonlinear constraints
setMethod("kktSOL", signature = "DEFNL", function(cpd, SolKkt, W, kktslv, refine = FALSE){
    ## If iterative refinement
    if(refine){
        WVar <- SolKkt
        SolKkt <- kktSOL(cpd, SolKkt, W, kktslv, refine = FALSE)
        for(i in 0:1){
            WVar2 <- WVar
            WVar2 <- kktRES(cpd, SolKkt, WVar2, W)
            SolKkt <- kktSOL(cpd, WVar2, W, kktslv, refine = FALSE)
        }
        return(SolKkt)
    }
    ## Solving KKT-sytem
    SolKkt@s <- lapply(1:cpd@k, function(j) uinv(SolKkt@s[[j]], W[[j]]@W[["lambda"]]))
    Ws3 <- lapply(1:cpd@k, function(j) usnt(SolKkt@s[[j]]@u, W[[j]], inv = FALSE, trans = TRUE))
    SolKkt@z <- lapply(1:cpd@k, function(j){
        SolKkt@z[[j]]@u <- SolKkt@z[[j]]@u - Ws3[[j]]@u
        SolKkt@z[[j]]
    })
    ans <- with(kktslv@items, kktslv@f(x = SolKkt@x, y = SolKkt@y, z = SolKkt@z))
    SolKkt@x <- ans$x 
    SolKkt@y <- ans$y
    SolKkt@z <- ans$z
    SolKkt@s <- lapply(1:cpd@k, function(j){
        SolKkt@s[[j]]@u <- SolKkt@s[[j]]@u - SolKkt@z[[j]]@u
        SolKkt@s[[j]]
    })
    return(SolKkt)
})
##
## Method for determining a solution of a KKT-system for convex programs with nonlinear constraints
setMethod("kktSOL", signature = "DEFCP", function(cpd, SolKkt, W, kktslv, refine = FALSE){
    ## If iterative refinement
    if(refine){
        WVar <- SolKkt
        SolKkt <- kktSOL(cpd, SolKkt, W, kktslv, refine = FALSE)
        for(i in 0:1){
            WVar2 <- WVar
            WVar2 <- kktRES(cpd, SolKkt, WVar2, W)
            SolKkt <- kktSOL(cpd, WVar2, W, kktslv, refine = FALSE)
        }
        return(SolKkt)
    }
    ## Solving KKT-sytem
    SolKkt@s <- lapply(1:cpd@k, function(j) uinv(SolKkt@s[[j]], W[[j]]@W[["lambda"]]))
    Ws3 <- lapply(1:cpd@k, function(j) usnt(SolKkt@s[[j]]@u, W[[j]], inv = FALSE, trans = TRUE))
    SolKkt@z <- lapply(1:cpd@k, function(j){
        SolKkt@z[[j]]@u <- matrix(SolKkt@z[[j]]@u - Ws3[[j]]@u)
        SolKkt@z[[j]]
    })
    ans <- with(kktslv@items, kktslv@f(x = SolKkt@x, y = SolKkt@y, z = SolKkt@z))
    SolKkt@x <- ans$x 
    SolKkt@y <- ans$y
    SolKkt@z <- ans$z
    SolKkt@s <- lapply(1:cpd@k, function(j){
        SolKkt@s[[j]]@u <- SolKkt@s[[j]]@u - SolKkt@z[[j]]@u
        SolKkt@s[[j]]
    })
    return(SolKkt)
})
