##
## cps-method for 'DEFLP'
setMethod("cps", signature = "DEFLP", function(cpd){
    ##
    ## setting constants & variables
    ##
    ctrl <- cpd@ctrl
    n <- cpd@n
    idx <- 1:cpd@k
    cc <- lapply(cpd@conecon, function(cc) cc@vclass)
    idxNSC <- integer(0)
    idxPSD <- which("PSDV" == cc)
    if(length(idxPSD) < 1){
        idxNSC <- idx
    } else {
        idxNSC <- idx[-idxPSD]
    }
    m <- sum(unlist(lapply(cpd@conecon, function(cc) cc@dims)))
    sdv <- rep(NA, 3)
    names(sdv) <- c("dg", "dgi", "lg")
    cvgdvals <- rep(NA, 6)
    names(cvgdvals) <- c("pobj", "dobj", "pinf", "dinf", "gap", "k/t")
    CurSol <- new("CPS")
    CurPdv <- SolKkt2 <- SolKkt1 <- new("PDV")
    hL <- lapply(cpd@conecon, function(x) x@h)
    hrz <- rz <- zeroList <- lapply(hL, function(h){
        h@u <- 0 * h@u
        h
    })
    kktSLV <- switch(ctrl@method,
                        solve = kktSOLVE(cpd)
                     )
    ##
    ## Computing initial values of PDV / CPS and scalings
    ##
    W <- lapply(cpd@conecon, ntsc)
    kktslv <- kktSLV(W, cpd)
    resx0 <- max(1, sqrt(udot(cpd@q)))
    resy0 <- max(1, sqrt(udot(drop(cpd@b))))
    resz0 <- max(1.0, sum(unlist(lapply(hL, function(z) unrm2(z)))))
    ## Primal start
    InitPrim <- with(kktslv@items, kktslv@f(x = rep(0, n), y = cpd@b, z = hL))
    InitPrim$s <- lapply(InitPrim$z, function(z){
        z@u <- -z@u
        z})
    ts <- max(unlist(lapply(InitPrim$s, function(s) umss(s))))
    ## Dual start
    InitDual <- with(kktslv@items, kktslv@f(x = -cpd@q, y = 0 * cpd@b, z = zeroList))
    tz <- max(unlist(lapply(InitDual$z, function(z) umss(z))))
    ## Initial point
    CurPdv@x <- InitPrim$x
    CurPdv@y <- InitDual$y
    CurPdv@s <- InitPrim$s
    CurPdv@z <- InitDual$z
    
    nrms <- sum(unlist(lapply(CurPdv@s, function(s) unrm2(s))))
    nrmz <- sum(unlist(lapply(CurPdv@z, function(z) unrm2(z))))
    ## Initial point optimal?
    gap <- sum(sapply(idx, function(j) udot(CurPdv@s[[j]], CurPdv@z[[j]])))
    pcost <- udot(CurPdv@x, cpd@q)
    dcost <- -drop(crossprod(cpd@b, CurPdv@y)) - sum(sapply(idx, function(j) udot(hL[[j]], CurPdv@z[[j]])))
    rgap <- NULL
    if(pcost < 0.0) rgap <- gap / -pcost
    if(dcost > 0.0) rgap <- gap / dcost
    checkRgap <- if(!is.null(rgap)){
        rgap <= ctrl@reltol
    } else {
        FALSE
    }
    if((ts <= 0.0) & (tz <= 0.0) &
       (gap <= ctrl@abstol || checkRgap)){
        resx <- unrm2(rdual(CurPdv, cpd))
        resy <- unrm2(rprim(CurPdv, cpd))
        resc <- rcent(CurPdv, cpd)
        resz <- sum(unlist(lapply(resc, unrm2)))
        pres <- max(resy / resy0, resz / resz0)
        dres <- resx / resx0
        CurSol@x <- CurPdv@x
        CurSol@y <- CurPdv@y
        CurSol@s <- CurPdv@s
        CurSol@z <- CurPdv@z
        CurSol@pobj <- pcost
        CurSol@dobj <- dcost
        CurSol@dgap <- gap
        if(!is.null(rgap)) CurSol@rdgap <- rgap
        CurSol@certp <- pres
        CurSol@certd <- dres
        CurSol@pslack <- -ts
        CurSol@dslack <- -tz
        CurSol@niter <- 0
        CurSol@status <- "optimal"
        if(ctrl@trace) cat("Optimal solution found.\n")
        return(CurSol)
    }
    if(ts >= -1e-8 * max(1.0, nrms)){
        CurPdv@s <- lapply(CurPdv@s, umsa, alpha = ts, init = TRUE)
    }
    
    if(tz >= -1e-8 * max(1.0, nrmz)){
        CurPdv@z <- lapply(CurPdv@z, umsa, alpha = tz, init = TRUE)
    }
        
    gap <- sum(sapply(idx, function(j) udot(CurPdv@s[[j]], CurPdv@z[[j]])))
    ##
    ## Start iterations
    ##
    for(i in 0:(ctrl@maxiters + 1)){
        ## Evaluate residuals, gap and stopping criteria
        ## Dual residuals
        hrx <- drop(-drop(crossprod(cpd@A, CurPdv@y)) - Reduce("+", lapply(idx, function(j) crossprod(cpd@conecon[[j]]@G, CurPdv@z[[j]]@u))))
        hresx <- sqrt(udot(hrx))
        rx <- hrx - cpd@q * CurPdv@tau
        resx <- sqrt(udot(rx)) / CurPdv@tau
        ## Primal residuals
        hry <- drop(cpd@A %*% CurPdv@x)
        hresy <- sqrt(udot(hry))
        ry <- hry - drop(cpd@b * CurPdv@tau)
        resy <- sqrt(udot(ry)) / CurPdv@tau
        ## Centrality residuals
        hrz <- lapply(idx, function(j){
            hrz[[j]]@u <- CurPdv@s[[j]]@u + cpd@conecon[[j]]@G %*% CurPdv@x
            hrz[[j]]
        })
        hresz <- sum(unlist(lapply(hrz, function(r) unrm2(r))))
        rz <- zeroList
        rz <- lapply(idx, function(j){
            rz[[j]]@u <- hrz[[j]]@u - hL[[j]]@u * CurPdv@tau
            rz[[j]]
        })
        resz <- sum(unlist(lapply(rz, function(r) unrm2(r)))) / CurPdv@tau
        ## self-dual residuals
        hz <- drop(sum(sapply(idx, function(j)
                              udot(hL[[j]], CurPdv@z[[j]]))))
        by <- drop(crossprod(cpd@b, CurPdv@y))
        cx <- udot(CurPdv@x, cpd@q)
        rt <- CurPdv@kappa + cx + by + hz
        ## Statistics for stopping criteria
        pcost <- cx / CurPdv@tau
        dcost <- -(by + hz) / CurPdv@tau
        rgap <- NULL
        if(pcost < 0.0) rgap <- gap / -pcost
        if(dcost > 0.0) rgap <- gap / dcost
        pres <- max(resy / resy0, resz / resz0)   
        dres <- resx / resx0
        if(hz + by < 0){
            pinfres <- hresx / resx0 / (-hz - by)
        } else {
            pinfres <- NULL
        }
        if(cx < 0){
            dinfres <- max(hresy / resy0, hresz / resz0) / (-cx)
        } else {
            dinfres <- NULL
        }
        ## tracing status quo of IPM
        if(ctrl@trace){
            cat(paste("Iteration:", i, "\n"))
            cvgdvals[1:6] <- signif(c(pcost,
                                      dcost,
                                      pres,
                                      dres,
                                      gap,
                                      CurPdv@kappa / CurPdv@tau))
            print(cvgdvals)
        }
        ## Checking convergence /infeasibilities
        checkRgap <- if(!is.null(rgap)){
            rgap <= ctrl@reltol
        } else {
            FALSE
        }
        if((pres <= ctrl@feastol) & (dres <= ctrl@feastol) &
           (gap <= ctrl@abstol || checkRgap)){
            CurSol@x <- CurPdv@x / CurPdv@tau
            CurSol@y <- CurPdv@y / CurPdv@tau
            CurSol@s <- lapply(CurPdv@s, function(s){
                s@u <- s@u / CurPdv@tau
                s
            })
            CurSol@z <- lapply(CurPdv@z, function(z){
                z@u <- z@u / CurPdv@tau
                z
            })
            ts <- max(unlist(lapply(CurSol@s, function(s) umss(s))))
            tz <- max(unlist(lapply(CurSol@z, function(z) umss(z))))
            CurSol@pobj <- pcost
            CurSol@dobj <- dcost
            CurSol@dgap <- gap
            if(!is.null(rgap)) CurSol@rdgap <- rgap
            CurSol@certp <- pres
            CurSol@certd <- dres
            CurSol@pslack <- -ts
            CurSol@dslack <- -tz
            CurSol@niter <- i
            CurSol@status <- "optimal"
            if(ctrl@trace) cat("Optimal solution found.\n")
            return(CurSol)
        } else if((!is.null(pinfres)) && (pinfres <= ctrl@feastol)){
            denom <- -hz - by
            CurSol@x <- numeric(0)
            CurSol@y <- CurPdv@y / denom
            CurSol@s <- list()
            CurSol@z <- lapply(CurPdv@z, function(z){
                z@u <- z@u / denom
                z
            })
            tz <- max(unlist(lapply(CurPdv@z, function(z) umss(z))))
            CurSol@pobj <- numeric(0)
            CurSol@dobj <- 1.0
            CurSol@dgap <- numeric(0)
            CurSol@rdgap <- numeric(0)
            CurSol@certp <- pinfres
            CurSol@certd <- numeric(0)
            CurSol@pslack <- numeric(0)
            CurSol@dslack <- -tz
            CurSol@niter <- i
            CurSol@status <- "primal infeasible"            
            if(ctrl@trace) cat("Certificate of primal infeasibility found.\n")
            return(CurSol)
        } else if((!is.null(dinfres)) && (dinfres <= ctrl@feastol)){
            denom <- -cx
            CurSol@x <- CurPdv@x / denom
            CurSol@y <- numeric(0)
            CurSol@s <- lapply(CurPdv@s, function(s){
                s@u <- s@u / denom
                s
            })
            CurSol@z <- list()
            ts <- max(unlist(lapply(CurSol@s, function(s) umss(s))))
            CurSol@pobj <- -1.0
            CurSol@dobj <- numeric(0)
            CurSol@dgap <- numeric(0)
            CurSol@rdgap <- numeric(0)
            CurSol@certp <- numeric(0)
            CurSol@certd <- dinfres
            CurSol@pslack <- -ts
            CurSol@dslack <- numeric(0)
            CurSol@niter <- i
            CurSol@status <- "dual infeasible"           
            if(ctrl@trace) cat("Certificate of dual infeasibility found.\n")
            return(CurSol)
        }
        ##
        ## Compute initial scalings
        ##
        if(i == 0){
            W <- sapply(idx, function(j) ntsc(CurPdv@s[[j]], CurPdv@z[[j]]))
            dg <- sqrt(CurPdv@kappa / CurPdv@tau)
            dgi <- sqrt(CurPdv@tau / CurPdv@kappa)
            lg <- sqrt(CurPdv@tau * CurPdv@kappa)
            sdv[1:3] <- c(dg, dgi, lg)
        }
        lambdasq <- lapply(idx, function(j) uprd(W[[j]]@W[["lambda"]]))
        lgsq <- sdv["lg"]^2
        ##
        ## Solution step 1 (same for affine and combined solution)
        ##
        kktslv <- kktSLV(W, cpd)
        ans1 <- try(with(kktslv@items, kktslv@f(x = -cpd@q, y = cpd@b, z = hL)))
        if(class(ans1) == "try-error"){
            CurSol@x <- CurPdv@x / CurPdv@tau
            CurSol@y <- CurPdv@y / CurPdv@tau
            CurSol@s <- lapply(CurSol@s, function(s){
                s@u <- s@u / CurPdv@tau
                s
            })
            CurSol@z <- lapply(CurSol@z, function(z){
                z@u <- z@u / CurPdv@tau
                z
            })
            ts <- max(unlist(lapply(CurSol@s, function(s) umss(s))))
            tz <- max(unlist(lapply(CurSol@z, function(z) umss(z))))
            CurSol@pobj <- pcost
            CurSol@dobj <- dcost
            CurSol@dgap <- gap
            if(!is.null(rgap)) CurSol@rdgap <- rgap
            CurSol@certp <- pres
            CurSol@certd <- dres
            CurSol@pslack <- -ts
            CurSol@dslack <- -tz
            CurSol@niter <- i
            CurSol@status <- "unknown"
            if(ctrl@trace) cat("Terminated (singular KKT matrix).\n")
            return(CurSol)
        }
        SolKkt1@x <- sdv["dgi"] * ans1$x
        SolKkt1@y <- sdv["dgi"] * ans1$y
        SolKkt1@z <- lapply(ans1$z, function(z){
            z@u <- sdv["dgi"] * z@u
            z})
        WhL <- lapply(idx, function(j) usnt(hL[[j]]@u, W[[j]], trans = TRUE, inv = TRUE))
        mu <- sum(unlist(lapply(idx, function(j) drop(crossprod(W[[j]]@W[["lambda"]]@u))))) / (m + 1)
        sigma <- 0
        ##
        ## Solving for affine and combined direction in two-round for-loop
        ##
        for(ii in 0:1){
            SolKkt2@s <- lambdasq
            SolKkt2@kappa <- lgsq
            if(ii == 1){
                SolKkt2@s <- lapply(idx, function(j){
                    SolKkt2@s[[j]]@u <- SolKkt2@s[[j]]@u + dsdz[[j]]@u - uone(SolKkt2@s[[j]])@u * sigma * mu
                    SolKkt2@s[[j]]                
                })
                SolKkt2@kappa <- SolKkt2@kappa + dkdt - sigma * mu
            }
            SolKkt2@x <- (1 - sigma) * rx
            SolKkt2@y <- (1 - sigma) * ry
            SolKkt2@z <- lapply(rz, function(z){
                z@u <- (1 - sigma) * z@u
                z
            })
            SolKkt2@tau <- (1 - sigma) * rt
            SolKkt2@x
            SolKkt2 <- kktSOL(cpd, SolKkt2, SolKkt1, W, WhL, sdv, kktslv, refine = ctrl@refine)
            SolKkt2@x
            if(ii == 0){
                ## ds o dz for Mehrotra correction
                dsdz <- lapply(idx, function(j) uprd(SolKkt2@s[[j]], SolKkt2@z[[j]]))
                ## dkappa * dtau
                dkdt <- SolKkt2@kappa * SolKkt2@tau
            }
            SolKkt2@s <- lapply(idx, function(j) uslb(u = SolKkt2@s[[j]], lambda = W[[j]]@W[["lambda"]]))
            SolKkt2@z <- lapply(idx, function(j) uslb(u = SolKkt2@z[[j]], lambda = W[[j]]@W[["lambda"]]))
            MaxStepS <- lapply(SolKkt2@s, umss)
            MaxStepZ <- lapply(SolKkt2@z, umss)
            ts <- max(unlist(lapply(MaxStepS, function(x) x$ms)))
            tz <- max(unlist(lapply(MaxStepZ, function(x) x$ms)))
            tt <- -SolKkt2@tau / sdv["lg"]
            tk <- -SolKkt2@kappa / sdv["lg"]
            tm <- max(c(0, ts, tz, tt, tk))
            if(tm == 0.0){
                step <- 1.0
            } else {
                if(ii == 0){
                    step <- min(1.0, 1.0 / tm)
                } else {
                    step <- min(1.0, ctrl@stepadj / tm)
                }
            }
            if(ii == 0) sigma <- (1.0 - step)^3
        }
        ##
        ## Updating x, y; s and z (in current scaling)
        ##
        CurPdv@x <- CurPdv@x + step * SolKkt2@x
        CurPdv@y <- CurPdv@y + step * SolKkt2@y

        if(length(idxNSC) > 0){
            for(j in idxNSC){
                SolKkt2@s[[j]] <- umsa(SolKkt2@s[[j]], alpha = step, init = FALSE)
                SolKkt2@z[[j]] <- umsa(SolKkt2@z[[j]], alpha = step, init = FALSE)
                SolKkt2@s[[j]] <- uslb(u = SolKkt2@s[[j]], lambda = W[[j]]@W[["lambda"]], inv = TRUE)
                SolKkt2@z[[j]] <- uslb(u = SolKkt2@z[[j]], lambda = W[[j]]@W[["lambda"]], inv = TRUE)
            }
        }

        if(length(idxPSD) > 0){
            for(j in idxPSD){
                s <- Matrix(MaxStepS[[j]]$evd$vectors)
                dim(s) <- c(SolKkt2@s[[j]]@dims^2, 1)
                SolKkt2@s[[j]]@u <- s
                sigs <- MaxStepS[[j]]$evd$values
                z <- Matrix(MaxStepZ[[j]]$evd$vectors)
                dim(z) <- c(SolKkt2@z[[j]]@dims^2, 1)
                SolKkt2@z[[j]]@u <- z
                sigz <- MaxStepZ[[j]]$evd$values

                SolKkt2@s[[j]] <- uslb(u = SolKkt2@s[[j]], lambda = W[[j]]@W[["lambda"]], inv = TRUE)
                SolKkt2@z[[j]] <- uslb(u = SolKkt2@z[[j]], lambda = W[[j]]@W[["lambda"]], inv = TRUE)

                SolKkt2@s[[j]] <- umsa(SolKkt2@s[[j]], alpha = step, init = FALSE, sigma = sigs, lambda = W[[j]]@W[["lambda"]])
                SolKkt2@z[[j]] <- umsa(SolKkt2@z[[j]], alpha = step, init = FALSE, sigma = sigz, lambda = W[[j]]@W[["lambda"]])                
            }
        }
        ##
        ## updating lambda and scaling
        ##
        W <- lapply(idx, function(j){
            ntsu(W = W[[j]], s = SolKkt2@s[[j]], z =  SolKkt2@z[[j]])
        })
    
        sdv["dg"] <- sdv["dg"] * sqrt(1 - step * tk) / sqrt(1 - step * tt) 
        sdv["dgi"] <- 1 / sdv["dg"]
        sdv["lg"] <- sdv["lg"] * sqrt(1 - step * tt) * sqrt(1 - step * tk)

        CurPdv@s <- lapply(idx, function(j){
            usnt(W[[j]]@W[["lambda"]]@u, W[[j]], inv = FALSE, trans = TRUE)
        })

        CurPdv@z <- lapply(idx, function(j){
            usnt(W[[j]]@W[["lambda"]]@u, W[[j]], inv = TRUE, trans = FALSE)
        })

        CurPdv@kappa <- sdv["lg"] / sdv["dgi"]
        CurPdv@tau <- sdv["lg"] * sdv["dgi"]

        gap <- (sqrt(sum(unlist(lapply(idx, function(j) drop(crossprod(W[[j]]@W[["lambda"]]@u)))))) / CurPdv@tau)^2
    } ## End for-loop
    
    ## Preparing CPS-object for non-convergence  
    
    CurSol@x <- CurPdv@x / CurPdv@tau
    CurSol@y <- CurPdv@y / CurPdv@tau
    CurSol@s <- lapply(CurPdv@s, function(s){
        s@u <- s@u / CurPdv@tau
        s
    })
    CurSol@z <- lapply(CurPdv@z, function(z){
        z@u <- z@u / CurPdv@tau
        z
    })
    ts <- max(unlist(lapply(CurSol@s, function(s) umss(s))))
    tz <- max(unlist(lapply(CurSol@z, function(z) umss(z))))
    CurSol@pobj <- pcost
    CurSol@dobj <- dcost
    CurSol@dgap <- gap
    if(!is.null(rgap)) CurSol@rdgap <- rgap
    CurSol@certp <- pres
    CurSol@certd <- dres
    CurSol@pslack <- -ts
    CurSol@dslack <- -tz
    CurSol@niter <- i
    CurSol@status <- "unknown"
    if(ctrl@trace){
        cat(paste("\n\n** Optimal solution could not be determined in", ctrl@maxiters, "iterations. **\n"))
    }
    return(CurSol)
})
##
## cps-method for 'DEFQP'
setMethod("cps", signature = "DEFQP", function(cpd){
    ## setting constants & variables
    ctrl <- cpd@ctrl
    n <- cpd@n
    CurSol <- new("CPS")
    CurPdv <- SolKkt <- new("PDV")
    if(cpd@k < 1){ ## no inequality constraints
        Pinv <- try(solve(cpd@P))
        if(class(Pinv) == "try-error") stop("Problem is unbounded below.\n")
        if(nrow(cpd@A) < 1){ ## neither equality consztraints
            CurSol@x <- drop(Pinv %*% cpd@q)
            CurSol@pobj <- pobj(CurSol@x, cpd)
            CurSol@status <- "optimal"
            return(CurSol)
        } else { ## only equality constraints
            PinvAt <- tcrossprod(Pinv, cpd@A)
            Pinvq <- Pinv %*% cpd@q
            S <- -cpd@A %*% PinvAt
            Sinv <- try(solve(S))
            if(class(Sinv) == "try-error") stop("Inversion of Schur complement failed.\n")
            CurSol@y <- CurPdv@y <- drop(Sinv %*% (cpd@A %*% Pinvq + cpd@b))
            CurSol@x <- CurPdv@x <- drop(Pinv %*% (-crossprod(cpd@A, CurSol@y) - cpd@q))
            CurSol@pobj <- pobj(CurPdv, cpd)
            CurSol@dobj <- dobj(CurPdv, cpd)
            ## primal infeasibilty
            nomin <- unrm2(drop(cpd@A %*% CurSol@x - cpd@b))
            denom <- max(1, unrm2(drop(cpd@b)))
            CurSol@certp <- nomin / denom
            nomin <- unrm2(drop(cpd@P %*% CurSol@x + crossprod(cpd@A, CurSol@y) + cpd@q))
            denom <- max(1, unrm2(cpd@q))
            CurSol@certd <- nomin / denom
            if((CurSol@certp <= ctrl@feastol) & (CurSol@certd <= ctrl@feastol)){
                CurSol@status <- "optimal"
            } else {
                CurSol@status <- "unknown"
            }
            return(CurSol)
        }
    }
    ##
    ## QPs with inequality constraints
    ##
    idx <- 1:cpd@k
    cc <- lapply(cpd@conecon, function(cc) cc@vclass)
    idxNSC <- integer(0)
    idxPSD <- which("PSDV" == cc)
    if(length(idxPSD) < 1){
        idxNSC <- idx
    } else {
        idxNSC <- idx[-idxPSD]
    }
    m <- sum(unlist(lapply(cpd@conecon, function(cc) cc@dims)))
    cvgdvals <- rep(NA, 5)
    names(cvgdvals) <- c("pobj", "dobj", "pinf", "dinf", "gap")
    hL <- lapply(cpd@conecon, function(x) x@h)
    hrz <- rz <- zeroList <- lapply(hL, function(h){
        h@u <- 0 * h@u
        h
    })
    kktSLV <- switch(ctrl@method,
                        solve = kktSOLVE(cpd)
                     )
    ##
    ## Computing initial values of PDV / CPS and scalings
    ##
    W <- lapply(cpd@conecon, ntsc)
    kktslv <- kktSLV(W, cpd)
    resx0 <- max(1, sqrt(udot(cpd@q)))
    resy0 <- max(1, sqrt(udot(drop(cpd@b))))
    resz0 <- max(1.0, sum(unlist(lapply(hL, function(z) unrm2(z)))))
    ## Initial point
    InitPdv <- with(kktslv@items, kktslv@f(x = -cpd@q, y = cpd@b, z = hL))
    CurPdv@x <- InitPdv$x
    CurPdv@y <- InitPdv$y
    CurPdv@z <- InitPdv$z
    CurPdv@s <- lapply(CurPdv@z, function(z){
        z@u <- -z@u
        z})

    ts <- max(unlist(lapply(CurPdv@s, function(s) umss(s))))    
    nrms <- sum(unlist(lapply(CurPdv@s, function(s) unrm2(s))))
    if(ts >= -1e-8 * max(1.0, nrms)){
        CurPdv@s <- lapply(CurPdv@s, umsa, alpha = ts, init = TRUE)
    }
    tz <- max(unlist(lapply(CurPdv@z, function(z) umss(z))))
    nrmz <- sum(unlist(lapply(CurPdv@z, function(z) unrm2(z)))) 
    if(tz >= -1e-8 * max(1.0, nrmz)){
        CurPdv@z <- lapply(CurPdv@z, umsa, alpha = tz, init = TRUE)
    }

    gap <- sum(sapply(idx, function(j) udot(CurPdv@s[[j]], CurPdv@z[[j]])))

    ##
    ## Start iterations
    ##
    for(i in 0:(ctrl@maxiters + 1)){
        ## Dual Residuals
        rx <- rdual(CurPdv, cpd)
        resx <- sqrt(udot(rx))
        ## Primal Residuals
        ry <- rprim(CurPdv, cpd)
        resy <- sqrt(udot(ry))
        ## Central Residuals
        rz <- rcent(CurPdv, cpd)
        resz <- sum(unlist(lapply(rz, function(cc) unrm2(cc))))
        ## Statistics for stopping criteria
        pcost <- pobj(CurPdv, cpd)
        dcost <- pcost + udot(CurPdv@y, ry) + sum(unlist(lapply(idx, function(j) udot(CurPdv@z[[j]], rz[[j]])))) - gap
        rgap <- NULL
        if(pcost < 0.0) rgap <- gap / -pcost
        if(dcost > 0.0) rgap <- gap / dcost
        pres <- max(resy / resy0, resz / resz0)   
        dres <- resx / resx0
        ## tracing status quo of IPM
        if(ctrl@trace){
            cat(paste("Iteration:", i, "\n"))
            cvgdvals[1:5] <- signif(c(pcost,
                                      dcost,
                                      pres,
                                      dres,
                                      gap))
            print(cvgdvals)
        }
        ## Checking convergence
        checkRgap <- if(!is.null(rgap)){
            rgap <= ctrl@reltol
        } else {
            FALSE
        }
        if((pres <= ctrl@feastol) & (dres <= ctrl@feastol) &
           (gap <= ctrl@abstol || checkRgap)){
            CurSol@x <- CurPdv@x 
            CurSol@y <- CurPdv@y
            CurSol@s <- CurPdv@s
            CurSol@z <- CurPdv@z
            ts <- max(unlist(lapply(CurSol@s, function(s) umss(s))))
            tz <- max(unlist(lapply(CurSol@z, function(z) umss(z))))
            CurSol@pobj <- pcost
            CurSol@dobj <- dcost
            CurSol@dgap <- gap
            if(!is.null(rgap)) CurSol@rdgap <- rgap
            CurSol@certp <- pres
            CurSol@certd <- dres
            CurSol@pslack <- -ts
            CurSol@dslack <- -tz
            CurSol@niter <- i
            CurSol@status <- "optimal"
            if(ctrl@trace) cat("Optimal solution found.\n")
            return(CurSol)
        }
        ##
        ## Compute initial scalings
        ##
        if(i == 0){
            W <- sapply(idx, function(j) ntsc(CurPdv@s[[j]], CurPdv@z[[j]]))
        }
        lambdasq <- lapply(idx, function(j) uprd(W[[j]]@W[["lambda"]]))
        kktslv <- kktSLV(W, cpd)
        mu <- gap / m
        sigma <- 0
        ##
        ## Solving for affine and combined direction in two-round for-loop
        ##
        for(ii in 0:1){
            if(length(idxNSC) > 0){
                SolKkt@s <- lapply(idxNSC, function(j){
                    s <- lambdasq[[j]]
                    s@u <- -lambdasq[[j]]@u + uone(CurPdv@s[[j]])@u * sigma * mu
                    s
                })
            }
            if(length(idxPSD) > 0){
                SolKkt@s <- lapply(idxNSC, function(j){
                    s <- lambdasq[[j]]
                    s@u <- -lambdasq[[j]]@u + uone(SolKkt@s[[j]])@u * sigma * mu
                    s
                })
            }
            SolKkt@x <- -rx
            SolKkt@y <- -ry
            SolKkt@z <- lapply(rz, function(z){
                z@u <- -z@u
                z
            })
            SolKkt <- kktSOL(cpd, SolKkt, W, kktslv, refine = ctrl@refine)

            ## ds o dz for Mehrotra correction
            dsdz <- sum(unlist(lapply(idx, function(j) udot(SolKkt@s[[j]], SolKkt@z[[j]]))))

            SolKkt@s <- lapply(idx, function(j) uslb(u = SolKkt@s[[j]], lambda = W[[j]]@W[["lambda"]]))
            SolKkt@z <- lapply(idx, function(j) uslb(u = SolKkt@z[[j]], lambda = W[[j]]@W[["lambda"]]))

            MaxStepS <- lapply(SolKkt@s, umss)
            MaxStepZ <- lapply(SolKkt@z, umss)

            ts <- max(unlist(lapply(MaxStepS, function(x) x$ms)))
            tz <- max(unlist(lapply(MaxStepZ, function(x) x$ms)))

            tm <- max(c(0, ts, tz))

            if(tm == 0.0){
                step <- 1.0
            } else {
                if(ii == 0){
                    step <- min(1.0, 1.0 / tm)
                } else {
                    step <- min(1.0, ctrl@stepadj / tm)
                }
            }
            if(ii == 0){
                sigma <- min(1.0, max(0.0, 1.0 - step + dsdz / gap * step^2))^3 
            }
        }

        ##
        ## Updating x, y; s and z (in current scaling)
        ##
        CurPdv@x <- CurPdv@x + step * SolKkt@x
        CurPdv@y <- CurPdv@y + step * SolKkt@y    

        if(length(idxNSC) > 0){
            for(j in idxNSC){
                SolKkt@s[[j]] <- umsa(SolKkt@s[[j]], alpha = step, init = FALSE)
                SolKkt@z[[j]] <- umsa(SolKkt@z[[j]], alpha = step, init = FALSE)
                SolKkt@s[[j]] <- uslb(u = SolKkt@s[[j]], lambda = W[[j]]@W[["lambda"]], inv = TRUE)
                SolKkt@z[[j]] <- uslb(u = SolKkt@z[[j]], lambda = W[[j]]@W[["lambda"]], inv = TRUE)
            }
        }

        if(length(idxPSD) > 0){
            for(j in idxPSD){
                s <- Matrix(MaxStepS[[j]]$evd$vectors)
                dim(s) <- c(SolKkt@s[[j]]@dims^2, 1)
                SolKkt@s[[j]]@u <- s
                sigs <- MaxStepS[[j]]$evd$values
                z <- Matrix(MaxStepZ[[j]]$evd$vectors)
                dim(z) <- c(SolKkt@z[[j]]@dims^2, 1)
                SolKkt@z[[j]]@u <- z
                sigz <- MaxStepZ[[j]]$evd$values

                SolKkt@s[[j]] <- uslb(u = SolKkt@s[[j]], lambda = W[[j]]@W[["lambda"]], inv = TRUE)
                SolKkt@z[[j]] <- uslb(u = SolKkt@z[[j]], lambda = W[[j]]@W[["lambda"]], inv = TRUE)

                SolKkt@s[[j]] <- umsa(SolKkt@s[[j]], alpha = step, init = FALSE, sigma = sigs, lambda = W[[j]]@W[["lambda"]])
                SolKkt@z[[j]] <- umsa(SolKkt@z[[j]], alpha = step, init = FALSE, sigma = sigz, lambda = W[[j]]@W[["lambda"]])                
            }
        }
        ##
        ## updating lambda and scaling
        ##
        W <- lapply(idx, function(j){
            ntsu(W = W[[j]], s = SolKkt@s[[j]], z =  SolKkt@z[[j]])
        })
    
        CurPdv@s <- lapply(idx, function(j){
            usnt(W[[j]]@W[["lambda"]]@u, W[[j]], inv = FALSE, trans = TRUE)
        })

        CurPdv@z <- lapply(idx, function(j){
            usnt(W[[j]]@W[["lambda"]]@u, W[[j]], inv = TRUE, trans = FALSE)
        })

        gap <- sum(unlist(lapply(idx, function(j) udot(W[[j]]@W[["lambda"]], W[[j]]@W[["lambda"]]))))       
    }
    
    ## Preparing CPS-object for non-convergence  
    
    CurSol@x <- CurPdv@x
    CurSol@y <- CurPdv@y
    CurSol@s <- CurPdv@s
    CurSol@z <- CurPdv@z
    ts <- max(unlist(lapply(CurSol@s, function(s) umss(s))))
    tz <- max(unlist(lapply(CurSol@z, function(z) umss(z))))
    CurSol@pobj <- pcost
    CurSol@dobj <- dcost
    CurSol@dgap <- gap
    if(!is.null(rgap)) CurSol@rdgap <- rgap
    CurSol@certp <- pres
    CurSol@certd <- dres
    CurSol@pslack <- -ts
    CurSol@dslack <- -tz
    CurSol@niter <- i
    CurSol@status <- "unknown"
    if(ctrl@trace){
        cat(paste("\n\n** Optimal solution could not be determined in", ctrl@maxiters, "iterations. **\n"))
    }
    return(CurSol)       
})

