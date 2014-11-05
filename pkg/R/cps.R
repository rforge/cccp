##
## cps-method for 'DEFLP'
setMethod("cps", signature = "DEFLP", function(cpd){
    ##
    ## setting constants & variables
    ##
    ctrl <- cpd@ctrl
    n <- cpd@n
    idx <- 1:cpd@k
    cc <- lapply(cpd@cList, function(cc) cc@vclass)
    idxNSC <- integer(0)
    idxPSD <- which("PSDV" == cc)
    if(length(idxPSD) < 1){
        idxNSC <- idx
    } else {
        idxNSC <- idx[-idxPSD]
    }
    m <- sum(unlist(lapply(cpd@cList, function(cc) cc@dims)))
    sdv <- rep(NA, 3)
    names(sdv) <- c("dg", "dgi", "lg")
    cvgdvals <- rep(NA, 6)
    names(cvgdvals) <- c("pobj", "dobj", "pinf", "dinf", "gap", "k/t")
    CurSol <- new("CPS")
    CurPdv <- SolKkt2 <- SolKkt1 <- new("PDV")
    hL <- lapply(cpd@cList, function(x) x@h)
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
    W <- lapply(cpd@cList, ntsc)
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
        hrx <- drop(-drop(crossprod(cpd@A, CurPdv@y)) - Reduce("+", lapply(idx, function(j) crossprod(cpd@cList[[j]]@G, CurPdv@z[[j]]@u))))
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
            hrz[[j]]@u <- CurPdv@s[[j]]@u + cpd@cList[[j]]@G %*% CurPdv@x
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
            SolKkt2 <- kktSOL(cpd, SolKkt2, SolKkt1, W, WhL, sdv, kktslv, refine = ctrl@refine)
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
                s <- matrix(MaxStepS[[j]]$evd$vectors)
                dim(s) <- c(SolKkt2@s[[j]]@dims^2, 1)
                SolKkt2@s[[j]]@u <- s
                sigs <- MaxStepS[[j]]$evd$values
                z <- matrix(MaxStepZ[[j]]$evd$vectors)
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
    cc <- lapply(cpd@cList, function(cc) cc@vclass)
    idxNSC <- integer(0)
    idxPSD <- which("PSDV" == cc)
    if(length(idxPSD) < 1){
        idxNSC <- idx
    } else {
        idxNSC <- idx[-idxPSD]
    }
    m <- sum(unlist(lapply(cpd@cList, function(cc) cc@dims)))
    cvgdvals <- rep(NA, 5)
    names(cvgdvals) <- c("pobj", "dobj", "pinf", "dinf", "gap")
    hL <- lapply(cpd@cList, function(x) x@h)
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
    W <- lapply(cpd@cList, ntsc)
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
            SolKkt@s <- lapply(idx, function(j){
                s <- lambdasq[[j]]
                s@u <- -lambdasq[[j]]@u + uone(CurPdv@s[[j]])@u * sigma * mu
                s
            })            
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
                s <- matrix(MaxStepS[[j]]$evd$vectors)
                dim(s) <- c(SolKkt@s[[j]]@dims^2, 1)
                SolKkt@s[[j]]@u <- s
                sigs <- MaxStepS[[j]]$evd$values
                z <- matrix(MaxStepZ[[j]]$evd$vectors)
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
##
## cps-method for 'DEFNL'
setMethod("cps", signature = "DEFNL", function(cpd){
    ## setting constants & variables
    ctrl <- cpd@ctrl
    n <- cpd@n
    mnl <- cpd@mnl
    m <- sum(unlist(lapply(cpd@cList, function(cc) cc@dims)))
    MaxRelaxedIters <- ctrl@maxreliter
    RelaxedIters <- 0L
    kktSLV <- switch(ctrl@method,
                        solve = kktSOLVE(cpd)
                     )
    idx <- 1:cpd@k
    cc <- lapply(cpd@cList, function(cc) cc@vclass)
    idxNSC <- integer(0)
    idxPSD <- which("PSDV" == cc)
    if(length(idxPSD) < 1){
        idxNSC <- idx
    } else {
        idxNSC <- idx[-idxPSD]
    }
    ConeCon <- ifelse(cpd@k > 1, TRUE, FALSE)
    cvgdvals <- rep(NA, 5)
    names(cvgdvals) <- c("pobj", "dobj", "pinf", "dinf", "gap")
    CurSol <- new("CPS")
    CurPdv <- NewPdv <- SolKkt <- new("PDV")
    CurPdv@x <- cpd@x0
    CurPdv@y <- rep(0, nrow(cpd@A))
    CurPdv@s <- CurPdv@z <- lapply(cpd@cList, initp)
    ##
    ## Start iterations
    ##
    for(i in 0:(ctrl@maxiters + 1)){
        ## Setting Df to 'G' slot of NLFC
        cpd@cList[[1]]@G <- do.call("rbind", lapply(cpd@nlfList, grad, x = CurPdv@x))
        ## Setting f to 'h' slot of NLFC
        cpd@cList[[1]]@h@u <- matrix(unlist(lapply(
            cpd@nlfList, function(nlf) nlf(CurPdv@x))),
                                     nrow = cpd@mnl, ncol = 1)
        ## Computing Hessian and setting to cpd@H
        cpd@H <- Reduce("+", lapply(1:cpd@mnl, function(j){
            matrix(CurPdv@z[[1]]@u[j, 1] * hessian(func = cpd@nlfList[[j]], x = CurPdv@x))
        }))
        ## Computing gap
        gap <- sum(sapply(idx, function(j) udot(CurPdv@s[[j]], CurPdv@z[[j]])))
        ## Computing residuals
        ## Dual Residuals
        rx <- rdual(CurPdv, cpd)
        resx <- sqrt(udot(rx))
        ## Primal Residuals
        ry <- rprim(CurPdv, cpd)
        resy <- sqrt(udot(ry))
        ## Central Residuals
        rz <- rcent(CurPdv, cpd)
        resznl <- unrm2(rz[[1]])
        resz <- sum(unlist(lapply(rz, function(cc) unrm2(cc))))
        ## Statistics for stopping criteria
        pcost <- pobj(CurPdv, cpd)
        dcost <- pcost + udot(CurPdv@y, ry) + sum(unlist(lapply(idx, function(j)
                                                                udot(CurPdv@z[[j]], rz[[j]])))) - gap
        rgap <- NULL
        if(pcost < 0.0) rgap <- gap / -pcost
        if(dcost > 0.0) rgap <- gap / dcost
        pres <- sqrt(resy^2 + resz^2)
        dres <- resx
        if(i == 0){
            resx0 <- max(1.0, resx)
            resznl0 <- max(1.0, resznl)
            pres0 <- max(1.0, pres)
            dres0 <- max(1.0, dres)
            gap0 <- gap
            theta1 <- 1.0 / gap0
            theta2 <- 1.0 / resx0
            theta3 <- 1.0 / resznl0
        }
        phi <- theta1 * gap + theta2 * resx + theta3 * resznl
        pres <- pres / pres0
        dres <- dres / dres0

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
        sigma <- eta <- 0.0
        kktslv <- kktSLV(W, cpd)
        ##
        ## Solving for affine and combined direction in two-round for-loop
        ##
        for(ii in 0:1){
            mu <- gap / m
            SolKkt@s <- lapply(idx, function(j){
                s <- lambdasq[[j]]
                s@u <- -lambdasq[[j]]@u + uone(CurPdv@s[[j]])@u * sigma * mu
                s
            })
            SolKkt@x <- -(1 - eta) * rx
            SolKkt@y <- -(1 - eta) * ry
            SolKkt@z <- lapply(rz, function(z){
                z@u <- -(1 - eta) * z@u
                z
            })
            SolKkt <- try(kktSOL(cpd, SolKkt, W, kktslv, refine = ctrl@refine))
            if(class(SolKkt) == "try-error"){
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
                if(ctrl@trace) cat("Terminated (singular KKT matrix).\n")
                return(CurSol)
            }
            
            ## Inner product ds'*dz and unscaled steps used in line search.
            dsdz <- sum(unlist(lapply(idx, function(j) udot(SolKkt@s[[j]], SolKkt@z[[j]]))))
            ## unscaling slack-variables
            dsu <- lapply(idx, function(j) usnt(SolKkt@s[[j]]@u, W[[j]], inv = FALSE, trans = TRUE))
            dzu <- lapply(idx, function(j) usnt(SolKkt@z[[j]]@u, W[[j]], inv = TRUE, trans = FALSE))
            ## Maximum step to boundary
            SolKkt@s <- lapply(idx, function(j) uslb(u = SolKkt@s[[j]], lambda = W[[j]]@W[["lambda"]]))
            SolKkt@z <- lapply(idx, function(j) uslb(u = SolKkt@z[[j]], lambda = W[[j]]@W[["lambda"]]))

            MaxStepS <- lapply(SolKkt@s, umss)
            MaxStepZ <- lapply(SolKkt@z, umss)

            ts <- max(unlist(lapply(MaxStepS, function(x) Re(x$ms))))
            tz <- max(unlist(lapply(MaxStepZ, function(x) Re(x$ms))))
            tm <- max(c(0, ts, tz))

            if(tm == 0.0){
                step <- 1.0
            } else {
                step <- min(1.0, ctrl@stepadj / tm)
            }

            ## Backtracking until x is in the domain of f
            backtrack <- TRUE
            while(backtrack){
                x <- CurPdv@x + step * SolKkt@x
                Fval <- unlist(lapply(cpd@nlfList, function(f) f(x)))
                if(any(is.nan(Fval))){
                    newFval <- NULL
                } else {
                    newFval <- matrix(Fval, nrow = mnl, ncol = 1)
                    newDF <- matrix(do.call("rbind", lapply(cpd@nlfList, grad, x = x)))
                    backtrack = FALSE
                }
                step <- step * ctrl@beta
            }

            ## Merit function
            phi <- theta1 * gap + theta2 * resx + theta3 * resznl
            if(ii == 0){
                dphi <- -phi
            } else {
                dphi <- -theta1 * (1 - sigma) * gap - theta2 * (1 - eta) * resx - theta3 * (1 - eta) * resznl
            }

            ## Line search
            backtrack = TRUE
            while(backtrack){
                NewPdv@x <- CurPdv@x + step * SolKkt@x
                NewPdv@y <- CurPdv@y + step * SolKkt@y
                NewPdv@z <- lapply(idx, function(j){
                    ans <- dzu[[j]]
                    ans@u <- CurPdv@z[[j]]@u + step * dzu[[j]]@u
                    ans
                })
                NewPdv@s <- lapply(idx, function(j){
                    ans <- dsu[[j]]
                    ans@u <- CurPdv@s[[j]]@u + step * dsu[[j]]@u
                    ans
                })
                               
                cpd@cList[[1]]@h@u <- matrix(unlist(lapply(
                    cpd@nlfList, function(f) f(NewPdv@x))), ncol = 1) ## newf
                cpd@cList[[1]]@G <- do.call("rbind", lapply(
                    cpd@nlfList, grad, x = NewPdv@x)) ## newDf

                ## Residuals
                ## Dual Residuals
                newrx <- rdual(NewPdv, cpd)
                newresx <- sqrt(udot(newrx))
                ## Central Residuals of nonlinear constraints
                newrznl <- drop(NewPdv@s[[1]]@u + cpd@cList[[1]]@h@u)
                newresznl <- unrm2(newrznl)

                newgap <- (1.0 - (1.0 - sigma) * step) * gap + step^2 * dsdz
                newphi <- theta1 * newgap + theta2 * newresx + theta3 * newresznl

                if(ii == 0){
                    check1 <- newgap <= (1.0 - ctrl@alpha * step) * gap
                    check2 <- (RelaxedIters >= 0) && (RelaxedIters < MaxRelaxedIters)
                    check3 <- newphi <= phi + ctrl@alpha * step * dphi
                    if(check1 && (check2 || check3)){
                        backtrack = FALSE
                        sigma <- min(newgap / gap, (newgap / gap)^3)
                        eta <- 0.0
                    } else {
                        step <- step * ctrl@beta
                    }
                } else {
                    if(RelaxedIters == -1L || ((RelaxedIters == 0L) && (MaxRelaxedIters == 0L))){
                        ## Do a standard line search
                        check3 <- newphi <= phi + ctrl@alpha * step * dphi
                        if(check3) {
                            RelaxedIters <- 0L
                            backtrack <- FALSE
                        } else {
                            step <- step * ctrl@beta
                        }
                    } else if(RelaxedIters == 0 && (RelaxedIters < MaxRelaxedIters)){
                        check3 <- newphi <= phi + ctrl@alpha * step * dphi
                        if(check3){
                            ## Relaxed line serach gives sufficient decreaase
                            RelaxedIters <- 0L
                        } else {
                            ## save state
                            cmnl <- cpd@cList[[1]]
                            phi0 <- phi
                            dphi0 <- dphi
                            gap0 <- gap
                            step0 <- step
                            W0 <- W
                            CurPdv0 <- CurPdv
                            lambdasq0 <- lambdasq
                            dsdz0 <- dsdz
                            sigma0 <- sigma
                            eta0 <- eta
                            rx0 <- rx
                            ry0 <- ry
                            rz0 <- rz
                            RelaxedIters <- 1L
                        }
                    backtrack <- FALSE
                } else if((RelaxedIters >= 0L) && (RelaxedIters < MaxRelaxedIters) && (MaxRelaxedIters > 0L)){
                    if(newphi <= phi0 + ctrl@alpha * step0 * dphi0){
                        ## Relaxed line search gives sufficient decrease
                        RelaxedIters <- 0L
                    } else {
                        ## Relaxed line search
                        RelaxedIters <- RelaxedIters + 1L
                    }
                    backtrack <- FALSE
                } else if(RelaxedIters == MaxRelaxedIters && (MaxRelaxedIters > 0)){
                    if(newphi <= phi0 + ctrl@alpha * step0 * dphi0){
                        ## Series of relaxed line searches ends with
                        ## sufficient decrease w.r.t. phi0
                        backtrack <- FALSE
                        RelaxedIters <- 0L
                    } else if(newphi >= phi0){
                        ## Resume last saved line search
                        cpd@cList[[1]] <- cmnl
                        phi <- phi0
                        dphi <- dphi0
                        gap <- gap0
                        step <- step0
                        W <- W0
                        CurPdv <- CurPdv0
                        lambdasq <- lambdasq0            
                        dsdz <- dsdz0
                        sigma <- sigma0
                        eta <- eta0
                        RelaxedIters <- -1L
                    } else if(newphi <= phi + ctrl@alpha * step * dphi)
                        ## Series of relaxed line seraches ends with
                        ## insufficient decrease w.r.t. phi0
                        backtrack <- FALSE
                        RelaxedIters <- -1L
                    }
                }
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
                s <- matrix(Re(MaxStepS[[j]]$evd$vectors))
                dim(s) <- c(SolKkt@s[[j]]@dims^2, 1)
                SolKkt@s[[j]]@u <- s
                sigs <- Re(MaxStepS[[j]]$evd$values)
                z <- matrix(Re(MaxStepZ[[j]]$evd$vectors))
                dim(z) <- c(SolKkt@z[[j]]@dims^2, 1)
                SolKkt@z[[j]]@u <- z
                sigz <- Re(MaxStepZ[[j]]$evd$values)

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
