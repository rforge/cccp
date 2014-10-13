##
## Solving KKT-system for convex programs with nonlinear objective
kktCP <- function(cpd){
    kktSLV <- function(W, cpd){
        ldK <- sum(dim(cpd@A)) - 1
        K <- Matrix(0, nrow = ldK, ncol = ldK)
        K <- as(K, "dsyMatrix")
        n <- cpd@n
        ne <- n + 1L
        idx <- 1:cpd@k
        K[1:n, 1:n] <- cpd@H ## Hessian: either 0 or z(j) H_j for non-linear constraints
        K[-(1:n), 1:n] <- cpd@A[, -ne]
        K[1:n, -(1:n)] <- t(cpd@A[, -ne])
        We <- W
        if(cpd@mnl > 1){
            We[[1]]@W$dnl <- We[[1]]@W$dnl[-1] 
            We[[1]]@W$dnli <- We[[1]]@W$dnli[-1] 
        } else {
            idx <- idx[-1]
        }
        GGList <- lapply(cpd@cList, function(cc) cc@G)
        GList <- lapply(GGList, function(G) Matrix(G[, -ne]))
        if(cpd@mnl > 1){
            GList[[1]] <- GList[[1]][-1, ] ## removing first row pertinent to f0
        } 
        for(j in idx){
            WitG <- usnt(GList[[j]], We[[j]], trans = TRUE, inv = TRUE)
            WiWitG <- usnt(WitG@u, We[[j]], trans = FALSE, inv = TRUE)
            GtWiWitG <- crossprod(GList[[j]], WiWitG@u)
            K[1:n, 1:n] <- K[1:n, 1:n] + GtWiWitG
        }            
        kktslv <- function(x, y, z){
            if(!is.list(z)) z <- list(z)
            a <- z[[1]]@u[1] ## slack with respect to f0
            ux <- x[-ne]
            ux
            ux <- ux + x[ne] * GGList[[1]][1, -ne]
            uz <- z
            idx <- 1:cpd@k
            if(cpd@mnl > 1){
                uz[[1]]@u <- uz[[1]]@u[-1] ## removing slack pertinent to f0
            } else {
                idx <- idx[-1]
            }
            We <- W
            if(cpd@mnl > 1){
                We[[1]]@W@dnl <- We[[1]]@W$dnl[-1] ## removing scaling pertinent to f0
                We[[1]]@W@dnli <- We[[1]]@W$dnli[-1] ## removing inverse scaling pertinent to f0
            } 
            GList <- lapply(GGList, function(G) Matrix(G[, -ne]))
            if(cpd@mnl > 1){
                GList[[1]] <- GList[[1]][-1, ] ## removing first row pertinent to f0
            }

            if(length(idx) > 0){
                GtWiWitz <- Reduce("+", lapply(idx, function(j){
                    Witz <- usnt(uz[[j]]@u, We[[j]], inv = TRUE, trans = TRUE)
                    WiWitz <- usnt(Witz@u, We[[j]], inv = TRUE, trans = FALSE)
                    GtWiWitz <- crossprod(GList[[j]], WiWitz@u)
                }))
                RHS <- c(ux + drop(GtWiWitz), drop(y))
            } else {
                RHS <- c(ux, drop(y))
            }
            ans <- solve(K, RHS)
            x0 <- ans[1:n, 1] ## solution to 'x' except 't'
            y <- ans[-(1:n), 1]
            z <- lapply(idx, function(j){
                uz <- new(cpd@cList[[j]]@vclass, u = GList[[j]] %*% x0 - uz[[j]]@u, dims = cpd@cList[[j]]@dims)
                usnt(uz@u, We[[j]], inv = TRUE, trans = TRUE)
            })

            if(cpd@mnl > 1){
                z[[1]]@u <- Matrix(c(-x[ne] * W[[1]]@W$dnl[1], z[[1]]@u)) ## amending slack pertinent to f0 (scaled)
            } else {
                z <- c(new("NLFV", u = Matrix(-x[ne] * W[[1]]@W$dnl[1]), dims = 1L), z) ## adding slack (scaled) pertinent to f0
            }
            x1 <- crossprod(GGList[[1]][1, -ne], x0) + W[[1]]@W$dnl[1]^2 * x[ne] - a
            x <- c(x0, x1)
            ## Computing value for 't' and associated slack
            
            return(list(x = x, y = y, z = z))
        }
        return(new("KKTSLV", f = kktslv, items = list(K = K, W = W, GGList = GGList, ne = ne, mnl = cpd@mnl)))
    }
    return(kktSLV)
}
