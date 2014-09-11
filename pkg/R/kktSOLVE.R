##
## Solving KKT-system by LU factorization
kktSOLVE <- function(cpd){
    kktSLV <- function(W, cpd){
        ldK <- sum(dim(cpd@A))
        K <- Matrix(0, nrow = ldK, ncol = ldK)
        K <- as(K, "dsyMatrix")
        n <- cpd@n
        idx <- 1:cpd@k
        if(class(cpd) == "DEFQP") K[1:n, 1:n] <- cpd@P
        if(class(cpd) == "DEFNL") K[1:n, 1:n] <- cpd@H
        K[-(1:n), 1:n] <- cpd@A
        K[1:n, -(1:n)] <- t(cpd@A)
        for(j in idx){
            WitG <- usnt(cpd@cList[[j]]@G, W[[j]], trans = TRUE, inv = TRUE)
            WiWitG <- usnt(WitG@u, W[[j]], trans = FALSE, inv = TRUE)
            GtWiWitG <- crossprod(cpd@cList[[j]]@G, WiWitG@u) 
            K[1:n, 1:n] <- K[1:n, 1:n] + GtWiWitG
        }            
        kktslv <- function(x, y, z){
            if(!is.list(z)) z <- list(z)
            GtWiWitz <- Reduce("+", lapply(idx, function(j){
                Witz <- usnt(z[[j]]@u, W[[j]], inv = TRUE, trans = TRUE)
                WiWitz <- usnt(Witz@u, W[[j]], inv = TRUE, trans = FALSE)
                GtWiWitz <- crossprod(cpd@cList[[j]]@G, WiWitz@u)
                GtWiWitz
            }))
            RHS <- c(x + drop(GtWiWitz), drop(y))
            ans <- solve(K, RHS)
            x <- ans[1:n, 1]
            y <- ans[-(1:n), 1]
            z <- lapply(idx, function(j){
                uz <- new(cpd@cList[[j]]@vclass, u = cpd@cList[[j]]@G %*% x - z[[j]]@u, dims = cpd@cList[[j]]@dims)
                usnt(uz@u, W[[j]], inv = TRUE, trans = TRUE)
            })
            return(list(x = x, y = y, z = z))
        }
        return(new("KKTSLV", f = kktslv, items = list(K = K, W = W)))
    }
    return(kktSLV)
}
