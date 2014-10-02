##
## Coercion of DEFCP to DEFNL objects
setAs("DEFCP", "DEFNL", function(from){
    ##
    ## transforming CP to its epigraph form
    ## additional variable 't' is left-added to objective/constraints
    ##
    ## linear objective
    q <- c(1, rep(0, from@n))
    n <- from@n + 1L
    if(from@k > 0){
        cList <- lapply(from@cList, function(cc){
            cc@G <- cbind2(0, cc@G)
            cc
        })
    } else {
        cList <- list()
    }
    ## objective as nonlinear constraint
    f0e <- function(x){
        from@f0(x[-1]) - x[1]
    }
    ## including objective to other nonlinear constraints, if applicable
    if(from@mnl > 0){
        nlfeList <- lapply(from@nlfList, function(nlf){
            fnle <- function(x){
                nlf(x[-1])
            }
            fnle
        })
        nlfeList <- c(f0e, nlfeList)
        mnl <- from@mnl + 1L
        cList[[1]]@G <- rbind2(0, cList[[1]]@G)
        cList[[1]]@dims <- cList[[1]]@dims + 1L
        cList[[1]]@h@u <- rbind2(0, cList[[1]]@h@u)
        cList[[1]]@h@dims <- cList[[1]]@h@dims + 1L
    } else {
        nlfeList <- list(f0e)
        mnl <- 1L
        h <- new("NLFV", u = Matrix(0, nrow = mnl), dims = mnl)
        nlfc <- new("NLFC",
                    G = Matrix(0, nrow = mnl, ncol = n),
                    h = h,
                    dims = mnl,
                    vclass = "NLFV")
        cList <- c(nlfc, cList)
    }
    k <- length(cList)
    ## Amending LHS of equality constraints, if applicable
    if(dim(from@A)[1] == 0){
        A <- Matrix(0, nrow = 0, ncol = n)
    } else {
        A <- cbind2(0, from@A)
    }

    new("DEFNL",
        x0 = c(0.0, from@x0),
        q = q,
        nlfList = nlfeList,
        cList = cList,
        A = A,
        b = from@b,
        k = k,
        n = n,
        mnl = mnl,
        H = Matrix(0, nrow = n, ncol = n),
        ctrl = from@ctrl)
    })
