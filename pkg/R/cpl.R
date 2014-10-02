##
## Function for solving a linear objective with general convex and cone constraints
cpl <- function(x0, q, nlfList = list(), cList = list(), A = NULL, b = NULL, optctrl = ctrl()){
    x0 <- as.vector(x0)
    q <- as.vector(q)
    n <- length(q)
    k <- length(cList)
    if(!identical(length(x0), length(q))){
        stop("Length of initial point 'x' is not equal to the dimension of the objective.\n")
    }
    mnl <- length(nlfList)
    if(mnl < 1){
        warning("Empty list for non-linear convex constraints provided.\nSolving as LP with cone and equality constraints, if applicable.\n")
        ans <- cccp(P = NULL, q = q, cList = cList, A = A, b = b, optctrl = optctrl)
        return(ans)
    }
    if(length(cList) > 0){
        coneclasses <- unlist(lapply(cList, class))
        if(!all(coneclasses %in% c("NNOC", "SOCC", "PSDC"))){
            stop("List elements of cone constraints must be of either class 'NNOC', or 'SOCC', or 'PSDC'.\n")
        }
    }
    fDom <- unlist(lapply(nlfList, function(fcc) fcc(x0)))
    idxnan <- which(is.nan(fDom))
    if(any(idxnan)){
        stop(paste("Initial point 'x0' is not in the domain of nonlinear convex constraint(s): ", idxnan, ".\n", sep = ""))
    }
    ## creating initial NLFC object and adding to cList as first object
    h <- new("NLFV", u = Matrix(0, nrow = mnl), dims = mnl)
    nlfc <- new("NLFC", G = Matrix(0, nrow = mnl, ncol = n), h = h, dims = mnl, vclass = "NLFV")
    cList <- c(nlfc, cList)
    ## Equality constraints
    if(is.null(A)){
        A <- Matrix(0, nrow = 0, ncol = n)
    } 
    if(is.null(dim(A))){
        A <- Matrix(A, nrow = 1)
    } else {
        if(!any(extends(class(A)) %in% "Matrix"))  A <- as(A, "Matrix")
    }
    if(is.null(b)){
        b <- Matrix(0, nrow = 0, ncol = 1)
    }
    if(is.null(dim(b))){
        b <- Matrix(b, ncol = 1)
    } else {
        if(!any(extends(class(A)) %in% "Matrix"))  b <- as(b, "Matrix")
    }
    cpdef <- new("DEFNL",
                 x0 = x0,
                 q = q,
                 nlfList = nlfList,
                 cList = cList,
                 A = A,
                 b = b,
                 k = k + 1L,
                 n = n,
                 mnl = mnl,
                 H = Matrix(0, nrow = n, ncol = n),
                 ctrl = optctrl)
    cpsol <- cps(cpdef)
    return(cpsol)   
}
