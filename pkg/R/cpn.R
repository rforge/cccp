##
## Function for solving a convex objective with general convex and cone constraints
cpn <- function(x0, f0, nlfList = list(), cList = list(), A = NULL, b = NULL, optctrl = ctrl()){
    x0 <- as.vector(x0)
    n <- length(x0)
    k <- length(cList)
    ## Checking whether x0 is in the domain of nonlinear objective
    f0Dom <- is.nan(f0(x0))
    if(f0Dom){
        stop("Initial point 'x0' is not in the domain of nonlinear objective 'f0'.\n")
    }
    if(length(cList) > 0){
        coneclasses <- unlist(lapply(cList, class))
        if(!all(coneclasses %in% c("NNOC", "SOCC", "PSDC"))){
            stop("List elements of cone constraints must be of either class 'NNOC', or 'SOCC', or 'PSDC'.\n")
        }
    }
    mnl <- length(nlfList)
    if(mnl > 0){
        ## Checking whether x0 is in the domain of nonlinear constraint(s)
        fDom <- unlist(lapply(nlfList, function(fcc) fcc(x0)))
        idxnan <- which(is.nan(fDom))
        if(any(idxnan)){
            stop(paste("Initial point 'x0' is not in the domain of nonlinear constraint(s): ", idxnan, ".\n", sep = ""))
        }
        ## creating initial NLFC object and adding to cList as first object
        h <- new("NLFV", u = Matrix(0, nrow = mnl), dims = mnl)
        nlfc <- new("NLFC", G = Matrix(0, nrow = mnl, ncol = n), h = h, dims = mnl, vclass = "NLFV")
        cList <- c(nlfc, cList)
    }
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
    ## checking whether x0 satisfies equality constraints
    if(!is.null(A)){
        eq <- identical(as.vector(A %*% x0 - b), rep(0, nrow(A)))
        if(!eq){
            stop("Initial point 'x0' does not satisfy equality constraints.\n")
        }
    }
    cpdef <- new("DEFCP",
                 x0 = x0,
                 f0 = f0,
                 nlfList = nlfList,
                 cList = cList,
                 A = A,
                 b = b,
                 k = length(cList),
                 n = n,
                 mnl = mnl,
                 H = Matrix(0, nrow = n, ncol = n),
                 ctrl = optctrl)
    cpnl <- as(cpdef, "DEFNL")
    cpsol <- cps(cpnl)
    cpsol@x <- cpsol@x[-1]
    if(length(nlfList) > 0){
        cpsol@s[[1]]@u <- Matrix(cpsol@s[[1]]@u[-1, 1])
        cpsol@s[[1]]@dims <- cpsol@s[[1]]@dims - 1L
        cpsol@z[[1]]@u <- Matrix(cpsol@z[[1]]@u[-1, 1])
        cpsol@z[[1]]@dims <- cpsol@z[[1]]@dims - 1L
    } else {
        cpsol@s <- cpsol@s[-1]
        cpsol@z <- cpsol@z[-1]
    }
    return(cpsol)   
}
