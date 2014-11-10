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
        h <- new("NLFV", u = matrix(0, nrow = mnl), dims = mnl)
        nlfc <- new("NLFC", G = matrix(0, nrow = mnl, ncol = n), h = h, dims = mnl, vclass = "NLFV")
        cList <- c(nlfc, cList)
    }
    ## Equality constraints
    if(is.null(A)){
        A <- matrix(0, nrow = 0, ncol = n)
    } 
    if(is.null(dim(A))){
        A <- matrix(A, nrow = 1)
    } 
    if(is.null(b)){
        b <- matrix(0, nrow = 0, ncol = 1)
    }
    if(is.null(dim(b))){
        b <- matrix(b, ncol = 1)
    } 
    ## checking whether x0 satisfies equality constraints
    if(!is.null(A)){
        eq <- identical(as.vector(A %*% x0 - b), rep(0, nrow(A)))
        if(!eq){
            stop("Initial point 'x0' does not satisfy equality constraints.\n")
        }
    }
    cpdef <- new("DEFCP",
                 q = c(rep(0, n), 1),
                 x0 = x0,
                 f0 = f0,
                 nlfList = nlfList,
                 cList = cList,
                 A = A,
                 b = b,
                 k = length(cList),
                 n = n,
                 mnl = mnl,
                 H = matrix(0, nrow = n, ncol = n),
                 ctrl = optctrl)
    cpsol <- cps(cpdef)
    return(cpsol)   
}
