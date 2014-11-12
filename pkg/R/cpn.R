##
## Function for solving a convex objective with general convex and cone constraints
cpn <- function(x0, f0, g0, h0, nlfList = list(), nlgList = list(), nlhList = list(),
                cList = list(), A = NULL, b = NULL, optctrl = ctrl()){
    if(is.matrix(x0)){
        warning("Matrix provided for x0, extracting first column for argument 'x0'.\n")
        x0 <- x0[, 1]
    }
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
        ## Checking length of lists and elements
        if(mnl != length(nlgList)){
            stop("Length of lists for nonlinear functions and gradient functions do differ.\n")
        }
        if(mnl != length(nlhList)){
            stop("Length of lists for nonlinear functions and Hessian functions do differ.\n")
        }
        if(!all(unlist(lapply(nlfList, function(f) class(f) == "function")))){
            stop("Not all list elements in 'nlfList' are functions.\n")
        }
        if(!all(unlist(lapply(nlgList, function(f) class(f) == "function")))){
            stop("Not all list elements in 'nlgList' are functions.\n")
        }
        if(!all(unlist(lapply(nlhList, function(f) class(f) == "function")))){
            stop("Not all list elements in 'nlhList' are functions.\n")
        }
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
                 g0 = g0,
                 h0 = h0,
                 nlfList = nlfList,
                 nlgList = nlgList,
                 nlhList = nlhList,
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
