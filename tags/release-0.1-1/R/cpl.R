##
## Function for solving a linear objective with general convex and cone constraints
cpl <- function(x0, q, nlfList = list(), nlgList = list(), nlhList = list(), cList = list(), A = NULL, b = NULL, optctrl = ctrl()){
    if(is.matrix(x0)){
        warning("Matrix provided for x0, extracting first column for argument 'x0'.\n")
        x0 <- x0[, 1]
    }
    if(is.matrix(q)){
        warning("Matrix provided for q, extracting first column for argument 'q'.\n")
        q <- q[, 1]
    }
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
    if(!all(unlist(lapply(nlfList, function(f) class(f) == "function")))){
        stop("Not all list elements in 'nlfList' are functions.\n")
    }
    if(!all(unlist(lapply(nlgList, function(f) class(f) == "function")))){
        stop("Not all list elements in 'nlgList' are functions.\n")
    }
    if(!all(unlist(lapply(nlhList, function(f) class(f) == "function")))){
        stop("Not all list elements in 'nlhList' are functions.\n")
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
    if(length(nlfList) != length(nlgList)){
        stop("Length of lists for nonlinear functions and gradient functions do differ.\n")
    }
    if(length(nlfList) != length(nlhList)){
        stop("Length of lists for nonlinear functions and Hessian functions do differ.\n")
    }
    ## creating initial NLFC object and adding to cList as first object
    h <- new("NLFV", u = matrix(0, nrow = mnl), dims = mnl)
    nlfc <- new("NLFC", G = matrix(0, nrow = mnl, ncol = n), h = h, dims = mnl, vclass = "NLFV")
    cList <- c(nlfc, cList)
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
    cpdef <- new("DEFNL",
                 x0 = x0,
                 q = q,
                 nlfList = nlfList,
                 nlgList = nlgList,
                 nlhList = nlhList,
                 cList = cList,
                 A = A,
                 b = b,
                 k = k + 1L,
                 n = n,
                 mnl = mnl,
                 H = matrix(0, nrow = n, ncol = n),
                 ctrl = optctrl)
    cpsol <- cps(cpdef)
    return(cpsol)   
}
