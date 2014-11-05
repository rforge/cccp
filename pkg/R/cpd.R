## Function for creating an object of either S4-class 'DEFLP' or 'DEFQP'
cpd <- function(P, q, A = NULL, b = NULL, cList = list(), optctrl = ctrl()){
    k <- length(cList)
    if(is.null(P)){
        n <- length(q)
    } else {
        n <- ncol(P)
    }
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
    if(is.null(P)){
        if(k < 1){
            if(optctrl@trace){
                cat("LP in standard form: Adding non-negativity constraint(s).\n")
            }
            G <- -diag(n)
            h <- rep(0, n)
            cList <-  list(nnoc(G = G, h = h))
            k <- 1L
        }
        ans <- new("DEFLP",
                   q = q,
                   A = A,
                   b = b,
                   n = n,
                   k = k,
                   cList = cList,
                   ctrl = optctrl)
        return(ans)
    } else {
        ans <- new("DEFQP",
                   P = P,
                   q = q,
                   A = A,
                   b = b,
                   n = n,
                   k = k,
                   cList = cList,
                   ctrl = optctrl)
        return(ans)
    }
}
