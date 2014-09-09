##
## Function for minimizing a L1-norm ||Pu - q||_1
## This is a wrapper function for the LP-cps method
l1 <- function(P, q = NULL, optctrl = ctrl()){
    m <- nrow(P)
    n <- ncol(P)
    if(is.null(q)){
        q <- rep(0, m)
    } else {
        q <- as.vector(q)
    }
    ## Creating NNO-constraint
    target <- c(rep(0, n), rep(1, m))
    G <- Matrix(0, nrow = m + m, ncol = n + m)
    G[1:m, 1:n] <- P
    D <- Diagonal(m)
    G[1:m, -c(1:n)] <- -D
    G[-c(1:m), 1:n] <- -P
    G[-c(1:m), -c(1:n)] <- -D
    h <- Matrix(c(q, -q), nrow = m + m, ncol = 1)
    nno1 <- nnoc(G = G, h = h)
    ## Defining LP and solving
    l1cpd <- cpd(P = NULL, q = target, conecon = list(nno1), optctrl = optctrl)
    l1cps <- cps(l1cpd)
    l1cps
}
