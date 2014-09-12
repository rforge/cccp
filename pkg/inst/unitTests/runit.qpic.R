##
## Unit testing of QPIC from demo
test.QPIC <- function(){
    P <- 2 * Matrix(c(2, .5, .5, 1), nrow = 2, ncol = 2)
    q <- c(1.0, 1.0)
    G <- -diag(2)
    h <- rep(0, 2)
    nno1 <- nnoc(G = G, h = h)
    A <- Matrix(c(1.0, 1.0), nrow = 1, ncol = 2)
    b <- 1.0
    ans <- cccp(P = P, q = q, A = A, b = b, cList = list(nno1))
    checkTrue(ans@status == "optimal")
    checkTrue(ans@certp <= 1e-7)
    checkTrue(ans@certd <= 1e-7)
    return()
}
