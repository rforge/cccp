##
## Unit testing of LPIC from demo
test.LPIC <- function(){
    q <- c(-4, -5)
    G <- matrix(c(2, 1, -1, 0,
                  1, 2, 0, -1),
                nrow = 4, ncol = 2)
    h <- c(3, 3, 0, 0)
    nno1 <- nnoc(G = G, h = h)
    ans <- cccp(q = q, cList = list(nno1), optctrl = ctrl(trace = FALSE))
    checkTrue(ans@status == "optimal")
    checkTrue(ans@certp <= 1e-7)
    checkTrue(ans@certd <= 1e-7)
    return()
}
