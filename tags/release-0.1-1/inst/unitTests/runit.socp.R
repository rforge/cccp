##
## Unit testing of SOCP from demo
test.SOCP <- function(){
    q <- c(-2, 1, 5)
    ## First SOC
    F1 <- matrix(c(-13, 3, 5,
                   -12, 12, -6),
                 nrow = 2, ncol = 3, byrow = TRUE)
    g1 <- c(-3, -2)
    d1 <- c(-12, -6, 5)
    f1 <- -12
    soc1 <- socc(F = F1, g = g1, d = d1, f = f1)
    ## Second SOC
    F2 <- matrix(c(-3, 6, 2,
                   1, 9, 2,
                   -1, -19, 3),
                 nrow = 3, ncol = 3, byrow = TRUE)
    g2 <- c(0, 3, -42)
    d2 <- c(-3, 6, -10)
    f2 <- 27
    soc2 <- socc(F = F2, g = g2, d = d2, f = f2)
    ans <- cccp(q = q, cList = list(soc1, soc2), optctrl = ctrl(trace = FALSE))
    checkTrue(ans@status == "optimal")
    checkTrue(ans@certp <= 1e-7)
    checkTrue(ans@certd <= 1e-7)
    return()
}
