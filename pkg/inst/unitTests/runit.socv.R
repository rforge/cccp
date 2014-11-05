##
## Unit testing of operators and NT scaling
test.SOCV <- function(){
    set.seed(12345)
    ## Creating SOC-points
    p <- 10L
    s <- new("SOCV", u = matrix(c(p + 1, runif(p - 1))), dims = p)
    z <- new("SOCV", u = matrix(c(p + 1, runif(p - 1))), dims = p)
    onep <- uone(s)
    checkEquals(sum(onep@u * onep@u), 1)
    term1 <- sum(uprd(s, z)@u * onep@u)
    term2 <- sum(s@u * z@u)
    checkEquals(term1, term2)
    checkEquals(uprd(s, onep), s)
    checkEquals(uprd(onep, s), s, check.attributes = FALSE)
    checkEquals(uprd(s, s), uprd(s))
    W <- ntsc(s, z)
    checkEqualsNumeric(usnt(s@u, W, inv = TRUE)@u, usnt(z@u, W, inv = FALSE)@u)
    return()
}
