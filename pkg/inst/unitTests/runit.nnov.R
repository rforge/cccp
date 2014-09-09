##
## Unit testing of operators and NT scaling
test.NNOV <- function(){
    set.seed(12345)
    ## Creating NNO-points
    n <- 10L
    s <- new("NNOV", u = Matrix(runif(n)), dims = n)
    z <- new("NNOV", u = Matrix(runif(n)), dims = n)
    onep <- uone(s)
    checkEquals(sum(onep@u * onep@u), s@dims)
    term1 <- sum(uprd(s, z)@u * onep@u)
    term2 <- sum(s@u * z@u)
    checkEquals(term1, term2)  
    checkEquals(uprd(s, onep), s)
    checkEquals(uprd(onep, s), s)
    checkEquals(uprd(s, s), uprd(s))
    W <- ntsc(s, z)
    checkEqualsNumeric(usnt(s@u, W, inv = TRUE)@u, usnt(z@u, W, inv = FALSE)@u)
    return()
}
