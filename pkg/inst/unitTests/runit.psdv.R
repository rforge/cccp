##
## Unit testing of operators and NT scaling
test.PSDV <- function(){
    set.seed(12345)
    ## Creating PSD-points
    n <- 4L
    s <- crossprod(matrix(rnorm(n^2), ncol = n, nrow = n))
    dim(s) <- c(n^2, 1)
    s <- new("PSDV", u = s, dims = n)
    z <- crossprod(matrix(rnorm(n^2), ncol = n, nrow = n))
    dim(z) <- c(n^2, 1)
    z <- new("PSDV", u = z, dims = n)
    onep <- uone(s)
    checkEquals(sum(onep@u * onep@u), s@dims)
    term1 <- sum(uprd(s, z)@u * onep@u)
    term2 <- sum(s@u * z@u)
    checkEquals(term1, term2)
    checkEqualsNumeric(uprd(s, onep)@u, s@u)
    checkEqualsNumeric(uprd(onep, s)@u, s@u)
    checkEqualsNumeric(uprd(s, s)@u, uprd(s)@u)
    W <- ntsc(s, z)
    ws <- usnt(s@u, W, trans = TRUE, inv = TRUE)@u
    wz <- usnt(z@u, W, trans = FALSE, inv = FALSE)@u
    checkEqualsNumeric(ws, wz)
    return()
}
