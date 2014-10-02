##
## Unit testing of ACNT from demo
test.ACNT <- function(){
    ## Creating objective
    f0 <- function(x){
        -sum(log(x))
    }
    ## equality constraint
    A <- Matrix(c(1, 1, 2), nrow = 1)
    b <- Matrix(1, nrow = 1)
    ## initial (feasible!) point
    x0 = c(0.25, 0.25, 0.25)
    ## solving problem
    ans <- cpn(x0 = x0, f0 = f0, A = A, b = b)
    checkTrue(ans@status == "optimal")
    checkTrue(ans@certp <= 1e-7)
    checkTrue(ans@certd <= 1e-7)
    return()
}
