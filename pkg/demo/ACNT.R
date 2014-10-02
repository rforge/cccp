##
## Demo for solving an analytic centering problem with equality constraints
##
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
acnt <- cpn(x0 = x0, f0 = f0, A = A, b = b)
acnt
getx(acnt)
