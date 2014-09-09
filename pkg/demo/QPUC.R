##
## Demo for solving an unconstrained QP
##
## Creating QP
## Objective
n <- 4L
M <- Matrix(rnorm(n^2), nrow = n, ncol = n)
P <- crossprod(M)
q <- rnorm(n)
## Using main function of package
qpucsol <- cccp(P = P, q = q)
qpucsol
getx(qpucsol)
##
## Equivalently:
qpuccpd <- cpd(P = P, q = q)
cps(qpuccpd)
