##
## Demo for solving a Linear Program with nonlinear inequality constraints
## (Example taken from cvxopt's userguide)
##
## Creating LP
q <- c(-1, -1)
G <- -diag(2)
h <- c(0, 0)
nno1 <- nnoc(G = G, h = h)
f1 <- function(x){
    1 - x[1]^2 - x[2]^2
}
f2 <- function(x){
    -2 + x[1]^2 + x[2]^2
}
ans <- cpl(x0 = c(2, 3), q = q, nlfList = list(f1, f2), cList = list(nno1))
getx(ans)
