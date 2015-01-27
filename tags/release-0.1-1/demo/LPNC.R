##
## Demo for solving a Linear Program with nonlinear inequality constraints
## (Example taken from cvxopt's userguide)
##
require(numDeriv)
## Creating LP
q <- c(-1, -1)
G <- -diag(2)
h <- c(0, 0)
nno1 <- nnoc(G = G, h = h)
f1 <- function(x) 1 - x[1]^2 - x[2]^2
f2 <- function(x) -2 + x[1]^2 + x[2]^2
g1 <- function(x, func = f1) grad(func = func, x = x)
g2 <- function(x, func = f2) grad(func = func, x = x)
h1 <- function(x, func = f1) hessian(func = func, x = x)
h2 <- function(x, func = f2) hessian(func = func, x = x)
ans <- cpl(x0 = c(2, 3), q = q,
           nlfList = list(f1, f2),
           nlgList = list(g1, g2),
           nlhList = list(h1, h2),
           cList = list(nno1))
getx(ans)
