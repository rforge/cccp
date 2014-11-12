##
## Demo for solving a nonlinear convex program
## (Example taken from cvxopt's userguide)
##
require(numDeriv)
## Creating Problem
## objective, gradient and Hessian
f0 <- function(x) -log(1 - x[1]^2) - log(1 - x[2]^2) - log(1 - x[3]^2)
g0 <- function(x, func = f0) grad(func = func, x = x)
h0 <- function(x, func = f0) hessian(func = func, x = x)
## SOC
F1 <- diag(3)
g1 <- rep(0, 3)
d1 <- rep(0, 3)
f1 <- 1
soc1 <- socc(F = F1, g = g1, d = d1, f = f1)
## PSD
F1 <- matrix(c(-21, -11, 0, -11, 10, 8, 0, 8, 5), nrow = 3, ncol = 3)
F2 <- matrix(c(0, 10, 16, 10, -10, -10, 16, -10, 3), nrow = 3, ncol = 3)
F3 <- matrix(c(-5, 2, -17, 2, -6, 8, -17, 8, 6), nrow = 3, ncol = 3)
F0 <- matrix(c(20, 10, 40, 10, 80, 10, 40, 10, 15), nrow = 3, ncol = 3)
psd1 <- psdc(Flist = list(F1, F2, F3), F0 = F0)
## Solving problem
ans <- cpn(x0 = c(0.5, 0.5, 0.5), f0 = f0, g0 = g0, h0 = h0,
           cList = list(soc1, psd1))
ans
getx(ans)
