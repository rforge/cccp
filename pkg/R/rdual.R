##
## Dual residuals for LPs with equality and inequality constraints
setMethod("rdual", signature = c("PDV", "DEFLP"), function(pdv, cpd){
    Gz <- 0
    for(i in 1:cpd@k) Gz <- Gz + crossprod(cpd@cList[[i]]@G, pdv@z[[i]]@u)
    Ay <- crossprod(cpd@A, pdv@y)
    drop(Ay + Gz + cpd@q)
})
## Dual residuals for QPs with equality and inequality constraints
setMethod("rdual", signature = c("PDV", "DEFQP"), function(pdv, cpd){
    Gz <- 0
    for(i in 1:cpd@k) Gz <- Gz + crossprod(cpd@cList[[i]]@G, pdv@z[[i]]@u)
    Ay <- crossprod(cpd@A, pdv@y)
    drop(cpd@P %*% pdv@x + cpd@q + Gz + Ay)
})
##
## Dual residuals for LPs with nonlinear/cone inequality and equality constraints
setMethod("rdual", signature = c("PDV", "DEFNL"), function(pdv, cpd){
    Gz <- 0
    for(i in 1:cpd@k) Gz <- Gz + crossprod(cpd@cList[[i]]@G, pdv@z[[i]]@u)
    Ay <- crossprod(cpd@A, pdv@y)
    drop(Ay + Gz + cpd@q)
})
##
## Dual residuals for CPs with nonlinear/cone inequality and equality constraints
setMethod("rdual", signature = c("PDV", "DEFCP"), function(pdv, cpd){
    Gz <- 0
    for(i in 1:cpd@k) Gz <- Gz + crossprod(cpd@cList[[i]]@G, pdv@z[[i]]@u)
    Ay <- crossprod(cpd@A, pdv@y)
    drop(Ay + Gz + cpd@q)
})
