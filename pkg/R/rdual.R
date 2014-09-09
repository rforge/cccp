##
## Dual residuals for LPs with equality and inequality constraints
setMethod("rdual", signature = c("PDV", "DEFLP"), function(pdv, cpd){
    Gz <- 0
    for(i in 1:cpd@k) Gz <- Gz + crossprod(cpd@conecon[[i]]@G, pdv@z[[i]]@u)
    Ay <- crossprod(cpd@A, pdv@y)
    drop(Ay + Gz + cpd@q)
})
## Dual residuals for QPs with equality and inequality constraints
setMethod("rdual", signature = c("PDV", "DEFQP"), function(pdv, cpd){
    Gz <- 0
    for(i in 1:cpd@k) Gz <- Gz + crossprod(cpd@conecon[[i]]@G, pdv@z[[i]]@u)
    Ay <- crossprod(cpd@A, pdv@y)
    drop(cpd@P %*% pdv@x + cpd@q + Gz + Ay)
})
