##
## Methods for value of primal objective at point 'x' of 'DEFLP'
setMethod("pobj", signature = c("numeric", "DEFLP"), function(pdv, cpd){
    drop(crossprod(pdv, cpd@q))
})
setMethod("pobj", signature = c("PDV", "DEFLP"), function(pdv, cpd){
    drop(crossprod(pdv@x, cpd@q))
})
##
## Methods for value of primal objective at point 'x' of 'DEFQP'
setMethod("pobj", signature = c("numeric", "DEFQP"), function(pdv, cpd){
    drop(0.5 * pdv %*% cpd@P %*% pdv + crossprod(cpd@q, pdv))
})
setMethod("pobj", signature = c("PDV", "DEFQP"), function(pdv, cpd){
    drop(0.5 * pdv@x %*% cpd@P %*% pdv@x + crossprod(cpd@q, pdv@x))
})
##
## Methods for value of primal objective at point 'x' of 'DEFNL'
setMethod("pobj", signature = c("numeric", "DEFNL"), function(pdv, cpd){
    drop(crossprod(pdv, cpd@q))
})
setMethod("pobj", signature = c("PDV", "DEFNL"), function(pdv, cpd){
    drop(crossprod(pdv@x, cpd@q))
})
