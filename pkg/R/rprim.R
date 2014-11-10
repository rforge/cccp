##
## Primal residuals for LPs with equality constraints
setMethod("rprim", signature = c("PDV", "DEFLP"), function(pdv, cpd){
    drop(cpd@b - cpd@A %*% pdv@x)
})
##
## Primal residuals for QPs with equality constraints
setMethod("rprim", signature = c("PDV", "DEFQP"), function(pdv, cpd){
    drop(cpd@A %*% pdv@x - cpd@b)
})
##
## Primal residuals for LPs with nonlinear/cone inequality and equality constraints
setMethod("rprim", signature = c("PDV", "DEFNL"), function(pdv, cpd){
    drop(cpd@b - cpd@A %*% pdv@x)
})
##
## Primal residuals for LPs with nonlinear/cone inequality and equality constraints
setMethod("rprim", signature = c("PDV", "DEFCP"), function(pdv, cpd){
    drop(cpd@b - cpd@A %*% pdv@x)
})
