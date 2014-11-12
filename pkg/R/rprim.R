##
## Primal residuals for CPs with equality constraints
setMethod("rprim", signature = c("PDV", "CPD"), function(pdv, cpd){
    drop(cpd@b - cpd@A %*% pdv@x)
})
