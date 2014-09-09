##
## Dual residuals for LPs with inequality constraints
setMethod("rcent", signature = c("PDV", "CPD"), function(pdv, cpd){
    lapply(1:cpd@k, function(j){
        ans <- pdv@s[[j]]
        ans@u <- pdv@s[[j]]@u + cpd@conecon[[j]]@G %*% pdv@x - cpd@conecon[[j]]@h@u
        ans
    })
})
