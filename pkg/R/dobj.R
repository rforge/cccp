##
## Methods for value of dual objective at point 'x' 'y' of 'DEFLP'
setMethod("dobj", signature = c("PDV", "DEFLP"), function(pdv, cpd){
    term1 <- drop(crossprod(cpd@b, pdv@y))
    term2 <- sum(unlist(lapply(1:cpd@k, function(j) udot(pdv@z[[j]], cpd@cList[[j]]@h))))
    -term1 - term2
})
##
## Methods for value of dual objective at point 'x' 'y' of 'DEFQP'
setMethod("dobj", signature = c("PDV", "DEFQP"), function(pdv, cpd){
    term1 <- pobj(pdv, cpd)
    term2 <- 0
    if(cpd@k > 0){
        term2 <- Reduce("+", lapply(1:cpd@k, function(j){
            crossprod(pdv@z[[j]]@u, cpd@cList[[j]]@G %*% pdv@x - cpd@cList[[j]]@h@u)
        }))
    }
    term3 <- crossprod(pdv@y, cpd@A %*% pdv@x - cpd@b)
    drop(term1 + term2 + term3)
})
##
## Methods for value of dual objective at point 'x' 'y' of 'DEFNL'
setMethod("dobj", signature = c("PDV", "DEFNL"), function(pdv, cpd){
    term1 <- pobj(pdv, cpd)
    term2 <- 0
    if(cpd@k > 0){
        term2 <- Reduce("+", lapply(2:(cpd@k + 1), function(j){
            crossprod(pdv@z[[j]]@u, cpd@cList[[j]]@G %*% pdv@x - cpd@cList[[j]]@h@u)
        }))
    }
    term2 <- term2 + crossprod(pdv@z[[1]]@u, cpd@cList[[1]]@h@u)
    term3 <- crossprod(pdv@y, cpd@A %*% pdv@x - cpd@b)
    drop(term1 + term2 + term3)
})
##
## Methods for value of dual objective at point 'x' 'y' of 'DEFCP'
setMethod("dobj", signature = c("PDV", "DEFCP"), function(pdv, cpd){
    term1 <- pobj(pdv, cpd)
    term2 <- 0
    if(cpd@k > 1){
        term2 <- Reduce("+", lapply(2:cpd@k, function(j){
            crossprod(pdv@z[[j]]@u, cpd@cList[[j]]@G %*% pdv@x - cpd@cList[[j]]@h@u)
        }))
    }
    term2 <- term2 + crossprod(pdv@z[[1]]@u, cpd@cList[[1]]@h@u)
    term3 <- crossprod(pdv@y, cpd@A %*% pdv@x - cpd@b)
    drop(term1 + term2 + term3)
})
