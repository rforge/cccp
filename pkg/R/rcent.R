##
## Centrality residuals for LPs with inequality constraints
setMethod("rcent", signature = c("PDV", "DEFLP"), function(pdv, cpd){
    lapply(1:cpd@k, function(j){
        ans <- pdv@s[[j]]
        ans@u <- pdv@s[[j]]@u + cpd@cList[[j]]@G %*% pdv@x - cpd@cList[[j]]@h@u
        ans
    })
})
##
## Centrality residuals for LPs with inequality constraints
setMethod("rcent", signature = c("PDV", "DEFQP"), function(pdv, cpd){
    lapply(1:cpd@k, function(j){
        ans <- pdv@s[[j]]
        ans@u <- pdv@s[[j]]@u + cpd@cList[[j]]@G %*% pdv@x - cpd@cList[[j]]@h@u
        ans
    })
})
##
## Centrality residuals for LPs with nonlinear/cone inequality and equality constraints
setMethod("rcent", signature = c("PDV", "DEFNL"), function(pdv, cpd){
    ans1 <- pdv@s[[1]]
    ans1@u <- pdv@s[[1]]@u + cpd@cList[[1]]@h@u
    ans2 <- NULL
    if(cpd@k > 1){
        ans2 <- lapply(2:cpd@k, function(j){
            ans <- pdv@s[[j]]
            ans@u <- pdv@s[[j]]@u + cpd@cList[[j]]@G %*% pdv@x - cpd@cList[[j]]@h@u
            ans
        })
    }
    c(ans1, ans2)
})
##
## Centrality residuals for LPs with nonlinear/cone inequality and equality constraints
setMethod("rcent", signature = c("PDV", "DEFCP"), function(pdv, cpd){
    ans1 <- pdv@s[[1]]
    ans1@u <- pdv@s[[1]]@u + cpd@cList[[1]]@h@u
    ans2 <- NULL
    if(cpd@k > 1){
        ans2 <- lapply(2:cpd@k, function(j){
            ans <- pdv@s[[j]]
            ans@u <- pdv@s[[j]]@u + cpd@cList[[j]]@G %*% pdv@x - cpd@cList[[j]]@h@u
            ans
        })
    }
    c(ans1, ans2)
})
