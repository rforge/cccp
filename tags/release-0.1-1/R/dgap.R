##
## duality gap s' * z for PDV-objects
setMethod("dgap", signature = "PDV", function(pdv){
    k <- length(pdv@s)
    gap <- 0
    for(j in 1:k){
        gap <- gap + udot(pdv@s[[j]], pdv@z[[j]]) 
    }
    gap
})
