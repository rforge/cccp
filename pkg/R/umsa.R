##
## Adjusting LNL-variable by maximum step-size
setMethod("umsa", signature = "LNLV", function(u, alpha, init = TRUE){
    if(init){
        u@u <- u@u + (1 + alpha)
    } else {
        u@u <- uone(u)@u + alpha * u@u
    }
    u
})
##
## Adjusting SOC-variable by maximum step-size
setMethod("umsa", signature = "SOCV", function(u, alpha, init = TRUE){
    if(init){
        u@u[1] <- u@u[1] + (1 + alpha)
    } else {
        u@u <- uone(u)@u + alpha * u@u
    }
    u
})
##
## Adjusting PSD-variable by maximum step-size
setMethod("umsa", signature = "PSDV", function(u, alpha, init = TRUE, sigma = NULL, lambda = NULL){
    if(init){
        u@u <- u@u + uone(u)@u * (1 + alpha)
    } else {
        sig <- 1 + alpha * sigma
        dim(lambda@u) <- c(lambda@dims, lambda@dims)
        sig <- sig / diag(x = lambda@u)
        dim(u@u) <- c(u@dims, u@dims)
        for(j in 1:u@dims){
            u@u[, j] <- u@u[, j] * sqrt(sig[j])
        }
        dim(u@u) <- c(u@dims^2, 1)
    }
    u
})
