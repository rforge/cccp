##
## Validity function for S$-objects of class 'CTRL'
validCTRL <- function(object){
    if(!is.integer(object@maxiters)){
        return("\nThe count of maximal iterations must be an integer.\n")
    }
    if(object@maxiters < 1){
        return("\nThe count of maximal iterations must be positive and greater or equal to one.\n")
    }
    if(!is.null(dim(object@abstol)) | length(object@abstol) > 1){
        return("\nThe absolute tolerance for convergence must be a real scalar.\n")
    }
    if(!is.null(dim(object@reltol)) | length(object@reltol) > 1){
        return("\nThe relative tolerance for convergence must be a real scalar.\n")
    }
    if(!is.null(dim(object@feastol)) | length(object@feastol) > 1){
        return("\nThe feasabile tolerance for convergence must be a real scalar.\n")
    }
    if(object@abstol < 0 & object@reltol < 0){
        return("\nAt least one of 'reltol' and 'abstol' must be positive.\n")
    }
    if(object@feastol <= 0){
        return("\nThe convergence criteria for feasability must be positive.\n")
    }
    if(object@stepadj > 1){
        return("\nThe step size adjustment must be in the interval (0, 1].\n")
    }
    if(object@stepadj <= 0){
        return("\nThe step size adjustment must be in the interval (0, 1].\n")
    }
    TRUE
}
setValidity("CTRL", validCTRL)
