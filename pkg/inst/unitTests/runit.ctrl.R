##
## Unit testing of CTRL-objects
test.CTRL <- function(){
    ctrlObj <- ctrl()
    checkTrue(validObject(ctrlObj))
    checkTrue(class(ctrlObj) == "CTRL")
    checkException(ctrl(maxiter = -5.7))
    checkException(ctrl(maxiter = -5L))
    checkException(ctrl(maxiter = 0L))
    checkException(ctrl(abstol = c(1e-7, 1e-8)))
    checkException(ctrl(reltol = c(1e-7, 1e-8)))
    checkException(ctrl(feastol = c(1e-7, 1e-8)))
    checkException(ctrl(abstol = -1e-7, reltol = -1e-8))
    checkException(ctrl(feastol = -1e-7))
    checkException(ctrl(method = "lbm"))
    checkException(ctrl(stepadj = 2))    
    checkException(ctrl(stepadj = 0))        
    checkException(ctrl(trace = 2))    
    return()
}
