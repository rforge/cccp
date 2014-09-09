##
## Main function for defining and solving linear and quadratic programs with cone constraints
cccp <- function(P = NULL, q = NULL, A = NULL, b = NULL, conecon = list(),
                 optctrl = ctrl()){

    if(is.null(P) & is.null(q)){
        stop("At least, 'P' or 'q' must be provided for quadratic or linear objective.\n")
    }
    if(length(conecon) > 0){
        coneclasses <- unlist(lapply(conecon, class))
        if(!all(coneclasses %in% c("NNOC", "SOCC", "PSDC"))){
            stop("List elements of cone constraints must be of either class 'NNOC', or 'SOCC', or 'PSDC'.\n")
        }
    }
    cpdef <- cpd(P, q, A, b, conecon, optctrl)
    cpsol <- cps(cpdef)
    return(cpsol)
}
