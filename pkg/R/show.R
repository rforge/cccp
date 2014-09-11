##
## show-method for objects of S4-class 'QPD'
setMethod("show", signature = "CPD", function(object){
    title <- paste("*", object@title, "*")
    row <- paste(rep("*", nchar(title)), collapse = "")
    cat("\n")
    cat(row, "\n")
    cat(paste(title, "\n"))
    cat(row, "\n")
    cat("\n")
    cat(paste("Count of variables in objective:", object@n, "\n"))
    cat(paste("Count of equality constraints:", nrow(object@A), "\n"))
    if(class(object) == "DEFNL"){
        cat(paste("Count of nonlinear constraints:", object@mnl, "\n"))
    }
    cat(paste("Count of cone constraints:", object@k, "\n"))
    if(class(object) == "DEFNL"){
        cc <- unlist(lapply(object@cList, function(x) class(x)))
    } else {
        cc <- unlist(lapply(object@cList, function(x) class(x)))
    }
    cat("These consist of:\n")
    cat(paste("Constraints w.r.t. the nonnegative orthant:", max(0, sum(cc %in% "NNOC")), "\n"))
    cat(paste("Constraints w.r.t. the second-order cone:", max(0, sum(cc %in% "SOCC")), "\n"))
    cat(paste("Constraints w.r.t. the semidefinite cone:", max(0, sum(cc %in% "PSDC")), "\n"))
    cat("\n")
    cat("Use 'cps()' for finding a solution.\n")
})
##
## show-method for objects of S4-class 'CPS'
setMethod("show", signature = "CPS", function(object){
    title <- "* Solution of Convex Program *"
    row <- paste(rep("*", nchar(title)), collapse = "")
    cat("\n")
    cat(row, "\n")
    cat(title, "\n")
    cat(row, "\n")
    cat("\n")
    cat(paste("Value of primal objective:", signif(object@pobj), "\n"))
    if(length(object@dobj) > 0){
        cat(paste("Value of dual objective:", signif(object@dobj), "\n"))
    }
    if(length(object@dgap) > 0){
        cat(paste("Value of duality gap:", signif(object@dgap), "\n"))
    }
    if(length(object@rdgap) > 0){
        cat(paste("Value of relative duality gap:", signif(object@rdgap), "\n"))
    }
    if(length(object@certp) > 0){
        cat(paste("Certificate of primal infeasibility:", signif(object@certp), "\n"))
    }
    if(length(object@certd) > 0){
        cat(paste("Certificate of dual infeasibility:", signif(object@certd), "\n"))
    }
    if(length(object@pslack) > 0){
        cat(paste("Value of smallest primal slack:", signif(object@pslack), "\n"))
    }
    if(length(object@dslack) > 0){
        cat(paste("Value of smallest dual slack:", signif(object@pslack), "\n"))
    }
    cat(paste("Status of solution:", object@status, "\n"))
    cat(paste("Count of iterations:", object@niter, "\n\n"))
    cat("Solutions are contained in the slots 'x', 'y', 's' and 'z'.\n")
    cat("Use 'getx()', 'gety()', 'gets()' and 'getz()', respectively.\n")
})
