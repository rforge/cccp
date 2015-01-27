##
## Function for creating 'CTRL' objects 
ctrl <- function(maxiters = 100L, abstol = 1e-7, reltol = 1e-6, feastol = 1e-7,
                 refine = FALSE, stepadj = 0.95, alpha = 0.01, beta = 0.5, maxreliter = 8L,
                 trace = TRUE, method = c("solve")){
    method <- match.arg(method)
    new("CTRL",
        maxiters = as.integer(maxiters),
        abstol = abstol,
        reltol = reltol,
        feastol = feastol,
        refine = refine,
        stepadj = stepadj,
        alpha = alpha,
        beta = beta,
        maxreliter = maxreliter,
        trace = trace,
        method = method)
}
##
## Function for creating 'NNOC' objects 
nnoc <- function(G, h){
    h <- new("NNOV", u = matrix(h), dims = nrow(G))
    new("NNOC", G = G, h = h, dims = nrow(G), vclass = "NNOV")    
}
##
## Function for creating 'SOCC' objects 
socc <- function(F, g, d, f){
    G <- matrix(0, nrow = nrow(F) + 1, ncol = ncol(F))
    G[1, ] <- -d
    G[-1, ] <- -F
    h <- new("SOCV", u = matrix(c(f, g), ncol = 1), dims = nrow(G))   
    new("SOCC", G = G, h = h, dims = nrow(G), vclass = "SOCV")    
}
##
## Function for creating 'PSDC' objects 
psdc <- function(Flist, F0){
    m <- nrow(Flist[[1]])
    G <- do.call("cbind", lapply(Flist, function(x) c(x)))
    h <- new("PSDV", u = matrix(as(F0, "vector"), ncol = 1), dims = m)
    new("PSDC", G = G, h = h, dims = m, vclass = "PSDV")        
}
##
## Function for creating 'NLFC' objects 
nlfc <- function(G, h){
    h <- new("NLFV", u = matrix(h), dims = nrow(G))
    new("NLFC", G = G, h = h, dims = nrow(G), vclass = "NLFV")    
}
