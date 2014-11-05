##
## All S4-class definitions 
##
## S4-class for optimization control parameters 
setClass("CTRL", representation(maxiters = "integer",
                                abstol = "numeric",
                                reltol = "numeric",
                                feastol = "numeric",
                                refine = "logical",
                                stepadj = "numeric",
                                alpha = "numeric",
                                beta = "numeric",
                                maxreliter = "integer",
                                trace = "logical",
                                method = "character"))
## S4-class for linear program definition
setClass("DEFLP", representation(q = "vector",
                                 A = "matrix",
                                 b = "matrix",
                                 cList = "list",
                                 n = "integer",
                                 k = "integer",
                                 ctrl = "CTRL",
                                 title = "character"),
         prototype = list(title = "Linear Program"))
## S4-class for quadratic program definition
setClass("DEFQP", representation(P = "matrix",
                                 q = "vector",
                                 A = "matrix",
                                 b = "matrix",
                                 cList = "list",
                                 n = "integer",
                                 k = "integer",
                                 ctrl = "CTRL",
                                 title = "character"),
         prototype = list(title = "Quadratic Program"))
## S4-class for linar program with non-linear constraints definition
setClass("DEFNL", representation(x0 = "vector",
                                 q = "vector",
                                 nlfList = "list",
                                 cList = "list",
                                 A = "matrix",
                                 b = "matrix",
                                 k = "integer",
                                 n = "integer",
                                 mnl = "integer",
                                 ctrl = "CTRL",
                                 H = "matrix",
                                 title = "character"),
         prototype = list(title = "Linear Program with non-linear constraints"))
## S4-class for convex program with non-linear constraints definition
setClass("DEFCP", representation(q = "vector",
                                 x0 = "vector",
                                 f0 = "function",
                                 nlfList = "list",
                                 cList = "list",
                                 A = "matrix",
                                 b = "matrix",
                                 k = "integer",
                                 n = "integer",
                                 mnl = "integer",
                                 ctrl = "CTRL",
                                 H = "matrix",
                                 title = "character"),
         prototype = list(title = "Convex Program with non-linear constraints"))
## S4-class union of convex programs
setClassUnion("CPD", members = c("DEFLP", "DEFQP", "DEFNL", "DEFCP"))
##
## S4-class for solving KKT-system
setClass("KKTSLV", representation(f = "function",
                                  items = "list"))
##
## S4-class for primal-dual variables
setClass("PDV", representation(x = "numeric",
                               y = "numeric",
                               s = "list",
                               z = "list",
                               kappa = "numeric",
                               tau = "numeric"),
         prototype = list(
             x = numeric(0),
             y = numeric(0),
             s = list(),
             z = list(),
             kappa = 1.0,
             tau = 1.0))
##
## S4-class for solution of convex problem
setClass("CPS", representation(pobj = "numeric",
                               dobj = "numeric",
                               dgap = "numeric",
                               rdgap = "numeric",
                               certp = "numeric",
                               certd = "numeric",
                               pslack = "numeric",
                               dslack = "numeric",
                               status = "character",
                               niter = "integer"
                               ),
         prototype = list(
             x = numeric(0),
             y = numeric(0),
             s = list(),
             z = list(),
             pobj = numeric(0),
             dobj = numeric(0),
             dgap = numeric(0),
             rdgap = numeric(0),             
             certp = numeric(0),
             certd = numeric(0),
             pslack = numeric(0),
             dslack = numeric(0),
             status = "unknown",
             niter = 0L), contains = "PDV")
##
## S4-classes with respect to cone constraints
##
## NNO-variable
setClass("NNOV", representation(
    u = "matrix",
    dims = "integer")
         )
## NNO-constraint
setClass("NNOC", representation(
    G = "matrix",
    h = "NNOV",
    dims = "integer",
    vclass = "character")
         )
## NNO-scaling
setClass("NNOS", representation(
    W = "list")
         )
## SOC-variable
setClass("SOCV", representation(
    u = "matrix",
    dims = "integer")
         )
## SOC-constraint
setClass("SOCC", representation(
    G = "matrix",
    h = "SOCV",
    dims = "integer",
    vclass = "character")
         )
## SOC-scaling
setClass("SOCS", representation(
    W = "list")
         )
## PSD-variable 
setClass("PSDV", representation(
    u = "matrix",
    dims = "integer")
         )
## SOC-constraint
setClass("PSDC", representation(
    G = "matrix",
    h = "PSDV",
    dims = "integer",
    vclass = "character")
         )
## SOC-scaling
setClass("PSDS", representation(
    W = "list")
         )
## NLF-variable
setClass("NLFV", representation(
    u = "matrix",
    dims = "integer")
         )
## NLF-constraint
setClass("NLFC", representation(
    G = "matrix",
    h = "NLFV",
    dims = "integer",
    vclass = "character")
         )
## NLF-scaling
setClass("NLFS", representation(
    W = "list")
         )
## S4-class union of variables pertinent to linear and nonlinear constraints
setClassUnion("LNLV", members = c("NNOV", "NLFV"))
