# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> ----
# >>>> This whole file is commented out <<<<<<<< ----
# >>>> because for compatibility with other projects <<<<<<<< ----
# >>>> we chose the dlill/conveniencefunctions <<<<<<<< ----
# >>>> implementation of norm_indiv <<<<<<<< ----
# >>>> However, to keep this code readily available  <<<<<<<< ----
# >>>> we still have it in the main R folder <<<<<<<< ----
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> ----

# 
# # normIndiv ----
# 
# # @export
# # @rdname normL2
# normIndiv <- function(data,
#                       prd0,
#                       errmodel = NULL,
#                       forcings = NULL,
#                       iiv = NULL,
#                       conditional = NULL,
#                       fixed.grid,
#                       nauxtimes = 500,
#                       cores = 1,
#                       deriv = TRUE,
#                       attr.name = "data"
# ) {
#   
#   # .. 1 Conditions ----- #
#   x.conditions <- names(fixed.grid)[-(1:2)]
#   e.conditions <- names(attr(errmodel, "mappings"))
#   prd0.conditions <- dMod::getConditions(prd0)
#   data.conditions <- names(data)
#   
#   
#   if (!all(data.conditions %in% x.conditions))
#     stop("The prediction function does not provide predictions for all conditions in the data.")
#   
#   cn_prd <- NULL
#   if (length(prd0.conditions) > 0) {
#     cn_prd <- prd0.conditions[1]
#     warning("prd0 contains named conditions. Only ", cn_prd, "will be used.\n",
#             "Preferably, use a general prediction function with no specific conditions.")
#   }
#   
#   # .. 2 Simulation times ----- #
#   timesD <- sort(unique(c(0, do.call(c, lapply(data, function(d) d$time)))))
#   timesD <- sort(union(timesD, seq(min(timesD), max(timesD), length.out = nauxtimes)))
#   
#   # .. 3 Optimization options ----- #
#   if (Sys.info()["sysname"] == "Windows" & cores > 1) {
#     warning("Parallelization of conditions on Windows not yet implemented. Unsing only 1 core.")
#     cores <- 1
#   }
#   # [] If more optizers come, one might want to define global variable with derivative-based/free optimizers
#   useDerivs <- deriv
#   
#   # .. 4 BLOQ options ----- #
#   opt.BLOQ = c("M3", "M4NM", "M4BEAL", "M1")
#   opt.hessian = c(
#     ALOQ_part1 = TRUE, ALOQ_part2 = TRUE, ALOQ_part3 = TRUE,
#     BLOQ_part1 = TRUE, BLOQ_part2 = TRUE, BLOQ_part3 = TRUE,
#     PD = TRUE  # enforce Hessian to be positive semidefinite, by setting nearest negative eigenvalues to zero
#   )
#   
#   # .. 5 Construct est.grid ----- #
#   est.grid <- build_est.grid(prd0 = prd0, fixed.grid = fixed.grid, conditional = conditional, condition.grid = attr(data, "condition.grid"), iiv = iiv)
#   parameters <- getParameters_est.grid(est.grid)
#   
#   # .. 6 Controls ----- #
#   controls <- list(timesD = timesD, attr.name = attr.name,
#                    conditions = x.conditions,
#                    opt.BLOQ = opt.BLOQ, opt.hessian = opt.hessian)
#   
#   # -------------------------------------------------------------------------#
#   # Objective Function ---- #
#   # -------------------------------------------------------------------------#
#   myfn <- function(...,
#                    fixed = NULL,
#                    deriv=useDerivs,
#                    conditions = controls$conditions,
#                    # env = NULL,
#                    opt.BLOQ    = controls$opt.BLOQ,
#                    opt.hessian = controls$opt.hessian,
#                    returnResiduals = FALSE) {
#     
#     # .. 1 BLOQ options ----- #
#     opt.hessian <- c(opt.hessian, controls$opt.hessian[setdiff(names(controls$opt.hessian), names(opt.hessian))])
#     opt.BLOQ    <- opt.BLOQ[1]
#     
#     # .. 2 Passed parameters ----- #
#     arglist    <- list(...)
#     arglist    <- arglist[dMod::match.fnargs(arglist, "pars")]
#     parsouter  <- arglist[[1]]
#     fixedouter <- fixed
#     
#     if (!is.null(attr(parsouter, "deriv"))){
#       parsouter <- unclass(parsouter)
#       attr(parsouter, "deriv") <- NULL
#       # TODO
#       # Implement chain rule for this function to allow function concatenation with to the right:
#       # obj_indiv*p(pars) does not work
#       # warning("The 'deriv' attribute of the pars-argument to obj_data has been set to NULL")
#     }
#     
#     # .. 3 Evaluate norm ----- #
#     outlist <- parallel::mclapply(conditions, function(cn) {
#       # outlist <- lapply(conditions, function(cn) {
#       
#       # .... 1 Parameters ------ #
#       # Fill general pars with individual values
#       dummy  <- make_pars(parsouter = parsouter, fixedouter = fixedouter, condition = cn, est.grid = est.grid,  fixed.grid = fixed.grid)
#       pars0  <- dummy$pars
#       fixed0 <- dummy$fixed
#       
#       
#       # .... 2 Prediction ------ #
#       prediction <- prd0(times = controls$timesD, pars = pars0, fixed = fixed0, deriv = deriv, conditions = cn_prd)
#       
#       if (nrow(prediction[[1]]) < length(timesD))
#         warning("Integrator has problems reaching tmax. Try increasing nauxtimes")
#       
#       
#       # .... 3 Calculate residuals ------ #
#       
#       # Skip condition if the data does not provide any observation for the prediction of the prd function.
#       if (length(intersect(data[[cn]][["name"]], colnames(prediction[[1]]))) == 0){
#         mywrss <- init_empty_objlist(pars0, deriv = deriv)
#         mywrss <- rename_objlist(mywrss, cn, est.grid)
#         return(mywrss)
#       }
#       err <- NULL
#       if (!is.null(errmodel) && (is.null(e.conditions) | (cn %in% e.conditions)))
#         err <- errmodel(out = prediction[[1]], pars = dMod::getParameters(prediction[[1]]), conditions = cn) 
#       nout <- res(data[[cn]], prediction[[1]], err[[cn]])
#       mywrss <- nll(nout = nout, pars = pars0, deriv = deriv)
#       chisquare <- attr(mywrss, "chisquare")
#       # .... 4 Rename general parnames into individual parnames ------ #
#       mywrss <- rename_objlist(mywrss, cn, est.grid)
#       if (returnResiduals) attr(mywrss, "residuals") <- cbind(condition = cn, nout)
#       return(mywrss)
#       # })
#     }, mc.cores = cores)
#     
#     # Remove list elements which could not be evaluated
#     failed <- vapply(outlist, is.character, FUN.VALUE = FALSE)
#     if (any(failed))
#       warning("Objective function could not be evaluated for conditions", paste0(conditions[failed], collapse = ", "))
#     outlist <- outlist[!failed]
#     if (length(outlist) == 0)
#       stop("Objective function could not be evaluated for any of the conditions.")
#     
#     # Extract values, gradients and hessians
#     values    <- lapply(outlist, function(x) x[["value"]]+1e6*sum(failed))
#     gradients <- lapply(outlist, function(x) x[["gradient"]])
#     hessians  <- lapply(outlist, function(x) x[["hessian"]])
#     residuals <- do.call(rbind, lapply(outlist, function(x) attr(x, "residuals")))
#     
#     # Generate output list with combined value, gradient and hessian
#     out <- init_empty_objlist(parsouter, deriv = deriv)
#     out[["value"]] <- Reduce("+", values)
#     if (deriv) {
#       for (grad in gradients) out$gradient[names(grad)]                 <- out$gradient[names(grad)] + grad
#       for (hes in hessians)   out$hessian[rownames(hes), colnames(hes)] <- out$hessian[rownames(hes), colnames(hes)] + hes
#     }
#     
#     
#     nearPD2 <- function (x, corr = FALSE){
#       # Ensure symmetry
#       X <- 0.5*(x+t(x))
#       # Eigenvalue decomposition
#       e <- eigen(X)
#       # Determine eigenvalues
#       # If all >0 then do nothing
#       # If minimum smaller than -1e-4 then do nothing
#       # If minimum larger than -1e-4 then do the algorithm
#       eV <- e$values
#       if (min(eV) > 0) return(X)
#       # Adjust negative eigenvalues to 0
#       eV[eV<0] <- 0
#       Xout <- e$vectors %*% diag(eV) %*% ginv(e$vectors)
#       if (corr) diag(Xout) <- 1
#       return(Xout)
#     }
#     
#     if (opt.hessian["PD"] & deriv) {
#       dn <- dimnames(out$hessian)
#       out$hessian <- nearPD2(out$hessian)
#       dimnames(out$hessian) <- dn
#     }
#     
#     # Combine contributions and attach attributes
#     attr(out, controls$attr.name) <- out$value
#     
#     # Add residuals if available
#     out[["residuals"]] <- residuals
#     
#     # Add chisquare to outputs
#     attr(out, "chisquare") <- Reduce("+", c(0, do.call(c,lapply(outlist, attr, which = "chisquare"))))
#     
#     return(out)
#   }
#   
#   
#   # -------------------------------------------------------------------------#
#   # Other Indiv functions ---- #
#   # -------------------------------------------------------------------------#
#   myparfn <- Reduce("+", lapply(controls$conditions, function(cn) {
#     parfn(
#       p2p = function(pars, fixed, deriv) {
#         dummy <- make_pars(parsouter = pars, fixedouter = fixed, condition = cn, est.grid = est.grid,  fixed.grid = fixed.grid)
#         pars0  <- dummy$pars
#         fixed0 <- dummy$fixed
#         as.parvec(c(pars0, fixed0))
#       },
#       parameters = parameters,
#       condition = cn
#     )
#   }))
#   
#   class(myfn) <- c("objfn", "fn")
#   attr(myfn, "conditions") <- data.conditions
#   attr(myfn, "parameters") <- parameters
#   attr(myfn, "modelname")  <- dMod::modelname(prd0, errmodel)
#   attr(myfn, "p_indiv")    <- myparfn
#   
#   return(myfn)
#   
#   
# }
# 
# 
# # Grid-related functions ----
# 
# # Build the est.grid from prd, fixed.grid, conditional and condition.grid
# #
# # Performs indiviudalization and localization of parameters
# #
# # @param prd0 prdfn for 1 condition. Extracts parameters via dMod::getParameters
# # @param fixed.grid As specified in template_project_sysfit.R: data.frame(parname, partask, ids...)
# # @param conditional As specified in template_project_sysfit.R: data.frame(parname, covname, covvalue)
# # @param condition.grid from datalist
# #
# # @return est.grid data.frame(parname, ids...)
# #
# # @examples
# # prd0 <- dMod::P(c(a = "exp(a + ETA_a)", "b" = "b", d = "d"))
# # fixed.grid <- stats::setNames(data.frame(parname = "d", partask = "Cond_specific", "1" = 1, "2" = NA, "3" = NA, stringsAsFactors = FALSE), c("parname", "partask", as.character(1:3)))
# # conditional <- data.frame(parname = c("b", "b", "d"), covname = "SEX", covvalue = c("male", "female", "male"), stringsAsFactors = FALSE)
# # condition.grid <- data.frame(ID = 1:3, SEX = c("female", "male", "male"), stringsAsFactors = FALSE)
# # iiv <- "ETA_a"
# # build_est.grid(prd0, fixed.grid, conditional, condition.grid, iiv)
# build_est.grid <- function(prd0, fixed.grid, conditional, condition.grid, iiv = NULL) {
#   # [] unit test?
#   check_cond(conditional, condition.grid)
#   
#   parameters     <- dMod::getParameters(prd0)
#   parameters_est <- setdiff(parameters, fixed.grid$parname)
#   parameters_est <- union(parameters_est, conditional$parname)
#   ids            <- condition.grid$condition
#   
#   est.grid <- matrix(parameters_est, nrow = length(parameters_est), ncol = nrow(condition.grid), dimnames = list(parameters_est, ids))
#   # 1 ETAs - iiv
#   if (length(iiv))
#     est.grid[iiv,] <- matrix(paste0(t(est.grid[iiv,,drop = FALSE]), "_", ids), ncol = length(ids), byrow = TRUE)
#   # 2 LOCmodel
#   if (length(conditional$parname)){
#     for (x in seq_along(conditional$parname)) {
#       # 1 Replace est pars by localized pars
#       cn      <- conditional[x,, drop = FALSE]
#       ids_est <- as.character(condition.grid$condition[as.character(condition.grid[[cn$covname]]) == as.character(cn$covvalue)])
#       est.grid[cn$parname,ids_est] <- paste0(est.grid[cn$parname,ids_est], "_", cn$covname, "_", cn$covvalue)
#       
#       if (length(intersect(cn$parname, fixed.grid$parname))){
#         # 2 Replace fixed localized pars by NA
#         ids_fix <- vapply(fixed.grid[fixed.grid$parname == cn$parname, -c(1:2)], function(y) !is.na(y), TRUE)
#         ids_fix <- names(ids_fix)[ids_fix]
#         est.grid[cn$parname,ids_fix] <- NA
#       }
#     }
#   }
#   # 3 Output data.frame
#   # [] Handle partask properly
#   est.grid <- cbind(data.frame(parname = parameters_est, partask = parameters_est, stringsAsFactors = FALSE),
#                     est.grid, stringsAsFactors = FALSE)
#   rownames(est.grid) <- NULL
#   
#   check_grids(fixed.grid, est.grid)
#   
#   est.grid
# }
# 
# 
# # Extract parameter names from est.grid
# #
# # @param est.grid data.frame(parname, partask, ids...)
# #
# # @return character of outer parameter names
# #
# # @examples
# # eg_good <- stats::setNames(data.frame(parname = "d", partask = "Cond_specific", "1" = NA, "2" = "dummy", "3" = "dummy", stringsAsFactors = FALSE), c("parname", "partask", as.character(1:3)))
# # getParameters_est.grid(eg_good)
# getParameters_est.grid <- function(est.grid) {
#   parameters <- unique(unlist(est.grid[-(1:2)], use.names = FALSE))
#   parameters <- parameters[!is.na(parameters)]
#   parameters <- sort(parameters)
#   parameters
# }
# 
# 
# 
# 
# # Extract individualized parameters for one ID 
# # 
# # From the outer individualized parsouter and fixedouter, get the generalized parsouter0, fixedouter0 to supply to prd0
# # 
# # Catch the following two cases:
# # * Entry in est.grid is NA                 => Parameter value present in fixed.grid
# # * Parameter name not present in parsouter => Parameter value present in fixed_outer
# # 
# # @param parsouter,fixedouter Arguments which are passed to the function returned by [normIndiv]
# # @param condition The specific condition for which the parameters should be extracted
# # @param est.grid,fixed.grid Lookup tables.
# #
# # @md
# # @return list(pars, fixed) to put into prd0
# #
# # @importFrom stats setNames
# make_pars <- function(parsouter, fixedouter, condition, est.grid,  fixed.grid) {
#   
#   # 1 Generate lookup for outer parameters
#   outer_lookup <- stats::setNames(est.grid[[condition]], est.grid[["parname"]])
#   outer_lookup <- outer_lookup[!is.na(outer_lookup)]
#   
#   # 2 Generate vector of all supplied "outer pars" to look things up
#   outer.vec <- c(parsouter, fixedouter)
#   
#   # 3.1 Look up est parameters
#   pars_lookup <- outer_lookup[outer_lookup %in% names(parsouter)]
#   pars0 <- NULL
#   if (length(pars_lookup))
#     pars0 <- stats::setNames(parsouter[pars_lookup], names(pars_lookup))
#   
#   # 3.2 Look up fixed parameters
#   fixed_lookup <- outer_lookup[outer_lookup %in% names(fixedouter)]
#   fixedouter <- NULL
#   if (length(fixed_lookup))
#     fixedouter <- stats::setNames(fixedouter[fixed_lookup], names(fixed_lookup))
#   
#   # 4 Get fixed parameters from fixed.grid and combine with fixedouter
#   fixed0 <- NULL
#   if (nrow(fixed.grid)){
#     fixed0 <- stats::setNames(fixed.grid[[condition]], fixed.grid[["parname"]])
#     fixed0 <- fixed0[!is.na(fixed0)]
#   }
#   fixed0 <- c(fixed0, fixedouter)
#   
#   # 5 Output
#   list(pars = pars0, fixed = fixed0)
# }
# 
# 
# # Do some consistency checks on conditional and condition.grid
# #
# # @param conditional see [normL2]
# # @param condition.grid condition.grid from data
# # @md
# check_cond <- function(conditional, condition.grid) {
#   if (length(setdiff(conditional$covname, names(condition.grid))))
#     stop("The following covariates to individualize parameters are not available: ", paste0(setdiff(conditional$covname, names(condition.grid)), collapse = ", "))
#   
#   cov_values <- unique(unlist(condition.grid[conditional$covname], use.names = FALSE))
#   if (length(setdiff(conditional$covvalue, cov_values)))
#     stop("The following covariate values to individualize parameters are not available: ", paste0(setdiff(conditional$covvalue, cov_values), collapse = ", "))
# }
# 
# 
# 
# # Some consistency checks for fixed.grid and est.grid
# #
# # * Checks for localized parameters appearing in both grids that exactly one NA is in either of the grids
# # * More to come ...
# #
# # @param fixed.grid,est.grid data.frame(parname, partask, ids...)
# #
# # @return TRUE: All tests passed, else an error is thrown
# #
# # @examples
# # fg <- stats::setNames(data.frame(parname = "d", partask = "Cond_specific", "1" = 1, "2" = NA, "3" = NA, stringsAsFactors = FALSE), c("parname", "partask", as.character(1:3)))
# # eg_good <- stats::setNames(data.frame(parname = "d", partask = "Cond_specific", "1" = NA, "2" = "dummy", "3" = "dummy", stringsAsFactors = FALSE), c("parname", "partask", as.character(1:3)))
# # eg_bad1 <- stats::setNames(data.frame(parname = "d", partask = "Cond_specific", "1" = NA, "2" = NA, "3" = "dummy", stringsAsFactors = FALSE), c("parname", "partask", as.character(1:3)))
# # eg_bad2 <- stats::setNames(data.frame(parname = "d", partask = "Cond_specific", "1" = "dummy", "2" = "dummy", "3" = "dummy", stringsAsFactors = FALSE), c("parname", "partask", as.character(1:3)))
# #
# # check_grids(fg, eg_good)
# # check_grids(fg, eg_bad1)
# # check_grids(fg, eg_bad2)
# check_grids <- function(fixed.grid, est.grid) {
#   par_overlap <- intersect(fixed.grid$parname, est.grid$parname)
#   if (length(par_overlap)){
#     # NA in both grids
#     both_na_check <- lapply(stats::setNames(nm = par_overlap), function(x) {
#       f <- fixed.grid[fixed.grid$parname == x, -(1:2)]
#       e <- est.grid[est.grid$parname == x, -(1:2)]
#       two_nas <- c(names(f)[is.na(f)], names(e)[is.na(e)])
#       two_nas[duplicated(two_nas)]
#     })
#     both_na_check <- do.call(c, both_na_check)
#     if (length(both_na_check)) stop("Parameter is NA in both fixed.grid and est.grid. The following parameters and conditions are affected", "\n", "\n", deparse(both_na_check))
#     
#     # Not NA in both grids
#     none_na_check <- lapply(stats::setNames(nm = par_overlap), function(x) {
#       f <- fixed.grid[fixed.grid$parname == x, -(1:2)]
#       e <- est.grid[est.grid$parname == x, -(1:2)]
#       no_nas <- c(names(f)[!is.na(f)], names(e)[!is.na(e)])
#       no_nas[duplicated(no_nas)]
#     })
#     none_na_check <- do.call(c, none_na_check)
#     if (length(none_na_check)) stop("Parameter is not NA in both fixed.grid and est.grid. The following parameters and conditions affected", "\n", "\n", deparse(none_na_check))
#   }
# }
# 
# 
# 
# 
# # Indiv-objlist related functions ----
# 
# 
# 
# # Rename the derivatives of an objective function according to a lookup table
# #
# # @param myobjlist dMod objective list
# # @param condition a specific condition (ID)
# # @param est.grid  estimation grid
# #
# # @return objective list
# #
# # @author Daniel Lill
# rename_objlist <- function(myobjlist, condition, est.grid) {
#   est_lookup <- stats::setNames(est.grid[[condition]], est.grid[["parname"]])
#   
#   if (!is.null(myobjlist$gradient)){
#     grad_names <- stats::setNames(names(myobjlist$gradient), names(myobjlist$gradient))
#     
#     est_lookup_used <- names(est_lookup)[names(est_lookup) %in% grad_names]
#     grad_names[est_lookup_used] <- est_lookup[est_lookup_used]
#     
#     names(myobjlist$gradient) <- grad_names
#   }
#   if (!is.null(myobjlist$hessian))
#     dimnames(myobjlist$hessian) <- list(grad_names, grad_names)
#   
#   return(myobjlist)
# }
# 
# 
# 
# # Create an objlist with zeros as entries
# # @param pars named vector. Only names and length are used
# # @param deriv TRUE or FALSE
# # @examples
# # init_empty_objlist(setNames(rnorm(5), letters[1:5]))
# init_empty_objlist <- function(pars, deriv = TRUE) {
#   
#   if (!deriv)
#     return(dMod::objlist(0,NULL,NULL))
#   
#   dMod::objlist(value = 0,
#                 gradient = setNames(rep(0, length(pars)), names(pars)),
#                 hessian = matrix(0, nrow = length(pars), ncol = length(pars),
#                                  dimnames = list(names(pars), names(pars))))
# }
# 
# 
# 
# # Generalized Inverse of a Matrix
# #
# # Calculates the Moore-Penrose generalized inverse of a matrix X.
# #
# # @param X Matrix for which the Moore-Penrose inverse is required.
# # @param tol A relative tolerance to detect zero singular values.
# # @return A MP generalized inverse matrix for X.
# # @family Auxiliary
# # @export
# ginv <- function(X, tol = sqrt(.Machine$double.eps)) {
#   if(length(dim(X)) > 2L || !(is.numeric(X) || is.complex(X)))
#     stop("'X' must be a numeric or complex matrix")
#   if(!is.matrix(X)) X <- as.matrix(X)
#   Xsvd <- svd(X)
#   if(is.complex(X)) Xsvd$u <- Conj(Xsvd$u)
#   Positive <- Xsvd$d > max(tol * Xsvd$d[1L], 0)
#   if (all(Positive)) Xsvd$v %*% (1/Xsvd$d * t(Xsvd$u))
#   else if(!any(Positive)) array(0, dim(X)[2L:1L])
#   else Xsvd$v[, Positive, drop=FALSE] %*% ((1/Xsvd$d[Positive]) * t(Xsvd$u[, Positive, drop=FALSE]))
# }
