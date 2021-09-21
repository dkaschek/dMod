# -------------------------------------------------------------------------#
# normIndiv - internals ----
# -------------------------------------------------------------------------#

#' @export
unclass_parvec <- function(x) {setNames(unclass(x)[1:length(x)], names(x))}

#' Extract pars, fixed and parnames from grids for a given condition
#'
#' @param est.vec 
#' @param est.grid needs condition and ID
#' @param fix.grid 
#' @param ID 
#'
#' @return
#' @export
#'
#' @examples
#' pars <- c("k1_A" = 1, "k1_B" = 1, "k2_A" = 2, "k3" = 3, "k4" = 4)
#' fixed <- c("k2_B" = 2.5)
#' 
#' est.grid <- data.frame(ID = 1:2,
#'                        condition = c("A", "B"),
#'                        k1 = c("k1_A", "k1_B"),
#'                        k2 = c("k2_A", "k2_B"),
#'                        k3 = c("k3", NA),
#'                        k4 = c("k4", "k4"),
#'                        stringsAsFactors = FALSE)
#' fix.grid <- data.frame(ID = 1:2,
#'                        condition = c("A", "B"),
#'                        k3 = c(NA, 3.5),
#'                        k5 = c(5.1,5.2),
#'                        k6 = c(6,6),
#'                        stringsAsFactors = FALSE)
#' make_pars(pars, fixed, est.grid, fix.grid, 1)
#' make_pars(pars, fixed, est.grid, fix.grid, 2)
make_pars <- function(pars, fixed = NULL, est.grid, fix.grid, ID){
  
  i <- ID
  est.grid <- as.data.table(est.grid)
  fix.grid <- as.data.table(fix.grid)
  
  pars        <- unclass_parvec(pars)
  fixed       <- unclass_parvec(fixed)
  
  pars_outer  <- pars
  fixed_outer <- fixed
  
  # Lookup table for supplied grid.outer (entries in grid) to grid.inner (names of grid)
  parnames  <- unlist(est.grid[ID == i, !c("ID", "condition")])
  parnames <- parnames[!is.na(parnames)]
  
  # Get Parameters from grids
  # Look up names of all supplied
  supplied <- c(pars, fixed)
  supplied <- setNames(supplied[parnames], names(parnames))
  # Get fixed
  fixed <- unlist(fix.grid[ID == i, !c("ID", "condition")])
  fixed <- fixed[!is.na(fixed)]
  
  # Sort supplied "fixed" parameters back to fixed
  fixed <- c(fixed, supplied[parnames %in% names(fixed_outer)])
  pars <- supplied[!parnames %in% names(fixed_outer)]
  parnames_full <- parnames
  parnames <- parnames[!parnames %in% names(fixed_outer)]
  
  return(list(pars = unlist(pars), fixed = unlist(fixed), parnames = parnames, parnames_full = parnames_full))
}


#' Update attr(prediction, "deriv") to correct par.grid.outer names
#' 
#'
#' @param pred0 prd0(times,pars)[[1]]
#' @param pars pars
#' @param est.grid est.grid
#' @param cn name of condition
#'
#' @return pred0 with updated deriv argument
#' 
#' @export
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
renameDerivPars <- function(pred0, pars, est.grid, cn) {
  der <- attr(pred0, "deriv")
  dernm <- setdiff(colnames(der), "time")
  derpars <- do.call(rbind,strsplit(dernm, ".", fixed = TRUE))
  derpars <- data.table(derpars)
  setnames(derpars, c("y", "parinner"))
  
  eg <- data.table(est.grid)[condition == cn, !c("condition", "ID")]
  eg <- data.table(t(eg), keep.rownames = "parinner")
  setnames(eg, "V1", "parouter")
  
  derpars <- eg[derpars, on ="parinner"]
  derpars[,`:=`(dernmnew = paste0(y, ".", parouter))]
  
  colnames(der) <- c("time", derpars$dernmnew)
  der <- sumDuplicatedParsInDeriv(der)
  
  attr(pred0, "deriv") <- der
  pred0
}

#' Sum redundant columns
#'
#' Remove redundancies (happens when a parameter is duplicated and mapped to the same outer parameter). 
#' For example in est.grid = data.table(ID = 1, init_A = "Atot", Atot= "Atot", ...)
#'
#' @param der deriv of a single condition, as used in [renameDerivPars]
#'
#' @return der with redundant columns summed and duplicates removed
sumDuplicatedParsInDeriv <- function(der) {
  nm <- colnames(der)
  isDupe <- duplicated(nm)
  if (!any(isDupe)) {
    return(der)
  }
  nmDupe <- nm[isDupe]
  for (nmx in nmDupe) {
    for (i in 1:nrow(der))
      der[i,nmx] <- sum(der[i,nm == nmx])
  }
  der[, !isDupe, drop = FALSE]
}


#' Rename the gradient and hessian of an objlist
#'
#' @param objlist an objlist with grid.inner parnames 
#' @param parnames vector(grid.inner = grid.outer), e.g. as the `parnames` element from the result of [make_pars()]
#'
#' @return objlist with new names
#' @export
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#'
#' @examples
#' ol <- dMod:::init_empty_objlist(c(a = 1, b = 2))
#' parnames <- c(a = "c", b = "d")
#' renameDerivParsInObjlist(ol, parnames)
renameDerivParsInObjlist <- function(objlist, parnames) {
  objlist$gradient <- objlist$gradient[names(parnames)]
  names(objlist$gradient) <- unname(parnames)
  
  objlist$hessian <- objlist$hessian[names(parnames),names(parnames), drop = FALSE]
  dimnames(objlist$hessian) <- list(unname(parnames), unname(parnames))
  
  objlist <- sumDuplicatedParsInObjlist(objlist)
  
  objlist
}

#' Remove redundant outer names
#'
#' Remove redundancies (happens when a parameter is duplicated and mapped to the same outer parameter). 
#' For example in est.grid = data.table(ID = 1, init_A = "Atot", Atot= "Atot", ...)
#' 
#' @param objlist objlist with potentially duplicated names
#'
#' @return objlist with duplicated gradient and hessian elements summed and redundancies removed
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#'
#' @examples
#' ol <- dMod:::init_empty_objlist(c("S2" = 2, "S3" = 3, S2 = 3))
#' ol$gradient <- ol$gradient + 1:3
#' ol$hessian <- ol$hessian + 1:9
#' sumDuplicatedParsInObjlist(ol)
sumDuplicatedParsInObjlist <- function(ol) {
  attrs0 <- attributes(ol)
  attrs0 <- attrs0[setdiff(names(attrs0), c("names", "class"))]
  
  nm <- names(ol$gradient)
  isDupe <- duplicated(nm)
  if (!any(isDupe)) {
    return(ol)
  }
  nmDupe <- nm[isDupe]
  for (nmx in nmDupe) {
    ol$gradient[nmx] <- sum(ol$gradient[nm == nmx])
    for (i in 1:ncol(ol$hessian))
      ol$hessian[nmx,i] <- sum(ol$hessian[nm == nmx,i])
    for (i in 1:nrow(ol$hessian))
      ol$hessian[i,nmx] <- sum(ol$hessian[i,nm == nmx])
  }
  ol <- objlist(ol$value, ol$gradient[!isDupe], ol$hessian[!isDupe, !isDupe, drop = FALSE])
  
  attributes(ol) <- c(attributes(ol), attrs0)
  ol
}

#' Add single-valued parameters to the pargrids
#'
#' @param pars character or numeric
#' @param gridlist list(fix.grid, est.grid)
#' @param FLAGoverwrite Overwrite existing parameters in fix.grid or est.grid?
#' 
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#'
#' @return gridlist
#' @export 
#' 
#' @examples 
#' pars <- c(NewParSymbolic = "NewParSymbolic", NewParFixed = 1)
#' est.grid <- data.table(ID = 1:2,
#'                        condition = c("A", "B"),
#'                        k1 = c("k1_A", "k1_B"),
#'                        k2 = c("k2_A", "k2_B"),
#'                        k3 = c("k3", NA),
#'                        k4 = c("k4", "k4"),
#'                        stringsAsFactors = FALSE)
#' fix.grid <- data.table(ID = 1:2,
#'                        condition = c("A", "B"),
#'                        k3 = c(NA, 3.5),
#'                        k5 = c(5.1,5.2),
#'                        k6 = c(6,6),
#'                        stringsAsFactors = FALSE)
#' indiv_addGlobalParsToGridlist(c(NewParSymbolic = "NewParSymbolic", NewParFixed = 1), gridlist(est.grid = est.grid, fix.grid = fix.grid))
#' indiv_addGlobalParsToGridlist(c(k1 = 1), gridlist(est.grid = est.grid, fix.grid = fix.grid), FLAGoverwrite = FALSE) # nothing happens
#' indiv_addGlobalParsToGridlist(c(k1 = 1), gridlist(est.grid = est.grid, fix.grid = fix.grid), FLAGoverwrite = TRUE) # k1 is replaced and moved to fix.grid
indiv_addGlobalParsToGridlist <- function(pars, gridlist, FLAGoverwrite = FALSE) {
  # 1 Get grids
  est.grid <- gridlist$est.grid
  fix.grid <- gridlist$fix.grid
  # 2 Determine overwriting action: Cut down grids or pars
  if (FLAGoverwrite) {
    est.grid <- est.grid[,.SD,.SDcols = setdiff(names(est.grid), names(pars))]
    fix.grid <- fix.grid[,.SD,.SDcols = setdiff(names(fix.grid), names(pars))]
  } else {
    pars <- pars[setdiff(names(pars), names(est.grid))]
    pars <- pars[setdiff(names(pars), names(fix.grid))]
  }
  # 3 determine where to add the parameters to
  is_num  <- !is.na(as.numeric(pars))
  # 4 Add parameters to grid
  if (sum(!is_num)) est.grid <- data.table(est.grid, as.data.table(t(pars[!is_num])))
  if (sum( is_num)) fix.grid <- data.table(fix.grid, as.data.table(t(setNames(as.numeric(pars[ is_num]), names(pars[ is_num])))))
  
  list(est.grid = est.grid, fix.grid = fix.grid)
}


#' Add individualized parameters to grids
#' 
#' Can only add parameters, cannot update existing parameters
#'
#' @param pars data.table(ID, condition, pars...) Only one of ID, condition must be present
#' @param gridlist [dMod::gridlist()]
#'
#' @return modified gridlist
#' 
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#' @export
#'
#' @examples
#' pars <- data.table(ID = 1:2, newFix = c(1,2), newEst = c("a", "b"), newMix = c("a", 1))
#' est.grid <- data.table(ID = 1:2,
#'                        condition = c("A", "B"),
#'                        k1 = c("k1_A", "k1_B"),
#'                        k2 = c("k2_A", "k2_B"),
#'                        k3 = c("k3", NA),
#'                        k4 = c("k4", "k4"),
#'                        stringsAsFactors = FALSE)
#' fix.grid <- data.table(ID = 1:2,
#'                        condition = c("A", "B"),
#'                        k3 = c(NA, 3.5),
#'                        k5 = c(5.1,5.2),
#'                        k6 = c(6,6),
#'                        stringsAsFactors = FALSE)
#' gl <- gridlist(est.grid = est.grid, fix.grid = fix.grid)
#' indiv_addLocalParsToGridList(pars, gl)
indiv_addLocalParsToGridList <- function(pars, gridlist, FLAGoverwrite = FALSE) {
  est.grid <- gridlist$est.grid
  fix.grid <- gridlist$fix.grid
  
  # Determine which cols to join on. Assuming ID and condition are present in fix.grid and est.grid
  joincols <- intersect(c("ID", "condition"), names(pars))
  parscols <- setdiff(names(pars), joincols)
  setcolorder(pars, joincols)
  
  # 2 Determine overwriting action: Cut down grids or pars
  if (FLAGoverwrite) {
    est.grid <- est.grid[,.SD,.SDcols = c(joincols, setdiff(names(est.grid), names(pars)))]
    fix.grid <- fix.grid[,.SD,.SDcols = c(joincols, setdiff(names(fix.grid), names(pars)))]
  } else {
    pars <- pars[,.SD,.SDcols = c(joincols, setdiff(names(pars), names(est.grid)))]
    pars <- pars[,.SD,.SDcols = c(joincols, setdiff(names(pars), names(fix.grid)))]
  }
  
  # Split pars into fix and est, introduce NAs for mixed parameters
  pars_fix <- copy(pars)
  # power move: delete all symbolic columns, replace symbols with NA in remaining cols
  pars_fix[,(parscols) := lapply(.SD, function(x) {
    x <- as.numeric(x)
    if (all(is.na(x))) {
      return(NULL)
    } else x}), .SDcols = parscols]
  fix.grid <- pars_fix[fix.grid, on = joincols]
  
  pars_est <- copy(pars)
  # power move: delete all numeric columns, replace numbers with NA in remaining cols
  pars_est[,(parscols) := lapply(.SD, function(x) {
    numidx <- !is.na(as.numeric(x)); 
    if (all(numidx)) {
      return(NULL)
    } else replace(x, numidx, NA_character_)}), .SDcols = parscols]
  est.grid <- pars_est[est.grid, on = joincols]
  
  gridlist(est.grid = est.grid, fix.grid = fix.grid)
}


#' Create an objlist with zeros as entries
#' @param pars named vector. Only names and length are used
#' @param deriv TRUE or FALSE
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#' @examples
#' init_empty_objlist(setNames(rnorm(5), letters[1:5]))
#' init_empty_objlist(setNames(rnorm(5), letters[1:5]), FLAGchisquare = TRUE)
init_empty_objlist <- function(pars, deriv = TRUE, FLAGchisquare = FALSE) {
  
  gr <- if (deriv) setNames(rep(0, length(pars)), names(pars)) else NULL
  he <- if (deriv) matrix(0, nrow = length(pars), ncol = length(pars), dimnames = list(names(pars), names(pars))) else NULL
  
  out <- dMod::objlist(value = 0, gradient = gr, hessian = he)
  if (FLAGchisquare) attr(out, "chisquare") <- 0
  out
}



#' Collect grids in list
#' 
#' Ensure all tables are data.tables
#' 
#' @param est.grid,fix.grid data.tables, will be coerced to one
#'
#' @return list of the grid
#' @export
gridlist <- function(est.grid, fix.grid) {
  est.grid <- as.data.table(est.grid)
  fix.grid <- as.data.table(fix.grid)
  list(est.grid = est.grid, fix.grid = fix.grid)
}



# -------------------------------------------------------------------------#
# normIndiv - externals ----
# -------------------------------------------------------------------------#

#' Title
#'
#' @param prd0 
#' @param est.grid 
#' @param fix.grid 
#'
#' @return
#' @export
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
PRD_indiv <- function(prd0, est.grid, fix.grid) {
  
  if (!is.data.table(est.grid)) warning("est.grid was coerced to data.table (was", class(est.grid), ")")
  if (!is.data.table(fix.grid)) warning("fix.grid was coerced to data.table (was", class(fix.grid), ")")
  
  est.grid <- data.table(est.grid)
  fix.grid <- data.table(fix.grid)
  
  setkeyv(est.grid, c("ID", "condition"))
  setkeyv(fix.grid, c("ID", "condition"))
  
  # Title
  #
  # @param times 
  # @param pars 
  # @param fixed 
  # @param deriv 
  # @param conditions 
  # @param FLAGbrowserN 0: Don't debug, 1: Debug when there is an error, 2: always debug
  # @param FLAGverbose 
  # @param FLAGrenameDerivPars Needed for datapointL2_indiv, where I need derivs wrt the outer parameters. Don't remember what this FLAG could ever be used for except for being TRUE
  #
  # @return
  # @export
  #
  # @examples
  prd <- function(times, pars, fixed = NULL, deriv = FALSE, conditions = est.grid$condition, 
                  FLAGbrowserN = 0, 
                  FLAGverbose = FALSE,
                  FLAGrenameDerivPars = TRUE
  ) {
    out <- lapply(setNames(nm = conditions), function(cn) {
      if (FLAGbrowserN == 2) browser()
      ID <- est.grid[condition == cn, ID]
      if (FLAGverbose) cat(ID, cn, "\n", sep = " ---- ")
      
      dummy <- make_pars(pars, fixed, est.grid, fix.grid, ID)
      pars_ <- dummy$pars
      fixed_ <- dummy$fixed
      
      if (length(missingPars <- setdiff(getParameters(prd0), c("time", names(c(pars_, fixed_))))))
        stop("The following parameters are missing: ", paste0(missingPars, collapse = ", "))
      
      pred0 <-try(prd0(times, pars_, fixed = fixed_, deriv = deriv, conditions = NULL)[[1]])
      if (inherits(pred0, "try-error") && FLAGbrowserN == 1) {
        browser()
        # Try this code to debug your model
        # 1 Parameters
        pinner <- p(pars_, fixed = fixed_)
        compare(names(pinner[[1]]), getParameters(x)) #setdiff(y,x) should be empty!
        # 2 ode-model
        pinner_test <- setNames(runif(length(getParameters(x))),getParameters(x))
        x(times, pinner_test, deriv = FALSE)
      }
      if (deriv && FLAGrenameDerivPars) pred0 <- renameDerivPars(pred0, pars, est.grid, cn)
      pred0
    })
    as.prdlist(out)
  }
  class(prd) <- c("prdfn", "fn")
  prd
}

#' Title
#'
#' @param p0 
#' @param est.grid 
#' @param fix.grid 
#'
#' @return
#' @export
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
P_indiv <- function(p0, est.grid, fix.grid) {
  
  if (!is.data.table(est.grid)) warning("est.grid was coerced to data.table (was", class(est.grid), ")")
  if (!is.data.table(fix.grid)) warning("fix.grid was coerced to data.table (was", class(fix.grid), ")")
  
  est.grid <- data.table(est.grid)
  fix.grid <- data.table(fix.grid)
  setkeyv(est.grid, c("ID", "condition"))
  setkeyv(fix.grid, c("ID", "condition"))
  
  
  # @param FLAGbrowser 0: Don't debug, >= 1: debug
  prd <- function(pars, fixed = NULL, deriv = FALSE, conditions = est.grid$condition, 
                  FLAGbrowser = FALSE, 
                  FLAGverbose = FALSE) {
    out <- lapply(setNames(nm = conditions), function(cn) {
      if (FLAGbrowser) browser()
      ID <- est.grid[condition == cn, ID]
      if (FLAGverbose) cat(ID, cn, "\n", sep = " ---- ")
      dummy <- make_pars(pars, fixed, est.grid, fix.grid, ID)
      pars_ <- dummy$pars
      fixed_ <- dummy$fixed
      
      missingPars <- setdiff(getParameters(p0), names(c(pars_, fixed_)))
      if (length(missingPars))
        stop("The following parameters are missing: ", paste0(missingPars, collapse = ", "))
      p0(pars_, fixed = fixed_, deriv = deriv)[[1]]
    })
    out
  }
  class(prd) <- c("parfn", "fn")
  prd
}



#' Fast normL2 
#'
#' @param data 
#' @param prd0 
#' @param errmodel 
#' @param est.grid 
#' @param fix.grid 
#' @param times 
#' @param attr.name 
#'
#' @return objective function
#' @export
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#'
#' @importFrom parallel mclapply
normL2_indiv <- function (data, prd0, errmodel = NULL, est.grid, fix.grid, times = NULL, attr.name = "data") {
  
  if (!is.data.table(est.grid)) warning("est.grid was coerced to data.table (was", class(est.grid), ")")
  if (!is.data.table(fix.grid)) warning("fix.grid was coerced to data.table (was", class(fix.grid), ")")
  
  est.grid <- data.table(est.grid)
  fix.grid <- data.table(fix.grid)
  setkeyv(est.grid, c("ID", "condition"))
  setkeyv(fix.grid, c("ID", "condition"))
  
  
  timesD <- sort(unique(c(0, do.call(c, lapply(data, function(d) d$time)))))
  if (!is.null(times)) 
    timesD <- sort(union(times, timesD))
  x.conditions <- est.grid$condition
  data.conditions <- names(data)
  e.conditions <- names(attr(errmodel, "mappings"))
  controls <- list(times = timesD, attr.name = attr.name, conditions = intersect(x.conditions, 
                                                                                 data.conditions))
  force(errmodel)
  force(fix.grid)
  force(est.grid)
  
  myfn <- function(..., fixed = NULL, deriv = TRUE, conditions = controls$conditions, simcores = 1, 
                   FLAGbrowser = FALSE, 
                   FLAGbrowser2 = FALSE, 
                   FLAGverbose = FALSE,
                   FLAGNaNInfwarnings = FALSE,
                   FixedConditions = fixed.conditions) {
    arglist <- list(...)
    arglist <- arglist[match.fnargs(arglist, "pars")]
    pars <- arglist[[1]]
    
    objlists <- lapply(setNames(nm = conditions), function(cn) {
      if (FLAGbrowser) browser()
      
      ID <- est.grid[condition == cn, ID]
      if (FLAGverbose) cat(ID, cn, "\n", sep = " ---- ")
      dummy <- make_pars(pars, fixed, est.grid, fix.grid, ID)
      pars_ <- dummy$pars
      fixed_ <- dummy$fixed
      
      if (!length(pars_)) return(init_empty_objlist(pars, deriv = deriv, FLAGchisquare = TRUE)) # No pars_ can happen if one fits only condition specific parameters and in this condition there are none
      
      prediction <- try(prd0(times = controls$times, pars = pars_, fixed = fixed_, deriv = deriv))
      
      if (inherits(prediction, "try-error")) 
        stop("Prediction failed in \n>>>condition = ", cn, "\n>>>ID = ", ID, "\n\nTry iterating p(pars), (x*p)(pars), ... to find the problem.")
      
      prediction <- prediction[[1]]
      prediction <- check_and_sanitize_prediction(prediction, data, cn, FLAGNaNInfwarnings)
      
      err <- NULL
      if (any(is.na(data[[cn]]$sigma))) {
        err <- errmodel(out = prediction, pars = getParameters(prediction), conditions = cn, deriv=deriv)
        mywrss <- nll(res(data[[cn]], prediction, err[[1]]), deriv = deriv, pars = pars)
      } else {
        mywrss <- nll(res(data[[cn]], prediction), deriv = deriv, pars = pars)
      }
      
      if (deriv) mywrss <- renameDerivParsInObjlist(mywrss, dummy$parnames) 
      
      mywrss
    })
    
    # Sum all objlists
    out <- Reduce("+", c(list(init_empty_objlist(c(pars, fixed), deriv = deriv)), objlists))
    
    # Consider fixed: return only derivs wrt pouter
    out$gradient <- out$gradient[names(pars)]
    out$hessian <- out$hessian[names(pars), names(pars)]
    
    # Populate attributes
    attr(out, controls$attr.name) <- out$value
    ll_conditions <- data.frame(
      logl = vapply(setNames(objlists, conditions), function(.x) .x$value, 1),
      chi2 = vapply(setNames(objlists, conditions), function(.x) attr(.x, "chisquare"), 1))
    ll_sum <- data.frame(logl = sum(ll_conditions$logl),
                         chi2 = sum(ll_conditions$chi2))
    attributes(out) <- c(attributes(out), list(ll_cond_df = ll_conditions))
    attributes(out) <- c(attributes(out), list(ll_sum_df = ll_sum))
    # attr(out, "AIC") <- out$value + length(pars) * 2
    # attr(out, "BIC") <- out$value + length(pars) * log(nrow(as.data.frame(data)))
    return(out)
  }
  
  class(myfn) <- c("objfn", "fn")
  attr(myfn, "conditions") <- data.conditions
  attr(myfn, "parameters") <- attr(prd0, "parameters")
  attr(myfn, "modelname") <- modelname(prd0, errmodel)
  return(myfn)
}

#' DatapointL2 without access to stored predictions
#'
#' @param name 
#' @param time 
#' @param value 
#' @param sigma 
#' @param attr.name 
#' @param condition 
#' @param prd_indiv 
#'
#' @return
#' @export
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
datapointL2_indiv <- function (name, time, value, sigma = 1, attr.name = "validation", 
                               condition, prd_indiv) {
  controls <- list(mu = structure(name, names = value)[1], 
                   time = time[1], sigma = sigma[1], attr.name = attr.name)
  
  myfn <- function(..., fixed = NULL, deriv = TRUE, conditions = NULL, FLAGbrowser = FALSE, SIMOPT.times = seq(0,time, length.out = 100)) {
    
    if (FLAGbrowser)
      browser()
    
    mu <- controls$mu
    time <- controls$time
    sigma <- controls$sigma
    attr.name <- controls$attr.name
    arglist <- list(...)
    arglist <- arglist[match.fnargs(arglist, c("pars"))]
    
    times <- sort(c(unique(SIMOPT.times, time)))
    pouter <- arglist[[1]]
    prediction <- prd_indiv(times, pouter, fixed = fixed, condition = condition, deriv = deriv, FLAGrenameDeriv = TRUE)
    if (!is.null(conditions) && !condition %in% conditions) 
      return()
    
    if (is.null(conditions) && !condition %in% names(prediction)) 
      stop("datapointL2 requests unavailable condition. Call the objective function explicitly stating the conditions argument.")
    
    datapar <- setdiff(names(mu), names(fixed))
    parapar <- setdiff(names(pouter), c(datapar, names(fixed)))
    time.index <- which(prediction[[condition]][, "time"] == time)
    if (length(time.index) == 0) 
      stop("datapointL2() requests time point for which no prediction is available. Please add missing time point by the times argument in normL2()")
    withDeriv <- !is.null(attr(prediction[[condition]], 
                               "deriv"))
    pred <- prediction[[condition]][time.index, ]
    deriv <- NULL
    if (withDeriv) 
      deriv <- attr(prediction[[condition]], "deriv")[time.index, ]
    pred <- pred[mu]
    
    # prd_indiv doesn't return the same derivpars in all conditions 
    derivnm_split <- strsplit(names(deriv), "\\.")
    derivnm_split <- lapply(derivnm_split, function(x) x[2])
    derivnm_split <- do.call(c, derivnm_split)
    parapar <- parapar[parapar %in% derivnm_split]
    
    if (withDeriv) {
      mu.para <- intersect(paste(mu, parapar, sep = "."), 
                           names(deriv))
      deriv <- deriv[mu.para]
    }
    res <- as.numeric(pred - c(fixed, pouter)[names(mu)])
    val <- as.numeric((res/sigma)^2)
    gr <- NULL
    hs <- NULL
    

    if (withDeriv) {
      dres.dp <- structure(rep(0, length(pouter)), names = names(pouter))
      if (length(parapar) > 0) 
        dres.dp[parapar] <- as.numeric(deriv)
      if (length(datapar) > 0) 
        dres.dp[datapar] <- -1
      gr <- 2 * res * dres.dp/sigma^2
      hs <- 2 * outer(dres.dp, dres.dp, "*")/sigma^2
      colnames(hs) <- rownames(hs) <- names(pouter)
    }
    out <- objlist(value = val, gradient = gr, hessian = hs)
    attr(out, controls$attr.name) <- out$value
    # attr(out, "prediction") <- pred
    return(out)
  }
  class(myfn) <- c("objfn", "fn")
  attr(myfn, "conditions") <- condition
  attr(myfn, "parameters") <- value[1]
  return(myfn)
}


#' DatapointL2 without env access
#'
#' @param name 
#' @param time 
#' @param value 
#' @param sigma 
#' @param attr.name 
#' @param condition 
#' @param prd_indiv 
#'
#' @return
#' @export
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
timepointL2_indiv <- function(name, time, value, sigma = 1, attr.name = "timepointL2", 
                              condition, prd_indiv) {
  
  # [] mu needs to be numeric and time needs tocharacter
  controls <- list(mu = structure(name, names = value)[1], 
                   time = time[1], sigma = sigma[1], attr.name = attr.name)
  
  myfn <- function(..., fixed = NULL, deriv = TRUE, conditions = NULL, env = NULL) {
    
    mu        <- controls$mu
    time      <- controls$time
    timepar <- 
      sigma     <- controls$sigma
    attr.name <- controls$attr.name
    
    arglist <- list(...)
    arglist <- arglist[match.fnargs(arglist, c("times", "pars"))]
    # ensure time point has prediction
    times      <- arglist[[1]]
    times      <- sort(c(unique(times, time)))
    pouter     <- arglist[[2]]
    prediction <- prd_indiv(times, pouter, condition = condition, deriv = deriv)
    if (!is.null(conditions) && !condition %in% conditions) 
      return()
    
    if (is.null(conditions) && !condition %in% names(prediction)) 
      stop("datapointL2 requests unavailable condition. Call the objective function explicitly stating the conditions argument.")
    
    datapar    <- setdiff(names(mu), names(fixed))
    parapar    <- setdiff(names(pouter), c(datapar, names(fixed)))
    time.index <- which(prediction[[condition]][, "time"] == time)
    
    withDeriv <- !is.null(attr(prediction[[condition]], "deriv"))
    pred      <- prediction[[condition]][time.index, ]
    deriv     <- NULL
    if (withDeriv) 
      deriv <- attr(prediction[[condition]], "deriv")[time.index,]
    
    pred <- pred[mu]
    if (withDeriv) {
      mu.para <- intersect(paste(mu, parapar, sep = "."), names(deriv))
      deriv <- deriv[mu.para]
    }
    
    res <- as.numeric(pred - c(fixed, pouter)[names(mu)])
    val <- as.numeric((res/sigma)^2)
    gr <- NULL
    hs <- NULL
    if (withDeriv) {
      dres.dp <- structure(rep(0, length(pouter)), names = names(pouter))
      if (length(parapar) > 0) 
        dres.dp[parapar] <- as.numeric(deriv)
      if (length(datapar) > 0) 
        dres.dp[datapar] <- -1
      gr <- 2 * res * dres.dp/sigma^2
      hs <- 2 * outer(dres.dp, dres.dp, "*")/sigma^2
      colnames(hs) <- rownames(hs) <- names(pouter)
    }
    out <- objlist(value = val, gradient = gr, hessian = hs)
    attr(out, controls$attr.name) <- out$value
    attr(out, "prediction") <- pred
    attr(out, "env") <- env
    return(out)
  }
  class(myfn) <- c("objfn", "fn")
  attr(myfn, "conditions") <- condition
  attr(myfn, "parameters") <- value[1]
  return(myfn)
}

# -------------------------------------------------------------------------#
# Other helper functions ----
# -------------------------------------------------------------------------#

#' Determine which parameters need sensitivities
#'
#' @param est.grid est.grid = data.table(condition, par1,...,parN)
#' @param trafo symbolic base trafo
#' @param reactions Object of class [eqnlist()]
#'
#' @return character, to go into the estimate-argument of odemodel
#' @export
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
getParametersToEstimate <- function(est.grid, trafo, reactions) {
  egNames <- names(est.grid)[-1]
  reg <- paste0("(", paste0(egNames, collapse = "|"), ")")
  trNames <- names(trafo)[grep(reg, trafo)]
  odeNames <- intersect(getParameters(reactions), trNames)
  odeNames
}


#' Title
#'
#' @param x 
#' @param ... 
#'
#' @return
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @export
#' @md
getParameters.data.table <- function(x,...) {
  unique(unlist(x[,!c("condition", "ID")], use.names = FALSE))
}

#' Get Parameter mappings outer.estgrid - inner.estgrid
#'
#' @param x est.grid
#'
#' @return named character
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#' @export
#' 
#' @examples 
#' est.grid <- data.frame(ID = 1:2,
#'                        condition = c("A", "B"),
#'                        k1 = c("k1_A", "k1_B"),
#'                        k2 = c("k2_A", "k2_B"),
#'                        k3 = c("k3", NA),
#'                        k4 = c("k4", "k4"),
#'                        stringsAsFactors = FALSE)
#' getEstGridParameterMapping(est.grid)
getEstGridParameterMapping <- function(x) {
  nm <- setdiff(names(x),c("condition", "ID"))
  do.call(c, lapply(nm, function(n) setNames(unique(x[[n]]), rep(n, length(unique(x[[n]]))))))
}


#' Run some checks on the prediction in normL2_indiv
#' 
#' If prediction is NA for observables which are not observed in a condition, they don't matter.
#' In this case, replace NA by 0, such that the error model can be evaluated.
#' 
#' @param prediction prediction for condition cn
#' @param data datalist
#' @param cn condition name for which the prediction was made
#' @param FLAGNaNInfwarnings print warnings?
#'
#' @return the prediction with harmless NA's replaced by 0
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
check_and_sanitize_prediction <- function(prediction, data, cn, FLAGNaNInfwarnings) {
  if (any(is.na(prediction))){
    whichcols <- unique(which(is.na(prediction), arr.ind = TRUE)[,2])
    nm <- colnames(prediction)[whichcols]
    
    if (length(intersect(data[[cn]]$name, nm)))
      stop("Prediction is.na for observables present in data in condition ", cn, "\n",
           "The following observables are affected: ", paste0(intersect(data[[cn]]$name, nm), collapse = ", "))
    
    if (FLAGNaNInfwarnings)
      warning("NaN in condition ", cn , " for the following names: ", paste0(nm, collapse = ", "))
    prediction[is.na(prediction)] <- 0
    attr(prediction, "deriv")[is.infinite(attr(prediction, "deriv"))|is.na(attr(prediction, "deriv"))] <- 0
    attr(prediction, "sensitivities")[is.infinite(attr(prediction, "sensitivities"))|is.na(attr(prediction, "sensitivities"))] <- 0
  }
  if (any(is.infinite(prediction))){
    whichcols <- unique(which(is.infinite(prediction), arr.ind = TRUE)[,2])
    nm <- colnames(prediction)[whichcols]
    
    if (length(intersect(data[[cn]]$name, nm)))
      warning("Prediction is infinite for observables present in data in condition ", cn, "\n",
              "The following observables are affected: ", paste0(intersect(data[[cn]]$name, nm), collapse = ", "),
              "These values are set to zero")
    
    if (FLAGNaNInfwarnings)
      warning("Inf in condition ", cn , " for the following names: ", paste0(nm, collapse = ", "))
    
    prediction[is.infinite(prediction)] <- 0
    attr(prediction, "deriv")[is.infinite(attr(prediction, "deriv"))|is.na(attr(prediction, "deriv"))] <- 0
    attr(prediction, "sensitivities")[is.infinite(attr(prediction, "sensitivities"))|is.na(attr(prediction, "sensitivities"))] <- 0
  }
  prediction
}
