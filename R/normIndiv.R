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
  
  pars        <- unclass_parvec(pars)
  fixed       <- unclass_parvec(fixed)
  
  pars_outer  <- pars
  fixed_outer <- fixed
  
  # Match parameters to est.grid: Need to consider supplied "fixed" as well!s
  pars <- c(pars, fixed)
  parnames  <- unlist(est.grid[est.grid$ID == ID, setdiff(names(est.grid), c("ID", "condition"))])
  parnames <- parnames[!is.na(parnames)]
  pars <- setNames(pars[parnames], names(parnames))
  
  # Get Parameters from fix.grid
  fixed <- unlist(fix.grid[fix.grid$ID == ID, setdiff(names(fix.grid), c("ID", "condition"))])
  # remove NAs
  fixed <- fixed[!is.na(fixed)]
  
  # Sort supplied "fixed" parameters back to fixed
  fixed <- c(fixed, pars[parnames %in% names(fixed_outer)])
  pars <- pars[!parnames %in% names(fixed_outer)]
  parnames <- parnames[!parnames %in% names(fixed_outer)]
  
  return(list(pars = unlist(pars), fixed = unlist(fixed), parnames = parnames))
}


#' Title
#'
#' @param pred0 prd0(times,pars)[[1]]
#' @param pars pars
#' @param est.grid est.grid
#' @param cn name of condition
#'
#' @return pred0 with updated deriv argument
#' @export
#'
#' @examples
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
  attr(pred0, "deriv") <- der
  pred0
}


#' Title
#'
#' @param x 
#' @param ... 
#'
#' @return
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#'
#' @examples
getParameters.data.table <- function(x,...) {
  unique(unlist(x[,!c("condition", "ID")], use.names = FALSE))
}


add_pars_to_grid <- function(pars, grid, FLAGoverwrite) {
  pars <- as.data.table(t(pars))
  if (FLAGoverwrite) {
    # Replace
    grid <- data.table(grid[,!names(..pars_chr)], pars)
  } else {
    # Append setdiff
    parnm <- setdiff(names(pars), names(grid)) 
    if (length(parnm))
      grid <- data.table(grid, pars[,..parnm])
  }
  grid
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
add_pars_to_grids <- function(pars, gridlist, FLAGoverwrite = FALSE) {
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
  if (sum( is_num)) fix.grid <- data.table(fix.grid, as.data.table(t(pars[ is_num])))
  
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
#'
#' @examples
PRD_indiv <- function(prd0, est.grid, fix.grid) {
  
  # Title
  #
  # @param times 
  # @param pars 
  # @param fixed 
  # @param deriv 
  # @param conditions 
  # @param FLAGbrowser 0: Don't debug, 1: Debug when there is an error, 2: always debug
  # @param FLAGverbose 
  # @param FLAGrenameDerivPars Needed for datapointL2_indiv, where I need derivs wrt the outer parameters
  #
  # @return
  # @export
  #
  # @examples
  prd <- function(times, pars, fixed = NULL, deriv = FALSE, conditions = est.grid$condition, 
                  FLAGbrowser = 0, 
                  FLAGverbose = FALSE,
                  FLAGrenameDerivPars = FALSE
  ) {
    out <- lapply(setNames(nm = conditions), function(cn) {
      if (FLAGbrowser == 2) browser()
      ID <- est.grid$ID[est.grid$condition == cn]
      if (FLAGverbose) cat(ID, cn, "\n", sep = " ---- ")
      
      dummy <- make_pars(pars, fixed, est.grid, fix.grid, ID)
      pars_ <- dummy$pars
      fixed_ <- dummy$fixed
      
      if (length(setdiff(getParameters(prd0), names(c(pars_, fixed_)))))
        stop("The following parameters are missing: ", paste0(setdiff(getParameters(prd0), names(c(pars_, fixed_))), collapse = ", "))
      
      pred0 <-try(prd0(times, pars_, fixed = fixed_, deriv = deriv, conditions = NULL)[[1]])
      if (inherits(pred0, "try-error") && FLAGbrowser == 1) {
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
#'
#' @examples
P_indiv <- function(p0, est.grid, fix.grid) {
  
  # @param FLAGbrowser 0: Don't debug, >= 1: debug
  prd <- function(pars, fixed = NULL, deriv = FALSE, conditions = est.grid$condition, 
                  FLAGbrowser = FALSE, 
                  FLAGverbose = FALSE) {
    out <- lapply(setNames(nm = conditions), function(cn) {
      if (FLAGbrowser) browser()
      ID <- est.grid$ID[est.grid$condition == cn]
      if (FLAGverbose) cat(ID, cn, "\n", sep = " ---- ")
      dummy <- make_pars(pars, fixed, est.grid, fix.grid, ID)
      pars_ <- dummy$pars
      fixed_ <- dummy$fixed
      
      if (length(setdiff(getParameters(prd0), names(c(pars_, fixed_)))))
        stop("The following parameters are missing: ", paste0(setdiff(getParameters(prd0), names(c(pars_, fixed_))), collapse = ", "))
      p0(pars_, fixed = fixed_, deriv = deriv)[[1]]
    })
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
#'
#' @importFrom parallel mclapply
normL2_indiv <- function (data, prd0, errmodel = NULL, est.grid, fix.grid, times = NULL, attr.name = "data", fixed.conditions = NULL) {
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
    calc_objval <- function(cn) {
      
      if (FLAGbrowser) browser()
      
      ID <- est.grid$ID[est.grid$condition == cn]
      if (FLAGverbose) cat(ID, cn, "\n", sep = " ---- ")
      dummy <- make_pars(pars, fixed, est.grid, fix.grid, ID)
      pars_ <- dummy$pars
      fixed_ <- dummy$fixed
      
      timesD <- controls$times
      attr.name <- controls$attr.name
      
      myderiv <- deriv
      if (!is.null(FixedConditions) && cn %in% FixedConditions) myderiv <- FALSE
      
      prediction <- try(prd0(times = timesD, pars = pars_, fixed = fixed_, deriv = myderiv))
      
      if (inherits(prediction, "try-error"))
        stop("Prediction failed in condition = ", cn, ", ID = ", ID, ".
             Try iterating p(pars), (x*p)(pars), ... to find the problem.")
      
      prediction <- prediction[[1]]
      
      # [] refactor: put the following stuff into own function catch_nonproblematicNanInfs(prediciton, data, cn, FLAGNaNInfWarnings)
      whichcols <- nm <- NULL
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
          stop("Prediction is infinite for observables present in data in condition ", cn, "\n",
               "The following observables are affected: ", paste0(intersect(data[[cn]]$name, nm), collapse = ", "))
        
        if (FLAGNaNInfwarnings)
          warning("Inf in condition ", cn , " for the following names: ", paste0(nm, collapse = ", "))
        
        prediction[is.infinite(prediction)] <- 0
        attr(prediction, "deriv")[is.infinite(attr(prediction, "deriv"))|is.na(attr(prediction, "deriv"))] <- 0
        attr(prediction, "sensitivities")[is.infinite(attr(prediction, "sensitivities"))|is.na(attr(prediction, "sensitivities"))] <- 0
      }
      
      err <- NULL
      if (any(is.na(data[[cn]]$sigma))) {
        err <- errmodel(out = prediction, pars = getParameters(prediction), conditions = cn, deriv=myderiv)
        mywrss <- nll(res(data[[cn]], prediction, err[[1]]), deriv = deriv, pars = pars)
      } else {
        mywrss <- nll(res(data[[cn]], prediction), deriv = deriv, pars = pars)
      }
      if (myderiv) {
        if (FLAGbrowser2) browser()
        mywrss$gradient <- mywrss$gradient[names(dummy$parnames)]
        names(mywrss$gradient) <- unname(dummy$parnames)
        
        mywrss$hessian <- mywrss$hessian[names(dummy$parnames),names(dummy$parnames)]
        dimnames(mywrss$hessian) <- list(unname(dummy$parnames), unname(dummy$parnames))
      }
      
      # [] catch conditions with NA value, don't include them in obj-calculation and print out warning
      return(mywrss)
    }
    
    if (simcores == 1)
      objlists <- lapply(setNames(nm = conditions), calc_objval)
    if (simcores > 1)
      objlists <- parallel::mclapply(setNames(nm = conditions), calc_objval, mc.cores = simcores)
    
    # Sum all objlists
    out <- objlist(value = 0, 
                   gradient = structure(rep(0, length(names(pars))), names = names(pars)), 
                   hessian = matrix(0, nrow = length(names(pars)), ncol = length(names(pars)), 
                                    dimnames = list(names(pars), names(pars))))
    out$value <- do.call(sum, lapply(objlists, function(.x) .x$value))
    if (deriv) {
      for (gr in lapply(objlists, function(.x) .x$gradient))
        out$gradient[names(gr)] <- out$gradient[names(gr)] + gr
      for (hs in lapply(objlists, function(.x) .x$hessian))
        out$hessian[rownames(hs), colnames(hs)] <- out$hessian[rownames(hs), colnames(hs)] + hs
    }
    
    # consider fixed: return only derivs wrt pouter
    out$gradient <- out$gradient[names(pars)]
    out$hessian <- out$hessian[names(pars), names(pars)]
    
    attr(out, controls$attr.name) <- out$value
    attr(out, "condition_obj") <- vapply(objlists, function(.x) .x$value, 1)
    attr(out, "AIC") <- out$value + length(pars) * 2
    attr(out, "BIC") <- out$value + length(pars) * log(nrow(as.data.frame(data)))
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
#'
#' @examples
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
#'
#' @examples
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

