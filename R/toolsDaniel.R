#' Run an R expression in the background (only on UNIX)
#' 
#' @description Generate an R code of the expression that is copied via \code{scp}
#' to any machine (ssh-key needed). Then collect the results.
#' @param expr An R expression
#' @param filename Character, defining the filename of the temporary file
#' @param machine Character, e.g. \code{"localhost"} or \code{"knecht1.fdm.uni-freiburg.de"}
#' @param input Character vector, the objects in the workspace that are stored
#' into an R data file and copied to the remove machine.
#' @param output Character vector, the objects in the workspace after evaluation of the
#' code that are stored into a results file. If \code{NULL}, just objects with names
#' different from elements in \code{input} are stored.
#' @return List of functions \code{get()} and \code{purge()}. \code{get()} copies the result file
#' to the working directory and loads it into the workspace. \code{purge()} deletes the temporary folder
#' from the working directory and the remote machine.
#' @export
#' @examples
#' 
#' out <- runbg(expression(M <- matrix(9:1, 3, 3)), 
#'              filename = "testrunbg", 
#'              machine = "localhost", 
#'              output = c("M"))
#' out$get()
#' print(M)
#' out$purge()
runbg <- function(expr, filename = "tmp", machine = "localhost", input = ls(), output = NULL) {
  
  # Save current workspace
  save(list = input, file = paste0(filename, ".RData"))
  
  # Get loaded packages
  pack <- sapply(strsplit(search(), "package:", fixed = TRUE), function(v) v[2])
  pack <- pack[!is.na(pack)]
  pack <- paste(paste0("library(", pack, ")"), collapse = "\n")
  
  # Define outputs
  if(is.null(output))
    output <- "setdiff(.newobjects, .oldobjects)"
  else
    output <- paste0("c('", paste(output, collapse = "', '"), "')")
  
  # Write program into character
  program <- paste(
    pack,
    paste0("setwd('~/", filename, "_folder')"),
    paste0("load('", filename, ".RData')"),
    ".oldobjects <- ls()",
    as.character(expr),
    ".newobjects <- ls()",
    
    paste0("save(list =", output ,", file = '", filename, "_result.RData')"),
    sep = "\n"
  )
  
  # Write program code into file
  cat(program, file = paste0(filename, ".R"))
  
  # Copy files to temporal folder
  system(paste0("ssh ", machine, " mkdir ", filename, "_folder/"))
  system(paste0("scp ", getwd(), "/", filename, ".R* ", machine, ":", filename, "_folder/"))
  system(paste0("scp ", getwd(), "/*.so ", machine, ":", filename, "_folder/"))
  
  # Run in background
  system(paste0("ssh ", machine, " R CMD BATCH ", filename, "_folder/", filename, ".R --vanilla"), intern = FALSE, wait = FALSE)
  
  out <- structure(vector("list", 2), names = c("get", "purge"))
  
  out[[1]] <- function() {
    
    system(paste0("scp ", machine, ":", filename, "_folder/", filename, "_result.RData ./"))
    load(file = paste0(filename, "_result.RData"), envir = .GlobalEnv, verbose = TRUE)
    
  }
  
  out[[2]] <- function() {
    
    system(paste0("ssh ", machine, " rm -r ", filename, "_folder"))
    system(paste0("rm ", filename, ".R*"))
    
  }
  
  return(out)
  
}





#' Generate the model objects for use in Xs (models with sensitivities)
#' 
#' @param f Named character vector with the ODE
#' @param forcings Character vector with the names of the forcings
#' @param fixed Character vector with the names of parameters (initial values and dynamic) for which
#' no sensitivities are required (will speed up the integration).
#' @param modelname Character, the name of the C file being generated.
#' @param ... Further arguments being passed to funC.
#' @return list with \code{func} (ODE object) and \code{extended} (ODE+Sensitivities object)
#' @export
#' @import cOde
generateModel <- function(f, forcings=NULL, fixed=NULL, modelname = "f", ...) {
  
  modelname_s <- paste0(modelname, "_s")
  
  func <- cOde::funC(f, forcings = forcings, modelname = modelname , ...)
  s <- sensitivitiesSymb(f, 
                         states = setdiff(attr(func, "variables"), fixed), 
                         parameters = setdiff(attr(func, "parameters"), fixed), 
                         inputs = forcings,
                         reduce = TRUE)
  fs <- c(f, s)
  outputs <- attr(s, "outputs")
  extended <- cOde::funC(fs, forcings = forcings, outputs = outputs, modelname = modelname_s, ...)
  
  list(func = func, extended = extended)
  
}


#' Generate the model objects for use in Xv (models with input estimation)
#' 
#' @param f Named character vector with the ODE
#' @param observed Character vector of the observed states (subset of \code{names(f)})
#' @param inputs Character vector of the input function to be estimated as part of the algorithm
#' @param forcings Character vector of forcings (the input functions that are NOT estimated
#' as part of the algorithm but still occur as driving force of the ODE)
#' @param ... Further arguments being passed to funC.
#' @return list with \code{func_au} (Combined ODE for states and adjoint sensitivities) and 
#' \code{func_l} (ODE of the log-likelihood and its gradient)
generateModelIE <- function(f, observed, inputs, forcings, scale=1, modelname = "f", ...) {
  
  
  nf <- length(f)
  ni <- length(inputs)
  
  fparse <- getParseData(parse(text=f), keep.source = TRUE)
  variables <- names(f)
  symbols <- unique(fparse$text[fparse$token == "SYMBOL"])
  forcings.t <- paste(c(forcings, inputs), "t", sep=".")
  parameters <- symbols[!symbols%in%c(variables, c(forcings, inputs), forcings.t)]
    
  
  ## Adjoint equtions + Input
  au <- adjointSymb(f, observed, inputs, parameters)
  forcings_au <- c(attr(au, "forcings"), names(f), c(forcings, inputs))
  au[1:length(au)] <- replaceSymbols(names(au), paste(scale, names(au), sep="*"), au)
  au[1:length(au)] <- paste0("(", au[1:length(au)], ")/", scale)
  attr(au, "inputs") <- replaceSymbols(names(au), paste(scale, names(au), sep="*"), attr(au, "inputs"))
    
  
  ## Input estimation equations
  fa <- c(f, au)
  fa <- replaceSymbols(inputs, attr(au, "inputs"), fa)
  boundary <- data.frame(name = names(fa),
                         yini = c(rep(1, nf), rep(NA, nf)),
                         yend = c(rep(NA, nf), rep(0, nf)))
  forcings_fa <- c(attr(au, "forcings"), forcings)
  func_fa <- cOde::funC(fa, forcings_fa, jacobian=TRUE, boundary=boundary, modelname = modelname, fcontrol = "einspline", ...)
  attr(func_fa, "inputs") <- attr(au, "inputs")
  
  ## Log-likelihood
  l <- c(attr(au, "chi"), attr(au, "grad"))
  forcings_l <- c(names(au), attr(au, "forcings"), names(f), c(forcings, inputs))
  func_l <- cOde::funC(l, forcings_l, modelname = paste0(modelname, "_l"), fcontrol = "einspline", ...)
  
  list(func_fa = func_fa, func_l = func_l)
  
  
}

#' Return some useful forcing functions as strings
#' 
#' @param type Which function to be returned
#' @param parameters Named vector, character or numeric. Replace parameters by the corresponding valus
#' in \code{parameters}.
#' @return String with the function
#' @export
forcingsSymb <- function(type =c("Gauss", "Fermi", "1-Fermi", "MM", "Signal"), parameters = NULL) {
  
  type <- match.arg(type)
  fun <- switch(type,
                "Gauss"   = "(scale*exp(-(time-mu)^2/(2*tau^2))/(tau*2.506628))",
                "Fermi"   = "(scale/(exp((time-mu)/tau)+1))",
                "1-Fermi" = "(scale*exp((time-mu)/tau)/(exp((time-mu)/tau)+1))",
                "MM"      = "(slope*time/(1 + slope*time/vmax))",
                "Signal"  = "max1*max2*(1-exp(-time/tau1))*exp(-time*tau2)"
  )
  
  if(!is.null(parameters)) {
    fun <- replaceSymbols(names(parameters), parameters, fun)
  }
  
  return(fun)
  
}

prepareFluxReduction <- function(f) {
  
  descr <- attr(f, "description")
  rates <- attr(f, "rates")
  S     <- attr(f, "SMatrix")
  
  fluxPars <- paste0("fluxPar_", 1:length(rates))
  rates <- paste0(fluxPars, "*(", rates, ")")
  data <- cbind(Description = descr, Rate = rates, as.data.frame(S))
  
  
  f <- generateEquations(data)
  
  fluxeq <- getFluxEquations(f)
  
  fluxEval <- funC.algebraic(fluxeq$fluxVector)
  
  attr(f, "fluxVector") <- fluxeq$fluxVector
  attr(f, "fluxPars") <- fluxPars
  attr(f, "fluxEval") <- fluxEval
  
  return(f)
  
  
}

getFluxEquations <- function(f) {
  
  fluxes.data <- with(attributes(f), {
    
    states <- colnames(SMatrix)
    
    fluxes <- do.call(rbind, lapply(1:length(states), function(i) {
      contribution <- which(!is.na(SMatrix[,i] ))
      data.frame(rate = rates[contribution], state = states[i], sign = SMatrix[contribution, i])
    }))
    rownames(fluxes) <- NULL
    
    return(fluxes)
    
  })
  
  fluxes.vector <- with(as.list(fluxes.data), paste0(sign, "*(", rate, ")"))
  names(fluxes.vector) <- paste(fluxes.data$state, fluxes.data$rate, sep=".")
  fluxes.vector.extended <- c(time = "time", fluxes.vector)
  
  states <- attr(f, "species")
  fluxes.unique <-as.character(unique(fluxes.data$rate))
  symbols <- getSymbols(fluxes.unique, exclude = states)
  nSymbols <- length(symbols)
  
  linears <- lapply(1:length(fluxes.unique), function(i) {
    
    isLinear <- unlist(lapply(symbols, function(mysymbol) {
      dflux <- paste(deparse(D(parse(text = fluxes.unique[i]), mysymbol)), collapse="")
      ddflux <- paste(deparse(D(D(parse(text = fluxes.unique[i]), mysymbol), mysymbol)), collapse="")
      return(ddflux == "0" & dflux != "0")
    }))
    
    return(symbols[isLinear])
    
    
  })
  names(linears) <- fluxes.unique
  
  
  return(list(fluxData = fluxes.data, fluxVector = fluxes.vector, linearContributors = linears))
  
}



getZeroFluxes <- function(out, rtol = .05, atol = 0) {
  
  ## out must have the colnames "time", "state.fluxEquation"
  ## states are not allowed to contain a "."-symbol
  
  states <- sapply(strsplit(colnames(out)[-1], ".", fixed=TRUE), function(v) v[1])
  fluxes <- sapply(strsplit(colnames(out)[-1], ".", fixed=TRUE), function(v) paste(v[-1], collapse=""))
  unique.states <- unique(states)
  unique.fluxes <- unique(fluxes)
  
  out.groups <- lapply(unique.states, function(s) {
    
    selected <- which(states == s)
    out.selected <- matrix(out[, selected + 1], nrow=dim(out)[1])
    colnames(out.selected) <- fluxes[selected]
    
    # Get L1 norm of fluxes
    abssum <- apply(abs(out.selected), 2, sum)
    abssum.extended <- rep(0, length(unique.fluxes))
    names(abssum.extended) <- unique.fluxes
    abssum.extended[names(abssum)] <- abssum
    
    # Normalize with respect to the L1 norm of the state derivative (sum of all fluxes)
    state.dot <- apply(out.selected, 1, sum)
    norm.state.dot <- sum(abs(state.dot))
    
    abssum.normed <- abssum/norm.state.dot
    abssum.normed.extended <- rep(0, length(unique.fluxes))
    names(abssum.normed.extended) <- unique.fluxes
    abssum.normed.extended[names(abssum.normed)] <- abssum.normed
    
    return(list(abssum.extended, abssum.normed.extended))
    
    
  })
  
  out.groups.abs <- do.call(rbind, lapply(out.groups, function(g) g[[1]]))
  rownames(out.groups.abs) <- unique.states
  out.groups.rel <- do.call(rbind, lapply(out.groups, function(g) g[[2]]))
  rownames(out.groups.rel) <- unique.states
  
  zero.fluxes.abs <- unique.fluxes[apply(out.groups.abs, 2, function(v) all(v < atol))]
  zero.fluxes.rel <- unique.fluxes[apply(out.groups.rel, 2, function(v) all(v < rtol))]
  zero.fluxes <- c(zero.fluxes.abs, zero.fluxes.rel)
  non.zero.fluxes <- unique.fluxes[!unique.fluxes%in%zero.fluxes]
  
  
  #non.zero.fluxes.rel <- unique.fluxes[apply(out.groups, 2, function(v) any(v > rtol))]
  
  zero.parameters <- unlist(lapply(strsplit(zero.fluxes, "*", fixed=TRUE), function(v) v[1]))
  
  return(list(fluxes.abs = out.groups.abs, 
              fluxes.rel = out.groups.rel, 
              fluxes.zero = zero.fluxes, 
              fluxes.nonzero = non.zero.fluxes, 
              parameters.zero = zero.parameters))
  
}


normalizeData <- function(data) {
  
  names <- unique(data$name)
  data.normalized <- do.call(rbind, lapply(names, function(n) {
    sub <- data[data$name == n,]
    mean.value <- mean(sub$value)
    sub$value <- sub$value/mean.value
    sub$sigma <- sub$sigma/abs(mean.value)
    return(sub)
  }))
  return(data.normalized)
  
}



#' Soft L2 constraint on parameters
#' 
#' @param p Namec numeric, the parameter value
#' @param mu Named numeric, the prior values
#' @param sigma Named numeric of length of mu or numeric of length one.
#' @param fixed Named numeric with fixed parameter values (contribute to the prior value
#' but not to gradient and Hessian)
#' @return List of class \code{obj}, i.e. objective value, gradient and Hessian as list.
#' @seealso \link{wrss}
#' @details Computes the constraint value 
#' \deqn{\frac{1}{2}\left(\frac{p-\mu}{\sigma}\right)^2}{0.5*(p-mu)^2/sigma^2}
#' and its derivatives with respect to p.
#' @examples
#' p <- c(A = 1, B = 2, C = 3)
#' mu <- c(A = 0, B = 0)
#' sigma <- c(A = 0.1, B = 1)
#' constraintL2(p, mu, sigma)
#' @export
constraintL2 <- function(p, mu, sigma = 1, fixed=NULL) {

  ## Augment sigma if length = 1
  if(length(sigma) == 1) 
    sigma <- structure(rep(sigma, length(mu)), names = names(mu)) 
  
  ## Extract contribution of fixed pars and delete names for calculation of gr and hs  
  par.fixed <- intersect(names(mu), names(fixed))
  sumOfFixed <- 0
  if(!is.null(par.fixed)) sumOfFixed <- sum(0.5*((fixed[par.fixed] - mu[par.fixed])/sigma[par.fixed])^2)
  
                         
  # Compute prior value and derivatives
  par <- intersect(names(mu), names(p))
    
  val <- sum((0.5*((p[par]-mu[par])/sigma[par])^2)) + sumOfFixed
  gr <- rep(0, length(p)); names(gr) <- names(p)
  gr[par] <- ((p[par]-mu[par])/(sigma[par]^2))
  
  hs <- matrix(0, length(p), length(p), dimnames = list(names(p), names(p)))
  diag(hs)[par] <- 1/sigma[par]^2
  
  dP <- attr(p, "deriv")
  if(!is.null(dP)) {
    gr <- as.vector(gr%*%dP); names(gr) <- colnames(dP)
    hs <- t(dP)%*%hs%*%dP; colnames(hs) <- colnames(dP); rownames(hs) <- colnames(dP)
  }
  
  out <- list(value=val,gradient=gr,hessian=hs)
  class(out) <- c("obj", "list")
  
  return(out)
  
}




#' L2 objective function for validation data point
#' 
#' @param p Namec numeric, the parameter values
#' @param prediction Matrix with first column "time" and one column per predicted state. Can have
#' an attribute \code{deriv}, the matrix of sensitivities. If present, derivatives of the objective
#' function with respect to the parameters are returned.
#' @param mu Named character of length one. Has the structure \code{mu = c(parname = statename)}, where
#' \code{statename} is one of the column names of \code{prediction} and \code{parname} is one of the
#' names of \code{p}, allowing to treat the validation data point as a parameter.
#' @param time Numeric of length one. An existing time point in \code{prediction}.
#' @param sigma Numeric of length one. The uncertainty assumed for the validation data point.
#' @param fixed Named numeric with fixed parameter values (contribute to the prior value
#' but not to gradient and Hessian)
#' @return List of class \code{obj}, i.e. objective value, gradient and Hessian as list.
#' @seealso \link{wrss}, \link{constraintL2}
#' @details Computes the constraint value 
#' \deqn{\left(\frac{x(t)-\mu}{\sigma}\right)^2}{(pred-p[names(mu)])^2/sigma^2}
#' and its derivatives with respect to p.
#' @examples
#' \dontrun{
#' prediction <- matrix(c(0, 1), nrow = 1, dimnames = list(NULL, c("time", "A")))
#' derivs <- matrix(c(0, 1, 0.1), nrow = 1, dimnames = list(NULL, c("time", "A.A", "A.k1")))
#' attr(prediction, "deriv") <- derivs
#' p0 <- c(A = 1, k1 = 2)
#' mu <- c(newpoint = "A")
#' timepoint <- 0
#' 
#' datapointL2(p = c(p, newpoint = 2), prediction, mu, timepoint)
#' datapointL2(p = c(p, newpoint = 1), prediction, mu, timepoint)
#' datapointL2(p = c(p, newpoint = 0), prediction, mu, timepoint)
#' }
#' @export
datapointL2 <- function(p, prediction, mu, time = 0, sigma = 1, fixed = NULL) {
  
  
  # Only one data point is allowed
  mu <- mu[1]; time <- time[1]; sigma <- sigma[1]
  
  # Divide parameter into data point and rest
  datapar <- setdiff(names(mu), names(fixed))
  parapar <- setdiff(names(p), c(datapar, names(fixed)))
  
  
  # Get predictions and derivatives at time point
  time.index <- which(prediction[,"time"] == time)
  withDeriv <- !is.null(attr(prediction, "deriv"))
  pred <- prediction[time.index, ]
  deriv <- NULL
  if(withDeriv)
    deriv <- attr(prediction, "deriv")[time.index, ]
  
  # Reduce to name = mu
  pred <- pred[mu]
  if(withDeriv) {
    mu.para <- intersect(paste(mu, parapar, sep = "."), names(deriv))
    deriv <- deriv[mu.para]
  }
  
  # Compute prior value and derivatives
  res <- as.numeric(pred - c(fixed, p)[names(mu)])
  val <- as.numeric((res/sigma)^2)
  gr <- NULL
  hs <- NULL
  
  if(withDeriv) {
    dres.dp <- structure(rep(0, length(p)), names = names(p))
    if(length(parapar) > 0) dres.dp[parapar] <- as.numeric(deriv)
    if(length(datapar) > 0) dres.dp[datapar] <- -1
    gr <- 2*res*dres.dp/sigma^2
    hs <- 2*outer(dres.dp, dres.dp, "*")/sigma^2; colnames(hs) <- rownames(hs) <- names(p)
  }
  
  out <- list(value=val,gradient=gr,hessian=hs)
  class(out) <- c("obj", "list")
  
  return(out)
  
}

#' L2 objective function for prior value
#' 
#' @description As a prior function, it returns derivatives with respect to
#' the penalty parameter in addition to parameter derivatives.
#' 
#' @param p Namec numeric, the parameter value
#' @param mu Named numeric, the prior values
#' @param lambda Character of length one. The name of the penalty paramter in \code{p}.
#' @param fixed Named numeric with fixed parameter values (contribute to the prior value
#' but not to gradient and Hessian)
#' @return List of class \code{obj}, i.e. objective value, gradient and Hessian as list.
#' @seealso \link{wrss}, \link{constraintExp2}
#' @details Computes the constraint value 
#' \deqn{e^{\lambda} \| p-\mu \|^2}{exp(lambda)*sum((p-mu)^2)}
#' and its derivatives with respect to p and lambda.
#' @examples
#' p <- c(A = 1, B = 2, C = 3, lambda = 0)
#' mu <- c(A = 0, B = 0)
#' priorL2(p, mu, lambda = "lambda")
#' @export
priorL2 <- function(p, mu, lambda = "lambda", fixed = NULL) {
  
  ## Extract contribution of fixed pars and delete names for calculation of gr and hs  
  par.fixed <- intersect(names(mu), names(fixed))
  sumOfFixed <- 0
  if(!is.null(par.fixed)) sumOfFixed <- sum(exp(c(fixed, p)[lambda])*(fixed[par.fixed] - mu[par.fixed])^2)
  
  
  # Compute prior value and derivatives
  par <- intersect(names(mu), names(p))
  par0 <- setdiff(par, lambda)
  
  val <- sum(exp(c(fixed, p)[lambda]) * (p[par]-mu[par])^2) + sumOfFixed
  gr <- rep(0, length(p)); names(gr) <- names(p)
  gr[par] <- 2*exp(c(fixed, p)[lambda])*(p[par]-mu[par])
  if(lambda %in% names(p)) {
    gr[lambda] <- sum(exp(c(fixed, p)[lambda]) * (p[par0]-mu[par0])^2) + sum(exp(c(fixed, p)[lambda]) * (fixed[par.fixed] - mu[par.fixed])^2)
  }
  
  hs <- matrix(0, length(p), length(p), dimnames = list(names(p), names(p)))
  diag(hs)[par] <- 2*exp(c(fixed, p)[lambda])
  if(lambda %in% names(p)) {
    hs[lambda, lambda] <- gr[lambda] 
    hs[lambda, par0] <- hs[par0, lambda] <- gr[par0]
  }
  
  dP <- attr(p, "deriv")
  if(!is.null(dP)) {
    gr <- as.vector(gr%*%dP); names(gr) <- colnames(dP)
    hs <- t(dP)%*%hs%*%dP; colnames(hs) <- colnames(dP); rownames(hs) <- colnames(dP)
  }
  
  out <- list(value=val,gradient=gr,hessian=hs)
  class(out) <- c("obj", "list")
  
  return(out)
  
  
}



#' Add two lists element by element
#' 
#' @param out1 List of numerics or matrices
#' @param out2 List with the same structure as out1 (there will be no warning when mismatching)
#' @details If out1 has names, out2 is assumed to share these names. Each element of the list out1
#' is inspected. If it has a \code{names} attributed, it is used to do a matching between out1 and out2.
#' The same holds for the attributed \code{dimnames}. In all other cases, the "+" operator is applied
#' the corresponding elements of out1 and out2 as they are.
#' @return List of length of out1. 
#' @aliases summation
#' @export "+.obj"
#' @export
"+.obj" <- function(out1, out2) {
  
  allnames <- c(names(out1), names(out2))
  what <- allnames[duplicated(allnames)]
  what.names <- what
  if(is.null(what)) {
    what <- 1:min(c(length(out1), length(out2)))
    what.names <- NULL
  }
  
  out12 <- lapply(what, function(w) {
    sub1 <- out1[[w]]
    sub2 <- out2[[w]]
    n <- names(sub1)
    dn <- dimnames(sub1)
    if(!is.null(n) && !is.null(sub1) %% !is.null(sub2)) {
      #print("case1: sum of vectors")
      sub1[n] + sub2[n]
    } else if(!is.null(dn) && !is.null(sub1) && !is.null(sub2)) {
      #print("case2: sum of matrices")
      matrix(sub1[dn[[1]], dn[[2]]] + sub2[dn[[1]], dn[[2]]], 
             length(dn[[1]]), length(dn[[2]]), dimnames = list(dn[[1]], dn[[2]]))
    } else if(!is.null(sub1) && !is.null(sub2)) {
      #print("case3: sum of scalars")
      sub1 + sub2
    } else {
      #print("case4")
      NULL
    }
  })
  names(out12) <- what.names
  
  class(out12) <- c("obj", "list")
  
  return(out12)
}



#' @export
mssample <- function(center, samplefun = "rnorm", fits = 20, ...) {
  
  sample.matrix <- do.call(rbind, lapply(1:fits, function(i) 
    do.call(samplefun, c(list(n = length(center), mean = center), list(...)))))
  
  colnames(sample.matrix) <- names(center)
  
  return(sample.matrix)

}

#' @export
ggopen <- function(plot = last_plot(), command = "xdg-open", ...) {
  filename <- tempfile(pattern = "Rplot", fileext = ".pdf")
  ggsave(filename = filename, plot = plot, ...)
  system(command = paste(command, filename))
}
