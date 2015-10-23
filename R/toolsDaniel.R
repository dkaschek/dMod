


#' Run an R expression in the background (only on UNIX)
#' 
#' @description Generate an R code of the expression that is copied via \code{scp}
#' to any machine (ssh-key needed). Then collect the results.
#' @param ... Some R code
#' @param filename Character, defining the filename of the temporary file. Random
#' file name ist chosen if NULL.
#' @param machine Character, e.g. \code{"localhost"} or \code{"knecht1.fdm.uni-freiburg.de"}
#' @param input Character vector, the objects in the workspace that are stored
#' into an R data file and copied to the remove machine.
#' @param compile Logical. If \code{TRUE}, C files are copied and compiled on the remote machine.
#' Otherwise, the .so files are copied.
#' @return List of functions \code{check}, \code{get()} and \code{purge()}. 
#' \code{check()} checks, if the result is ready.
#' \code{get()} copies the result file
#' to the working directory and loads it into the workspace. 
#' \code{purge()} deletes the temporary folder
#' from the working directory and the remote machine.
#' @export
#' @examples
#' 
#' out <- runbg(machine = "localhost", {
#'          M <- matrix(9:1, 3, 3)}
#'        )
#' out$get()
#' print(.runbgOutput)
#' out$purge()
runbg <- function(..., machine = "localhost", filename = NULL, input = ls(.GlobalEnv), compile = FALSE) {
  
  
  expr <- as.expression(substitute(...))
  
  # Set file name
  if (is.null(filename))
    filename <- paste0("tmp_", paste(sample(c(0:9, letters), 5, replace = TRUE), collapse = ""))
  
  # Save current workspace
  save(list = input, file = paste0(filename, ".RData"))
  
  # Get loaded packages
  pack <- sapply(strsplit(search(), "package:", fixed = TRUE), function(v) v[2])
  pack <- pack[!is.na(pack)]
  pack <- paste(paste0("library(", pack, ")"), collapse = "\n")
  
  # Define outputs
  output <- ".runbgOutput"
  # if(is.null(output))
  #   output <- "setdiff(.newobjects, .oldobjects)"
  # else
  #   output <- paste0("c('", paste(output, collapse = "', '"), "')")
 
  compile.line <- NULL
  if(compile)
    compile.line <- "cfiles <- list.files(pattern = '.c$'); for(cf in cfiles) system(paste('R CMD SHLIB', cf))"
   
  # Write program into character
  program <- paste(
    pack,
    paste0("setwd('~/", filename, "_folder')"),
    compile.line,
    paste0("load('", filename, ".RData')"),
    #".oldobjects <- ls()",
    paste0(".runbgOutput <- ", as.character(expr)),
    #".newobjects <- ls()",
    
    paste0("save(", output ,", file = '", filename, "_result.RData')"),
    sep = "\n"
  )
  
  # Write program code into file
  cat(program, file = paste0(filename, ".R"))
  
  # Copy files to temporal folder
  system(paste0("ssh ", machine, " mkdir ", filename, "_folder/"))
  system(paste0("scp ", getwd(), "/", filename, ".R* ", machine, ":", filename, "_folder/"))
  if(compile) {
    system(paste0("scp ", getwd(), "/*.c ", machine, ":", filename, "_folder/"))
  } else {
    system(paste0("scp ", getwd(), "/*.so ", machine, ":", filename, "_folder/"))
  }
  
  # Run in background
  system(paste0("ssh ", machine, " R CMD BATCH ", filename, "_folder/", filename, ".R --vanilla"), intern = FALSE, wait = FALSE)
  
  out <- structure(vector("list", 3), names = c("check", "get", "purge"))
  
  out[[1]] <- function() {
    
    check.out <- suppressWarnings(
      system(paste0("ssh ", machine, " ls ", filename, "_folder/ | grep -x ", filename, "_result.RData"), 
             intern = TRUE))
    
    if(length(check.out) > 0) cat("Result is ready!\n") else cat("Not ready!\n")
    
  }
  
  out[[2]] <- function() {
    
    system(paste0("scp ", machine, ":", filename, "_folder/", filename, "_result.RData ./"))
    load(file = paste0(filename, "_result.RData"), envir = .GlobalEnv, verbose = TRUE)
    
  }
  
  out[[3]] <- function() {
    
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
