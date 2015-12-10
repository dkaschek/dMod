#' Translate data into time-continuous data representations
#' 
#' @param data \code{data.frame} with name (factor), time (numeric), value (numeric) and sigma (numeric)
#' @param tau Numeric. The uncertainty of a single data point is smeared by a Gaussian with variance tau².
#' @param ngrid The number of time points returned.
#' @param type Character, the type of interpolation, i.e. "linear", "spline", "logspline" or "smooth".
#' @return a \code{data.frame} with columns "name", "time" and "value".
data2forc <- function(data, tau=NULL, ngrid = 1e3, type="linear") {
  
  obs <- unique(data$name)
  obsD <- paste0(obs, "D")
  weightsD <- paste0("weight", obs, "D")
  tmin <- min(data$time)
  tmax <- max(data$time)
  times <- seq(tmin, tmax, len=ngrid)
  
  # data points
  dataout <- do.call(rbind, lapply(obs, function(o) {
    
    subdata <- subset(data, name == o)
    if(type=="linear") {
      t <- subdata$time
      x <- subdata$value
      xoft <- approxfun(t, x)
      xvalues <- xoft(times)
      
      xofts <- splinefun(times, xvalues)
      
      out <- rbind(data.frame(name = paste0(o, "D"),time = times, value = xvalues),
                   data.frame(name = paste0(o, "DDot"),time = times, value = xofts(times, deriv=1)))
    }
    if(type=="spline") {
      t <- subdata$time
      x <- subdata$value
      xoft <- splinefun(t, x)
      out <- rbind(data.frame(name = paste0(o, "D"),time = times, value = xoft(times)),
                   data.frame(name = paste0(o, "DDot"), time = times, value = xoft(times, deriv=1)))
    }
    if(type=="logspline") {
      t <- subdata$time
      
      if(any(subdata$value <= 0)) {
        # no log
        x <- subdata$value 
      } else {
        # log
        subdata$value[subdata$value < 1e-5*max(subdata$value)] <- 1e-5*max(subdata$value)
        x <- log(subdata$value)
      }
      
      xoft <- splinefun(t, x)
      
      if(any(subdata$value <= 0)) {
        # no log
        out <- rbind(data.frame(name = paste0(o, "D"),time = times, value = xoft(times)),
                     data.frame(name = paste0(o, "DDot"), time = times, value = splinefun(times, xoft(times))(times, deriv=1)))
      } else {
        # log
        out <- rbind(data.frame(name = paste0(o, "D"),time = times, value = exp(xoft(times))),
                     data.frame(name = paste0(o, "DDot"), time = times, value = splinefun(times, exp(xoft(times)))(times, deriv=1)))
      }
      
      
    }
    if (type == "smooth") {
      t <- subdata$time
      x <- subdata$value
      xoft <- function(newt, ...) predict(smooth.spline(t, x, df = 5), newt, ...)$y
      out <- rbind(data.frame(name = paste0(o, "D"), time = times, 
                              value = xoft(times)), data.frame(name = paste0(o, 
                                                                             "DDot"), time = times, value = xoft(times, deriv = 1)))
      
    }
    
    return(out)
    
  }))
  
  # sigma value  
  
  
  
  if(is.null(tau)) {
    Deltat <- (tmax-tmin)/(dim(data)[1]/length(obs))
    tau <- Deltat/5
  }
  
  if(tau == 0) {
    weightout <- do.call(rbind, lapply(obs, function(o) {
      
      subdata <- subset(data, name == o & !is.na(value) & !is.na(sigma))
      
      t <- subdata$time
      x <- log(1/subdata$sigma^2)
      xoft <- splinefun(t, x)
      
      
      rbind(data.frame(name = paste0("weight", o, "D"), time = times, value = exp(xoft(times))),
            data.frame(name = paste0("weight", o, "DDot"), time=times, value = splinefun(times, exp(xoft(times)))(times, deriv=1)))
      
    }))
    
  } else {
    weightout <- do.call(rbind, lapply(obs, function(o) {
      
      subdata <- subset(data, name == o & !is.na(value) & !is.na(sigma))
      support <- subdata$time
      supfactor <- rep(1, length(support))
      if(any(support==tmin)) supfactor[support==tmin] <- 2
      if(any(support==tmax)) supfactor[support==tmax] <- 2
      sigma <- subdata$sigma
      weights <- Reduce("+", lapply(1:length(support), function(i) supfactor[i]*dnorm(times, support[i], tau) * (1/sigma[i]^2)))
      weightsdot <- splinefun(times, weights)(times, deriv=1)
      rbind(data.frame(name = paste0("weight", o, "D"), time = times, value = weights),
            data.frame(name = paste0("weight", o, "DDot"), time = times, value = weightsdot))
    }))
    
  }
  
  out <- rbind(dataout, weightout)
  
  return(out)
  
  
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


similarity <- function(fitlist, x, times, fixed, deriv = FALSE) {
  
  combinations <- combn(1:nrow(fitlist), 2)
  out <- do.call(rbind, lapply(1:ncol(combinations), function(j) {
    
    i1 <- combinations[1, j]
    i2 <- combinations[2, j]
    
    par1 <- unlist(fitlist[i1, -1])
    par2 <- unlist(fitlist[i2, -1])
    
    pred1 <- wide2long(x(times, par1, fixed = fixed, deriv = deriv))
    pred2 <- wide2long(x(times, par2, fixed = fixed, deriv = deriv))
    
    rss <- sum(((pred1$value - pred2$value)/(pred1$value + pred2$value + 1))^2)

    data.frame(i1 = c(i1, i2), i2 = c(i2, i1), rss = rss)    
    
  }))
  
  return(out)
  
}


