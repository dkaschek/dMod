#' @export
modelNLME <- function(prdfn, covtable = NULL, cores = 1) {
  
  
  covnames <- names(covtable)
  parnames <- getParameters(prdfn)
  
  
  model <- function(time, name, ...) {
    
    pars <- as.data.frame(c(list(time, name), list(...)))
    
    names(pars) <- c("time", "name", parnames, covnames)
    
    id <- cumsum(Reduce("|", lapply(pars[-(1:2)], function(x) !duplicated(x))))
    pars <- split(pars, id)
    
    output <- parallel::mclapply(pars, function(sub) {
      timesD <- unique(sub$time)
      parsD <- unlist(sub[1, parnames])
      condition <- paste(unlist(sub[1, covnames]), collapse = "_")
      
      prediction <- prdfn(timesD, parsD, conditions = condition)[[1]]
      template <- data.frame(name = sub$name, time = sub$time, value = 0, sigma = 1)
      
      myres <- res(template, prediction)
      
      output <- myres$prediction
      deriv <- as.matrix(attr(myres, "deriv")[, -(1:2)])
      
      list(output, deriv)
    }, mc.cores = cores)
    
    gradient <- do.call(rbind, lapply(output, function(x) x[[2]]))
    output <- unlist(lapply(output, function(x) x[[1]]))
    
    attr(output, "gradient") <- gradient
    
    return(output)
    
  }
  
  return(model)
  
}


#' @export
modelSAEMIX <- function(prdfn, cores = 1) {
  
  parnames <- getParameters(prdfn)
  
  model <- function(psi, id, xidep) {
    
    pars <- split(as.data.frame(psi[id, ]), id)
    
    output <- do.call(c, parallel::mclapply(1:nrow(psi), function(i) {
      
      parsD <- unlist(psi[i,])
      names(parsD) <- parnames
      timesD <- as.numeric(xidep[id == i, 1])
      namesD <- as.character(xidep[id == i, 2])
      
      prediction <- prdfn(timesD, parsD, deriv = FALSE)[[1]]
      template <- data.frame(name = namesD, time = timesD, value = 0, sigma = 1)
      
      myres <- res(template, prediction)
      
      return(myres$prediction)
      
      
    }, mc.cores = cores))
    
    return(output)
    
    
  }
  
  return(model)
  
  
}