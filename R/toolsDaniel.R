


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
#' @param wait Wait until executed
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
runbg <- function(..., machine = "localhost", filename = NULL, input = ls(.GlobalEnv), compile = FALSE, wait = FALSE) {
  
  
  expr <- as.expression(substitute(...))
  
  # Set file name
  if (is.null(filename))
    filename <- paste0("tmp_", paste(sample(c(0:9, letters), 5, replace = TRUE), collapse = ""))
  
  # Initialize output
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
  
  
  
  
  # Check if filename exists and load last result (only if wait == TRUE)
  resultfile <- paste0(filename, "_result.RData")
  if(file.exists(resultfile) & wait) {
    load(file = resultfile, envir = .GlobalEnv, verbose = TRUE)
    return(out)
  }
  
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
  system(paste0("ssh ", machine, " mkdir ", filename, "_folder/"), ignore.stdout = TRUE, ignore.stderr = TRUE)
  system(paste0("ssh ", machine, " rm ", filename, "_folder/*"), ignore.stdout = TRUE, ignore.stderr = TRUE)
  system(paste0("scp ", getwd(), "/", filename, ".R* ", machine, ":", filename, "_folder/"))
  if(compile) {
    system(paste0("scp ", getwd(), "/*.c ", machine, ":", filename, "_folder/"))
  } else {
    system(paste0("scp ", getwd(), "/*.so ", machine, ":", filename, "_folder/"))
  }
  
  # Run in background
  system(paste0("ssh ", machine, " R CMD BATCH ", filename, "_folder/", filename, ".R --vanilla"), intern = FALSE, wait = wait)
  
  if(wait) {
    out$get()
    out$purge()
  } else {
    return(out)
  }
  
  
}

#' Generate the model objects for use in Xs (models with sensitivities)
#' 
#' @param f Named character vector with the ODE
#' @param forcings Character vector with the names of the forcings
#' @param fixed Character vector with the names of parameters (initial values and dynamic) for which
#' no sensitivities are required (will speed up the integration).
#' @param modelname Character, the name of the C file being generated.
#' @param verbose Print compiler output to R command line.
#' @param ... Further arguments being passed to funC.
#' @return list with \code{func} (ODE object) and \code{extended} (ODE+Sensitivities object)
#' @export
#' @import cOde
generateModel <- function(f, forcings=NULL, fixed=NULL, modelname = "f", verbose = FALSE, ...) {
  
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

#' Generate sample for multi-start fit
#' 
#' @param center named numeric, the center around we sample
#' @param samplefun character, indicating the random number generator,
#' defaults to \code{"rnorm"}.
#' @param fits length of the sample
#' @param ... arguments going to \code{samplefun}
#' @return matrix with the parameter samples
#' @export
mssample <- function(center, samplefun = "rnorm", fits = 20, ...) {
  
  sample.matrix <- do.call(rbind, lapply(1:fits, function(i) 
    do.call(samplefun, c(list(n = length(center), mean = center), list(...)))))
  
  colnames(sample.matrix) <- names(center)
  
  return(sample.matrix)

}

#' Open last plot in external pdf viewer
#' 
#' @description Convenience function to show last plot in an external viewer.
#' @param plot \code{ggplot2} plot object.
#' @param command character, indicatig which pdf viewer is started.
#' @param ... arguments going to \code{ggsave}.
#' @export
ggopen <- function(plot = last_plot(), command = "xdg-open", ...) {
  filename <- tempfile(pattern = "Rplot", fileext = ".pdf")
  ggsave(filename = filename, plot = plot, ...)
  system(command = paste(command, filename))
}
