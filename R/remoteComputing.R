
#' Detect number of free cores (on UNIX)
#' 
#' @description Read \code{/proc/loadavg} and subtract from the number of cores
#' @param machine character, e.g. "user@@localhost".
#' @export 
detectFreeCores <- function(machine = NULL) {
  
  if (!is.null(machine)) {
    nCores <- sapply(machine, function(m) {
      occupied <- as.numeric(strsplit(system(paste("ssh", m, "cat /proc/loadavg"), intern = TRUE), split = " ", fixed = TRUE)[[1]][1])  
      nCores <- as.numeric(system(paste("ssh", m, "nproc"), intern = TRUE))
      max(c(0, round(nCores - occupied)))
    })
  } else {
    occupied <- as.numeric(strsplit(system("cat /proc/loadavg", intern = TRUE), split = " ", fixed = TRUE)[[1]][1])      
    nCores <- as.numeric(system("nproc", intern = TRUE))
    nCores <- max(c(0, round(nCores - occupied)))
  }
  
  return(nCores)   
  
  
}

#' Run an R expression in the background (only on UNIX)
#' 
#' @description Generate an R code of the expression that is copied via \code{scp}
#' to any machine (ssh-key needed). Then collect the results.
#' @details \code{runbg()} generates a workspace from the \code{input} argument
#' and copies the workspace and all C files or .so files to the remote machines via
#' \code{scp}. This will only work if *an ssh-key had been generated and added
#' to the authorized keys on the remote machine*. On the remote machine, the script
#' will attempt to load all packages that had been loaded in the local R session.
#' This means that *all loaded packages must be present on the remote machine*. The
#' code snippet, i.e. the \code{...} argument, can include several intermediate results
#' but only the last call which is not redirected into a variable is returned via the
#' variable \code{.runbgOutput}, see example below.
#' @param ... Some R code
#' @param filename Character, defining the filename of the temporary file. Random
#' file name ist chosen if NULL.
#' @param machine Character vector, e.g. \code{"localhost"} or \code{"knecht1.fdm.uni-freiburg.de"}
#' or \code{c(localhost, localhost)}.
#' @param input Character vector, the objects in the workspace that are stored
#' into an R data file and copied to the remove machine.
#' @param compile Logical. If \code{TRUE}, C files are copied and compiled on the remote machine.
#' Otherwise, the .so files are copied.
#' @param wait Logical. Wait until executed. If \code{TRUE}, the code checks if the result file
#' is already present in which case it is loaded. If not present, \code{runbg()} starts, produces
#' the result and loads it as \code{.runbgOutput} directly into the workspace. If \code{wait = FALSE},
#' \code{runbg()} starts in the background and the result is only loaded into the workspace
#' when the \code{get()} function is called, see Value section. 
#' @return List of functions \code{check}, \code{get()} and \code{purge()}. 
#' \code{check()} checks, if the result is ready.
#' \code{get()} copies the result file
#' to the working directory and loads it into the workspace as an object called \code{.runbgOutput}. 
#' This object is a list named according to the machines that contains the results returned by each
#' machine.
#' \code{purge()} deletes the temporary folder
#' from the working directory and the remote machines.
#' @export
#' @examples
#' \dontrun{
#' out_job1 <- runbg({
#'          M <- matrix(rnorm(1e2), 10, 10)
#'          solve(M)
#'          }, machine = c("localhost", "localhost"), filename = "job1")
#' out_job1$check()          
#' out_job1$get()
#' result <- .runbgOutput
#' print(result)
#' out_job1$purge()
#' }
runbg <- function(..., machine = "localhost", filename = NULL, input = ls(.GlobalEnv), compile = FALSE, wait = FALSE) {
  
  
  expr <- as.expression(substitute(...))
  nmachines <- length(machine)
  
  # Set file name
  if (is.null(filename))
    filename <- paste0("tmp_", paste(sample(c(0:9, letters), 5, replace = TRUE), collapse = ""))
  
  filename0 <- filename
  filename <- paste(filename, 1:nmachines, sep = "_")
  
  # Initialize output
  out <- structure(vector("list", 3), names = c("check", "get", "purge"))
  out[[1]] <- function() {
    
    check.out <- sapply(1:nmachines, function(m) length(suppressWarnings(
      system(paste0("ssh ", machine[m], " ls ", filename[m], "_folder/ | grep -x ", filename[m], "_result.RData"), 
             intern = TRUE))))
    
    if (all(check.out) > 0) 
      cat("Result is ready!\n")
    else if (any(check.out) > 0)
      cat("Result from machines", paste(which(check.out > 0), collapse = ", "), "are ready.")
    else if (all(check.out) == 0)
      cat("Not ready!\n") 
      
  }
  
  out[[2]] <- function() {
    
    result <- structure(vector(mode = "list", length = nmachines), names = machine)
    for (m in 1:nmachines) {
      system(paste0("scp ", machine[m], ":", filename[m], "_folder/", filename[m], "_result.RData ./"), ignore.stdout = TRUE, ignore.stderr = TRUE)
      check <- try(load(file = paste0(filename[m], "_result.RData")), silent = TRUE) 
      if (!inherits("try-error", check)) result[[m]] <- .runbgOutput
    }
    
    .GlobalEnv$.runbgOutput <- result
    
  }
  
  out[[3]] <- function() {
    
    for (m in 1:nmachines) {
      system(paste0("ssh ", machine[m], " rm -r ", filename[m], "_folder"))
    }
    system(paste0("rm ", filename0, "*"))
  }
  
  
  # Check if filenames exist and load last result (only if wait == TRUE)
  resultfile <- paste(filename, "result.RData", sep = "_")
  if (all(file.exists(resultfile)) & wait) {
    
    for (m in 1:nmachines) {
      
      result <- structure(vector(mode = "list", length = nmachines), names = machine)
      load(file = resultfile[m])
      result[[m]] <- .runbgOutput
      
    }
    .GlobalEnv$.runbgOutput <- result
    return(out)
  }
  
  # Save current workspace
  save(list = input, file = paste0(filename0, ".RData"), envir = .GlobalEnv)
  
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
  if (compile)
    compile.line <- "cfiles <- list.files(pattern = '.c$'); for(cf in cfiles) system(paste('R CMD SHLIB', cf))"
   
  # Write program into character
  program <- lapply(1:nmachines, function(m) paste(
    pack,
    paste0("setwd('~/", filename[m], "_folder')"),
    "rm(list = ls())",
    compile.line,
    paste0("load('", filename0, ".RData')"),
    #".oldobjects <- ls()",
    paste0(".runbgOutput <- try(", as.character(expr), ")"),
    #".newobjects <- ls()",
    
    paste0("save(", output ,", file = '", filename[m], "_result.RData')"),
    sep = "\n"
  ))
  
  # Write program code into file
  for (m in 1:nmachines) cat(program[[m]], file = paste0(filename[m], ".R"))
  
  # Copy files to temporal folder
  for (m in 1:nmachines) {
    
    system(paste0("ssh ", machine[m], " mkdir ", filename[m], "_folder/"), ignore.stdout = TRUE, ignore.stderr = TRUE)
    system(paste0("ssh ", machine[m], " rm ", filename[m], "_folder/*"), ignore.stdout = TRUE, ignore.stderr = TRUE)
    system(paste0("scp ", getwd(), "/", filename0, ".RData* ", machine[m], ":", filename[m], "_folder/"))
    system(paste0("scp ", getwd(), "/", filename[m], ".R* ", machine[m], ":", filename[m], "_folder/"))
    if (compile) {
      system(paste0("scp ", getwd(), "/*.c ", machine[m], ":", filename[m], "_folder/"))
    } else {
      system(paste0("scp ", getwd(), "/*.so ", machine[m], ":", filename[m], "_folder/"))
    }
    
  }
  
  # Run in background
  for (m in 1:nmachines) system(paste0("ssh ", machine[m], " R CMD BATCH --vanilla ", filename[m], "_folder/", filename[m], ".R"), intern = FALSE, wait = wait)
  
  if (wait) {
    out$get()
    out$purge()
  } else {
    return(out)
  }
  
  
}

