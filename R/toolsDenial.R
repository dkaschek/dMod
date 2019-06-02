#' SaveShiny for dMod.frames
#'
#' @param dMod.frame dmod.frame
#' @param hypothesis hypothesis
#' @param reactions eqnlist which are the basis for prd
#' @param pubref character. pubmed link or similar
#' @param fixed 
#' @param projectname Name of the folder being created on fermi
#' 
#' @export
#' 
#' 
#' @example inst/examples/saveShiny_dMod.frame.R
saveShiny_dMod.frame <- function(dMod.frame, hypothesis = 1, 
                                 reactions = dMod.frame$reactions[[hypothesis]], pubref = "none", fixed = dMod.frame$fixed[[hypothesis]],
                                 projectname, appfolder = "/home/mfehling/ShinyApps/dModtoShiny/"
) {
  
  # TODO: Detect which names are present in dMod.frame. This way, stuff like "reactions" could be passed with the dMod.frame as well.
  # TODO: Find a way to compare multiple hypotheses. 
      # Maybe allow for the option to pass dMod.frame$profile as complete list
      # But this doesn't solve the issue of comparing different steps -> Problem of dModtoShiny?
  
  
  saveShiny(x = dMod.frame$prd[[hypothesis]], 
            errmodel = dMod.frame$e[[hypothesis]],
            data = dMod.frame$data[[hypothesis]],
            parameters = dMod.frame$parframes[[hypothesis]], 
            profiles = dMod.frame$profiles[[hypothesis]], 
            
            reactions = reactions,
            pubref = pubref,
            fixed = fixed
  )
  
  
  # make a subfolder on the server (dModtoShiny should contain app.R)
  system(paste0("ssh fermi mkdir ", appfolder, projectname))
  
  # copy the .RData file to the subfolder and also all needed .so files
  folder <- paste0("fermi:", appfolder, projectname, "/")
  
  system(paste0("scp input_shiny.RData ", folder, "."))
  
  models <-  
    lapply(do.call(c, dMod.frame[hypothesis, drop = F]), function(i) {
      mymodelname <- try(modelname(i), silent = T)
      if (!inherits(mymodelname, "try-error")) return(mymodelname)
      else return(NULL)
    }) 
  models <- unique(do.call(c,models))
  
  .so <- .Platform$dynlib.ext
  files <- paste0(outer(models, c("", "_s", "_sdcv", "_deriv"), paste0), .so)
  files <- files[file.exists(files)]
  
  for(file in files){
    system(paste0("scp ", file, " ", folder, "."))
  }
  
  
}




#' Get the indices of the n largest (not necessarily best) steps of a parframe
#'
#' @param myparframe parframe, result from mstrust
#' @param nsteps number of steps
#' @param tol tolerance for stepdetection
#'
#' @return indices of the largest steps
#' @export
#' 
#' @seealso \link{parframe}
#' 
#' @importFrom stats setNames
#' 
#' @example inst/examples/getSteps.R
getStepIndices <- function(myparframe, nsteps = 5, tol = 1) {
  steps <- stepDetect(myparframe$value, tol)
  steps <- steps[order(c(diff(steps), nrow(myparframe)-max(steps)), decreasing = T)][1:nsteps]
  steps <- unique(sort(c(1, steps))) #include the first step no matter what
  setNames(steps, paste0("index", steps))
}


#' Get the rows of the n largest steps of a parframe
#'
#' @param myparframe parframe, result from mstrust
#' @param nsteps number of steps
#' @param tol tolerance for stepdetection
#'
#' @return parframe subsetted to the n largest steps
#' @export
#' 
#' @seealso \link{parframe}
#' 
#' @example inst/examples/getSteps.R
getSteps <- function(myparframe, nsteps = 5, tol = 1) {
  steps <- steps0 <- stepDetect(myparframe$value, tol)
  steps <- steps[order(c(diff(steps), nrow(myparframe)-max(steps)), decreasing = T)][1:nsteps]
  steps <- unique(sort(c(1, steps))) #include the first step no matter what
  
  if(length(steps0) <= nsteps) {
    nsteps <- length(steps0)
  }
  steps0 <- c(steps0, nrow(myparframe))
  steps <- c(steps, nrow(myparframe))
  steps_indices <- which(steps0%in%steps)
  
  step_members <- lapply(1:nsteps, function(i) {
    steps0[steps_indices[i]]:(steps0[steps_indices[i]+1]-1)
  }) 
  step_members <- do.call(c,step_members)
  
  myparframe[step_members,]
}



