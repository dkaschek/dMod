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
#' @example inst/examples/saveShiny_dMod.frame.R
saveShiny_dMod.frame <- function(dMod.frame, hypothesis = 1, 
                                 reactions = dMod.frame$reactions[[hypothesis]], pubref = "none", fixed = dMod.frame$fixed[[hypothesis]],
                                 projectname
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
  system(paste0("ssh fermi mkdir /home/mfehling/ShinyApps/dModtoShiny/", projectname))
  
  # copy the .RData file to the subfolder and also all needed .so files
  folder <- paste0("fermi:/home/mfehling/ShinyApps/dModtoShiny/", projectname, "/")
  
  system(paste0("scp input_shiny.RData ", folder, "."))
  
  models <- do.call(c, dMod.frame[hypothesis, drop = F]) %>%
    map(function(i) {
      mymodelname <- try(modelname(i), silent = T)
      if (!inherits(mymodelname, "try-error")) return(mymodelname)
      else return(NULL)
    }) %>% 
    do.call(c,.) %>%
    unique()
  .so <- .Platform$dynlib.ext
  files <- paste0(outer(models, c("", "_s", "_sdcv", "_deriv"), paste0), .so)
  files <- files[file.exists(files)]
  
  for(file in files){
    system(paste0("scp ", file, " ", folder, "."))
  }
  
  
}
