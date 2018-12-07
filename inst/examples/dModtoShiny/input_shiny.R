# example script to produce a project folder for the shiny app - the app looks for all subfolders within the app-folder
# see CCD4 example folder

#' Save the necessary objects for dModtoShiny with the correct names
#'
#' @param x prediction function
#' @param parameters Parframe, result of mstrust()
#' @param fixed character vector with parameter names for fixed parameters 
#' @param data datalist
#' @param reactions eqnlist-object
#' @param profiles Parframe, result of profile() or a list of profiles
#' @param pubref Character, link to publication
#' @param errmodel obsfn
saveShiny <- function(reactions, 
                      x, 
                      parameters, 
                      fixed = NULL, 
                      data, 
                      profiles, 
                      pubref="none", 
                      errmodel = NULL){
  save(reactions, x, fixed, data, parameters, profiles, pubref, file = "input_shiny.RData")
}

#save all relevant data, etc needed for the app in a .RData file
saveShiny(x = (g*x*p), parameters = testtrust, data = mydata, reactions = f, profiles = myprofiles, pubref = "https://www.ncbi.nlm.nih.gov/pubmed/27811075")

#make a subfolder for this Project (dModtoShiny should contain app.R)
mymodel <- "Projectname"
system(paste0("ssh fermi mkdir /home/mfehling/ShinyApps/dModtoShiny/",mymodel))

# copy the .RData file to the subfolder and also all needed .so files
folder <- paste0("fermi:/home/mfehling/ShinyApps/dModtoShiny/",mymodel,"/")

system(paste0("scp input_shiny.RData ",folder,"."))

myso <- paste0(getLocalDLLs(),".so")
for(file in myso){
  system(paste0("scp ",myso," ",folder,"."))
}


# SaveShiny_dMod.frame ----
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
saveShiny_dMod.frame <- function(dMod.frame, hypothesis = 1, 
                                 reactions = dMod.frame$reactions[[hypothesis]], pubref = "none", fixed = dMod.frame$fixed[[hypothesis]],
                                 projectname
                                 ) {
  
  # TODO: Detect which names are present in dMod.frame. This way, stuff like "reactions" could be passed with the dMod.frame as well.
  
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
