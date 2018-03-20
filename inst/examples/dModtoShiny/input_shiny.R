# example script to produce a project folder for the shiny app - the app looks for all subfolders within the app-folder
# see CCD4 example folder

saveShiny <- function(x = y, parameters= myparameters, fixed = NULL, data = mydatalist, reactions = myeqnlist, profiles = myprofiles, pubref="none", errmodel = NULL){
  
  save(x, parameters,fixed, data, reactions, profiles, pubref, file = "input_shiny.RData")
  
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

