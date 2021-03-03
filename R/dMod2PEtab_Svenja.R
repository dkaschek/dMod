# PEtab export stuff

#' Export observables as .tsv 
#' 
#' @description This function exports observables from dMod into the PEtab observable file.
#'  
#' @param observables dMod observable as named vector or eqnvec
#'   
#' @return PEtab observable file as .tsv
#'   
#' @author Svenja Kemmer
#'   
#' @export
#' 
writeObservablesTSV <- function(observables, errors, modelname){
  ## Create observables table
  obsIDs <- attr(observables, "names")
  obsFormula <- as.character(observables)
  obsTrafo <- attr(observables, "obsscales")
  errFormula <- as.character(errors)
  errDistribution <- "normal"
  
  obsFormula[which(obsTrafo !="lin")] <- gsub("^log\\(|^log10\\(|\\)$", "", obsFormula[which(obsTrafo !="lin")])
  
  obsDT <- data.table(observableId = obsIDs, 
                      observableFormula = obsFormula, 
                      observableTransformation = obsTrafo, 
                      noiseFormula = errFormula,
                      noiseDistribution = errDistribution
  )
  
  fwrite(obsDT, file = file.path(.exportFolder, paste0("observables_", modelname, ".tsv")), sep = "\t")
}


#' Export Measurement Data as .tsv 
#' 
#' @description This function imports data from the PEtab data file as a data list and defines errors if an error model is required.
#'  
#' @param data PEtab data file as .tsv
#' @param observables observables as eqnvec
#'   
#' @return data as data list and errors (if required) as eqnvec.
#'   
#' @author Marcus Rosenblatt and Svenja Kemmer
#'
#' @export
#' 
writeMeasurementsTSV <- function(data, errors){
  mydata <- mydata %>% as.data.frame() %>% as.data.table()
  DT <- mydata[,list(observableId = name, 
                     simulationConditionId = condition,
                     measurement = value,
                     time = time)]
  myobs <- read.csv(file = observables, sep = "\t") %>% as.data.frame()
  obs <- myobs$observableId %>% as.character()
  errors <- NULL
  
}

#' Get est.grid and fixed.grid
#'
#' @param jobname job name on cluster
#' @param recover 
#' @param nodetype 
#'
#' @return
#' @export
#'
#' @examples
getParGrids <- function(mytrafo = trafo, mytrafoL = trafoL, mycondition.grid = condition.grid, SS_pars = names(mySS_eqns_est)){
  
  # .. condition.grid -----
  
  myconditions <- mycondition.grid$condition
  rownames(mycondition.grid) <- myconditions
  
  # .. fixed.grid - all conditions fixed -----
  
  # select all trafo elements containing numbers
  fixed_trafo <- mytrafo[suppressWarnings(which(!(mytrafo %>% as.numeric()) %>% is.na()))] %>% as.eqnvec()
  fixed_df <- data.frame(as.list(fixed_trafo), stringsAsFactors=FALSE) 
  
  fixed_df2 <- fixed_df %>% rbind(fixed_df[rep(1, (length(myconditions)-1)), ]) %>% 
    cbind(condition = myconditions, ID = 1:length(myconditions)) %>% 
    select(ID, condition, everything()) 
  
  fixedpars <- names(fixed_trafo)
  
  # .. est.grid -----
  
  estpars <- names(mytrafo)[!names(mytrafo) %in% c(names(fixed_trafo), SS_pars)]
  est_trafo <- mytrafo[estpars]
  trafo_counterparts <- getSymbols(mytrafo)
  est_df <- NULL
  
  for(cond in myconditions){
    conditrafo <- mytrafoL[[cond]][estpars]
    if(any(c(str_detect(conditrafo, "exp\\("),str_detect(conditrafo, "10\\^\\(")))) conditrafo <- gsub("exp\\(", "", conditrafo) %>% gsub("10\\^\\(", "", .) %>% gsub("\\(", "", .) %>% gsub("\\)", "", .)
    if(any(c(str_detect(mytrafo, "exp\\("),str_detect(mytrafo, "10\\^\\(")))) mytrafo <- gsub("exp\\(", "", mytrafo) %>% gsub("10\\^\\(", "", .) %>%gsub("\\(", "", .) %>% gsub("\\)", "", .)
    
    # check for mathematical parameter trafos
    myoperations <- c("/|\\+|\\*")
    if(any(grepl(myoperations, conditrafo))){
      myreplpars <- grep(myoperations, conditrafo, value = TRUE)
      myorigpars <- grep(myoperations, mytrafo, value = TRUE)
      
      addpars <- NULL
      for(i in names(myreplpars)){
        myreplpar <- conditrafo[i]
        myorigpar <- mytrafo[i]
        parsorig <- strsplit(myorigpar, split = myoperations)[[1]]
        parsrepl <- strsplit(myreplpar, split = myoperations)[[1]]
        names(parsrepl) <- parsorig
        # check weather pars are already present in addpars
        for(j in names(parsrepl)){
          if(!(j %in% names(addpars))) addpars <- c(addpars, parsrepl[j])
        }  
        # remove numbers from original trafo
        addpars <- addpars[!grepl("^-?[[:digit:]]", names(addpars))]
      }
      
      excludenames <- intersect(c(names(addpars), names(myreplpars)), names(conditrafo))
      # substitute corresponding pars
      reduced_trafo <- subset(conditrafo, !names(conditrafo) %in% excludenames)
      conditrafo <- c(reduced_trafo, addpars)
    }
    
    # set all numbers/numerics = NA
    trafo.mod <- suppressWarnings(do.call(c, lapply(conditrafo, function(x) {fgh <- as.numeric(x); 
    if(is.na(fgh)){x} else NA})))
    
    mypars <- names(trafo.mod[!is.na(trafo.mod)])
    
    # filter condition specific est pars
    specifictrafo <- trafo.mod[mypars]
    # filter general est pars
    generaltrafo <- trafo.mod[mypars] %>% names()
    # filter general fixed pars
    numbertrafo <- setdiff(names(trafo.mod), mypars)
    # filter specific fixed numbers with general fixed pars as names
    numbers <- conditrafo[numbertrafo] %>% unclass()
    
    est_row <- data.frame(id = c(generaltrafo, names(numbers)), var = c(specifictrafo, numbers)) %>% data.table::transpose()
    names(est_row) <- est_row[1,]
    est_row <- est_row[-1,] %>% mutate(condition = cond)
    
    est_df <- rbind(est_df, est_row)
  }
  
  # .. fixed.grid - partly fixed conditions -----
  
  # filter fixedpars from est.grid 
  partly_fixed_df <- est_df %>% select_if(function(x) any(grepl("^-?[[:digit:]]",x))) 
  # assign NA to est pars
  partly_fixed_df <- sapply(partly_fixed_df, function(x) {
    sapply(x, function(z){
      if(!grepl("^-?[[:digit:]]", z)) z <- "NA"
      else z <- as.numeric(z)
    })
  })
  
  # append to fixed.grid
  fixed.grid <- fixed_df2 %>% cbind(partly_fixed_df %>% as.data.frame(stringsAsFactors = F)) %>% as_tibble()
  
  # assign "dummy" to fixed pars
  est_df <- sapply(est_df, function(x) {
    sapply(x, function(z){
      if(grepl("^-?[[:digit:]]", z)) z <- "dummy"
      else z <- z
    })
  })
  
  est.grid  <- est_df %>% as.data.frame(stringsAsFactors = F) %>% mutate(ID = 1:length(myconditions)) %>% 
    select(ID, condition, everything()) %>% as_tibble()
  
  list(est.grid, fixed.grid)
}


