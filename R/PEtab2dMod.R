#' Import an SBML model and corresponding PEtab objects 
#' 
#' @description This function imports an SBML model and corresponding PEtab files, e.g. from the Benchmark collection.
#'  
#' @param modelname name of folder containing all PEtab files of the model to be imported. NULL if file paths are defined separately (see below).
#' @param path2model path to model folder
#' @param TestCases TRUE to load feature test cases
#' @param path2TestCases path to feature test case folder
#' @param compile if FALSE, g, ODEmodel and err are loaded from .RData (if present) and compilation time is saved
#' @param SBML_file SBML model as .xml
#' @param observable_file PEtab observable file as .tsv
#' @param condition_file PEtab condition file as .tsv
#' @param data_file PEtab data file as .tsv
#' @param parameter_file PEtab parameter file as .tsv
#'
#' @details Objects such as model equations, parameters or data are automatically assigned to the following standard variables and written to your current working directory (via <<-):
#' reactions, observables, errors, g, x, p0, err, obj, mydata, ODEmodel, condition.grid, trafoL, pouter, times.
#' Compiled objects (g, ODEmodel and err) are saved in .RData.
#' 
#' @return name of imported model
#'   
#' @author Marcus Rosenblatt and Svenja Kemmer
#' 
#' @export
#' 
importPEtabSBML <- function(modelname = "Boehm_JProteomeRes2014",
                            path2model = "BenchmarkModels/",
                            testCases = FALSE,
                            path2TestCases = "PEtabTests/",
                            compile = TRUE,
                            SBML_file = NULL,
                            observable_file = NULL,
                            condition_file = NULL,
                            data_file = NULL,
                            parameter_file = NULL
)
{
  
  ## Define path to SBML and PEtab files --------------------
  
  starttime <- Sys.time()
  if(testCases == FALSE){
    if(is.null(SBML_file))       SBML_file <- paste0(path2model, modelname, "/model_", modelname, ".xml")
    if(is.null(observable_file)) observable_file <- paste0(path2model, modelname, "/observables_", modelname, ".tsv")
    if(is.null(condition_file))  condition_file <- paste0(path2model, modelname, "/experimentalCondition_", modelname, ".tsv")
    if(is.null(data_file))       data_file <- paste0(path2model, modelname, "/measurementData_", modelname, ".tsv")
    if(is.null(parameter_file))  parameter_file <- paste0(path2model, modelname, "/parameters_", modelname, ".tsv")
  } else{
    SBML_file <- paste0(path2TestCases, modelname, "/_model.xml")
    observable_file <- paste0(path2TestCases, modelname, "/_observables.tsv")
    condition_file <- paste0(path2TestCases, modelname, "/_conditions.tsv")
    data_file <- paste0(path2TestCases, modelname, "/_measurements.tsv")
    parameter_file <- paste0(path2TestCases, modelname, "/_parameters.tsv")
  }
  mywd <- getwd()
  if(!file.exists(SBML_file)){cat(paste0("The file ",mywd,SBML_file, " does not exist. Please check spelling or provide the file name via the SBML_file argument.")); return(NULL)}
  if(!file.exists(observable_file)){cat(paste0("The file ",mywd,observable_file, " does not exist. Please check spelling or provide the file name via the observable_file argument.")); return(NULL)}
  if(!file.exists(condition_file)){cat(paste0("The file ",mywd,condition_file, " does not exist. Please check spelling or provide the file name via the condition_file argument.")); return(NULL)}
  if(!file.exists(data_file)){cat(paste0("The file ",mywd,data_file, " does not exist. Please check spelling or provide the file name via the data_file argument.")); return(NULL)}
  if(!file.exists(parameter_file)){cat(paste0("The file ",mywd,parameter_file, " does not exist. Please check spelling or provide the file name via the parameter_file argument.")); return(NULL)}
  
  if(is.null(modelname)) modelname <- "mymodel"
  ## Load shared objects --------------------
  
  dir.create(paste0(mywd,"/CompiledObjects/"), showWarnings = FALSE)
  setwd(paste0(mywd,"/CompiledObjects/"))
  files_loaded <- FALSE
  if(compile == FALSE & file.exists(paste0(modelname,".RData"))){
    load(paste0(modelname,".RData"))
    files_loaded <- TRUE
  } 
  setwd(mywd)
  
  ## Model Definition - Equations --------------------
  
  cat("Reading SBML file ...\n")
  mylist <- getReactionsSBML(SBML_file, condition_file)
  myreactions <- mylist$reactions
  myreactions_orig <- mylist$reactions_orig
  myevents <- mylist$events
  mypreeqEvents <- mylist$preeqEvents
  mystates <- mylist$mystates
  reactions <<- myreactions
  
  
  ## Model Definition - Observables --------------------
  
  cat("Reading observables ...\n")
  myobservables <- getObservablesSBML(observable_file)
  observables <<- myobservables
  
  
  cat("Compiling observable function ...\n")
  if(!files_loaded) {
    setwd(paste0(mywd,"/CompiledObjects/"))
    myg <- Y(myobservables, myreactions, compile=TRUE, modelname=paste0("g_",modelname))
    setwd(mywd)
  }
  g <<- myg
  
  
  ## Get Data ------------
  
  cat("Reading data file ...\n")
  mydataSBML <- getDataPEtabSBML(data_file, observable_file)
  mydata <- mydataSBML$data
  mydata <<- mydata
  
  
  ## Model Generation ---------------------
  
  cat("Compiling ODE model ...\n")
  
  if(!files_loaded) {
    setwd(paste0(mywd,"/CompiledObjects/"))
    myodemodel <- odemodel(myreactions, forcings = NULL, events = myevents, fixed=NULL, modelname = paste0("odemodel_", modelname), jacobian = "inz.lsodes", compile = TRUE)
    setwd(mywd)
  }
  ODEmodel <<- myodemodel
  
  
  ## Check and define error model ------------ 
  
  cat("Check and compile error model ...\n")
  myerrors <- mydataSBML$errors
  errors <<- myerrors
  
  myerr <- NULL
  if(!files_loaded) {
    if(!is_empty(getSymbols(myerrors))){
      setwd(paste0(mywd,"/CompiledObjects/"))
      myerr <- Y(myerrors, f = c(as.eqnvec(myreactions), myobservables), states = names(myobservables), attach.input = FALSE, compile = TRUE, modelname = paste0("errfn_", modelname))
      setwd(mywd)
    }
  }
  err <<- myerr
  
  ## Define constraints, initials, parameters and compartments --------------
  
  cat("Reading parameters and initials ...\n")
  myparameters <- getParametersSBML(parameter_file, SBML_file)
  myconstraints <- myparameters$constraints
  SBMLfixedpars <- myparameters$SBMLfixedpars
  myfit_values <- myparameters$pouter
  myinitialsSBML <- getInitialsSBML(SBML_file, condition_file)
  mycompartments <- myinitialsSBML$compartments
  myinitials <- myinitialsSBML$initials
  
  ## Parameter transformations -----------
  
  # Generate condition.grid
  grid <- getConditionsSBML(conditions = condition_file, data = data_file) 
  mypreeqCons <- grid$preeqCons
  mycondition.grid <- grid$condition_grid
  
  if(!is.null(SBMLfixedpars)){
    for (i in 1:length(SBMLfixedpars)) {
      if(!names(SBMLfixedpars)[i] %in% names(mycondition.grid))  mycondition.grid[names(SBMLfixedpars)[i]] <- SBMLfixedpars[i]
    } 
  }
  condi_pars <- names(mycondition.grid)[!names(mycondition.grid) %in% c("conditionName","conditionId")]
  condition.grid <<- mycondition.grid
  
  cat("Generate parameter transformations ...\n")
  myinnerpars <- unique(c(getParameters(myodemodel), getParameters(myg), getSymbols(myerrors)))
  names(myinnerpars) <- myinnerpars
  trafo <- as.eqnvec(myinnerpars, names = myinnerpars)
  trafo <- replaceSymbols(names(mycompartments), mycompartments, trafo)
  # only overwrite intial if it's not defined in condition.grid
  for (i in 1:length(myinitials)) {
    if(!names(myinitials)[i] %in% names(mycondition.grid)){
      trafo <- replaceSymbols(names(myinitials)[i], myinitials[i], trafo)
    }
  }
  trafo <- replaceSymbols(names(myconstraints), myconstraints, trafo)
  
  # branch trafo for different conditions
  mytrafoL <- branch(trafo, table=mycondition.grid)
  # set preequilibration event initials to corresponding values
  if(!is.null(mypreeqEvents)){
    mypreeqEvents2replace <- filter(mypreeqEvents, !var%in%mystates)
    mytrafoL <- repar("x~y", mytrafoL , x = unique(mypreeqEvents2replace$var), y = attr(mypreeqEvents2replace, "initials"))
  }
  # set remaining event initial to 0
  mytrafoL <- repar("x~0", mytrafoL , x = setdiff(unique(myevents$var), unique(mypreeqEvents$var)))  
  
  # condition-specific assignment of parameters from condition grid
  if(length(condi_pars) > 0){
    for (j in 1:length(names(mytrafoL))) {
      for (i in 1:length(condi_pars)) {
        mytrafoL[[j]] <- repar(x~y, mytrafoL[[j]], x=condi_pars[i], y=mycondition.grid[j,condi_pars[i]])
      }
    }
  }
  
  
  # transform parameters according to scale defined in the parameter PEtab file
  parscales <- attr(myfit_values,"parscale")
  mynames <- names(parscales)
  for(i in 1:length(parscales)){
    par <- parscales[i]
    par[par=="lin"] <- ""
    par[par=="log10"] <- "10**"
    par[par=="log"] <- "exp"
    parameter <- mynames[i]
    mytrafoL <- repar(paste0("x~",par,"(x)"), mytrafoL, x = parameter)
  }
  trafoL <<- mytrafoL
  
  
  ## Specify prediction functions ------
  
  # Get numeric steady state for preequilibration conditions
  if(!is.null(mypreeqCons)){
    myf <- as.eqnvec(myreactions_orig)[mystates]
    cq <- conservedQuantities(myreactions_orig$smatrix)
    if(!is.null(cq)){
      for(i in 1:nrow(cq)){
        myf[getSymbols(cq)[1]] <- paste0(as.character(conservedQuantities(myreactions_orig$smatrix)[1,]),"-1")
      }
    }
    setwd(paste0(mywd,"/CompiledObjects/"))
    pSS <- P(myf, condition = "c0", method = "implicit", compile = TRUE, modelname = paste0("preeq_", modelname))
    setwd(mywd)
  } else pSS <- NULL
  
  cat("Generate prediction function ...\n")
  tolerances <- 1e-7
  myp0 <- myx <- NULL
  for (C in names(mytrafoL)) {
    if(C%in%mypreeqCons){
      myp0 <- myp0 + pSS*P(mytrafoL[[C]], condition = C)
    } else {
      myp0 <- myp0 + P(mytrafoL[[C]], condition = C)
    }
    myx <- myx + Xs(myodemodel, optionsOde = list(method = "lsoda", rtol = tolerances, atol = tolerances, maxsteps = 5000),
                    optionsSens = list(method = "lsodes", lrw=200000, rtol = tolerances, atol = tolerances),
                    condition = C)
  }
  
  p0 <<- myp0
  x <<- myx
  
  ## Generate objective function and initial parameter set -------
  
  myouterpars <- getParameters(myp0)
  mypouter <- structure(rep(NA,length(myouterpars)), names = myouterpars)
  
  common <- intersect(names(mypouter),names(myfit_values))
  mypouter[common] <- myfit_values[common]
  attr(mypouter, "parscales") <- parscales[which(names(myfit_values)%in%names(mypouter))]
  attr(mypouter, "lowerBound") <- attr(myfit_values,"lowerBound")[which(names(myfit_values)%in%names(mypouter))]
  attr(mypouter, "upperBound") <- attr(myfit_values,"upperBound")[which(names(myfit_values)%in%names(mypouter))]
  pouter <<- mypouter
  
  
  ## Define objective function -------
  
  cat("Generate objective function ...\n")
  if(!is.null(myerrors)){
    myobj <- normL2(mydata, myg*myx*myp0, errmodel = myerr) #+ constraintL2(prior, sigma=16)
  } else myobj <- normL2(mydata, myg*myx*myp0)
  obj <<- myobj
  
  mytimes <- seq(0,max(do.call(c, lapply(1:length(mydata), function(i) max(mydata[[i]]$time)))), len=501)
  times <<- mytimes
  
  if(!files_loaded){
    setwd(paste0(mywd,"/CompiledObjects/"))
    save(list = c("myg","myodemodel","myerr"),file = paste0(modelname,".RData"))
    setwd(mywd)
  }
  
  model_name <<- modelname
  
  endtime <- Sys.time()
  mytimediff <- as.numeric(difftime(endtime, starttime, unit="secs"))
  if(mytimediff > 3600) cat(paste0(modelname, " imported in ",as.character(format(as.numeric(difftime(endtime, starttime, unit="hours")), digits=3)), " hours.\n")) else
    if(mytimediff > 60) cat(paste0(modelname, " imported in ",as.character(format(as.numeric(difftime(endtime, starttime, unit="mins")), digits=3)), " minutes.\n")) else
      cat(paste0(modelname, " imported in ",as.character(format(as.numeric(difftime(endtime, starttime, unit="secs")), digits=3)), " seconds.\n"))
  return(modelname)
}



#' Fit a model imported via importPEtabSBML 
#' 
#' @description A wrapper function to use \link{mstrust} with imported PEtabSBML models. Some reasonable standard arguments for mstrust are used. Results of mstrust are written to Results folder.
#'  
#' @param objfun Objective function to be minimized as created by \link{importPEtabSBML}.
#' @param nrfits numeric, Number of fits to be performed
#' @param nrcores numeric, Number of cores to be used
#' @param useBounds  boolean, if TRUE, parameter bounds are taken as provided in PEtab format, if FALSE no parameter bounds are applied
#'   
#' @return parframe with the parameter estimated of the multi-start optimization
#'   
#' @author Marcus Rosenblatt and Svenja Kemmer
#' 
#' @export
#'   
fitModelPEtabSBML <- function(objfun=obj, nrfits=4, nrcores=4, useBounds=TRUE){
  prior <- structure(rep(0,length(pouter)))
  names(prior) <- names(pouter)
  mywd <- getwd()
  dir.create(paste0(mywd,"/Test/mstrust/"), showWarnings = FALSE)
  if(useBounds) out <- mstrust(objfun=objfun, center=msParframe(prior, n = nrfits+1, seed=47)[-1,], studyname=model_name, rinit = 0.1, rmax = 10,
                               fits = nrfits, cores = nrcores, samplefun = "rnorm", resultPath = "Test/mstrust/",
                               parlower = attr(pouter, "lowerBound"), parupper=attr(pouter, "upperBound"),
                               stats = FALSE, narrowing = NULL, iterlim=400, sd = 3)
  else out <- mstrust(objfun=objfun, center=msParframe(prior, n = nrfits, seed=47), studyname=model_name, rinit = 0.1, rmax = 10,
                      fits = nrfits, cores = nrcores, samplefun = "rnorm", resultPath = "Test/mstrust/",
                      stats = FALSE, narrowing = NULL, iterlim=400, sd = 3)
  if(any(lapply(out, function(fgh) fgh$converged)==TRUE)) return(as.parframe(out)) else {cat("No fit converged."); return(NULL)}
}



#' Test PEtabSBML import
#'
#' @description This function imports, evaluates and tests the PEtab model.
#'
#' @param models model names to test
#'
#' @return evaluation data.frame
#'
#' @author Marcus Rosenblatt and Svenja Kemmer
#'
#' @export
#'
testPEtabSBML <- function(models = c(
  #"Boehm_JProteomeRes2014"
  # "Fujita_SciSignal2010",
  # "Borghans_BiophysChem1997",
  # "Elowitz_Nature2000",
  # "Sneyd_PNAS2002",
  # "Crauste_CellSystems2017",
  # "Schwen_PONE2014",
  # "Raia_CancerResearch2011",
  # "Zheng_PNAS2012",
  # "Beer_MolBioSystems2014",
  # "Brannmark_JBC2010",
  # "Bruno_JExpBio2016",
  # "Chen_MSB2009",
  # "Fiedler_BMC2016",
  # "Weber_BMC2015",
  # "Swameye_PNAS2003"
  # "Bachmann_MSB2011"
  # "Lucarelli_CellSystems2018",
  "0001",
  "0002",
  "0003",
  "0004",
  "0005",
  "0006",
  "0007",
  "0008",
  "0009",
  "0010",
  "0011",
  "0012",
  "0013",
  "0014",
  "0015",
  "0016"
), testFit = TRUE, timelimit = 5000, testCases = FALSE) {
  try_with_time_limit <- function(expr, cpu = Inf, elapsed = Inf) {
    y <- try(
      {
        setTimeLimit(cpu, elapsed)
        expr
      },
      silent = TRUE
    )
    if (inherits(y, "try-error")) NULL else y
  }
  cat(green("Start test function...\n"))
  mywd <- getwd()
  teststarttime <- Sys.time()
  output <- NULL
  predictions <- NULL
  for (model in models) {
    setwd(mywd)
    importtest <- F
    plottest <- F
    bestfit <- NA
    cat(blue(paste0("Testing ", model, "\n")))
    fgh <- try_with_time_limit(
      {
        test <- try(importPEtabSBML(model, compile = T, testCases = testCases), silent = T)
        if (inherits(test, "try-error")) "import error" else test
      },
      timelimit
    )
    if (fgh == "import error") {
      cat(yellow("Import error or time limit exceeded for", model, "\n\n\n"))
      output <- rbind(output, data.frame(
        modelname = model, import = importtest,
        fitting_time = NA, plot = plottest, chi2 = NA, LL = NA, bestfit = NA, difference = NA
      ))
    } else {
      importtest <- T
      testobj <- try(obj(pouter))
      if (inherits(testobj, "try-error")) {
        cat(red("Warning: Error in calculation of objective function.\n"))
        output <- rbind(output, data.frame(
          modelname = model, import = importtest,
          fitting_time = NA, plot = plottest, chi2 = NA, LL = NA, bestfit = NA, difference = NA
        ))
      } else {
        if (is.numeric(testobj$value)) {
          cat(green("Calculation of objective function successful.\n"))
          if (testCases){
            # calculate predictions for trajectory comparison
            mysimulations <- read.csv(paste0("PEtabTests/", model, "/_simulations.tsv"), sep = "\t")
            simu_time <- unique(mysimulations$time)
            prediction <- (g*x*p0)(simu_time, pouter)
            predictions <- rbind(predictions, data.frame(
              modelname = model, pred = prediction, obs.transformation = NA
            ))
            # append observable scale to predictions
            for (i in 1:length(observables)) {
              scale <- attr(observables, "obsscales")[i]
              predictions <- predictions %>% mutate(obs.transformation = ifelse(modelname == model & pred.name == names(observables)[i], scale, obs.transformation))
            }
          }
        } else {
          cat(red("Warning: obj(pouter) is not numeric.\n"))
        }
        # objLL <- mynormL2(mydata, g * x * p0, outputLL = T)
        # testLL <- try(-0.5 * objLL(pouter)$value)
        # if (inherits(testLL, "try-error")) testLL <- NA
        if (testFit) {
          fitstarttime <- Sys.time()
          myframe <- fitModelPEtabSBML(nrfits = 20)
          fitendtime <- Sys.time()
          if (is.parframe(myframe) & !is.null(myframe)) {
            if (is.numeric(obj(myframe[1, ])$value)) {
              cat(green("Fit test successful.\n"))
              bestfit <- obj(myframe[1, ])$value
            } else {
              cat(red("Warning: obj(myframe) is not numeric.\n"))
            }
          } else {
            cat(red("Warning: Fit test not successful..\n"))
          }
          mytimediff <- as.numeric(difftime(fitendtime, fitstarttime, unit = "secs"))
          if (mytimediff > 3600) {
            cat(green(paste0("Fitting done in ", as.character(format(as.numeric(difftime(fitendtime, fitstarttime, unit = "hours")), digits = 3)), " hours.\n")))
          } else
            if (mytimediff > 60) {
              cat(green(paste0("Fitting done in ", as.character(format(as.numeric(difftime(fitendtime, fitstarttime, unit = "mins")), digits = 3)), " minutes.\n")))
            } else {
              cat(green(paste0("Fitting done in ", as.character(format(as.numeric(difftime(fitendtime, fitstarttime, unit = "secs")), digits = 3)), " seconds.\n")))
            }
        }
        pdf(file = paste0("Test/", model, "_plotAll.pdf"))
        plotPEtabSBML()
        dev.off()
        pdf(file = paste0("Test/", model, "_plotTargetsObserved.pdf"))
        plotPEtabSBML(name %in% names(observables))
        dev.off()
        pdf(file = paste0("Test/", model, "_plotConditionsObserved.pdf"))
        plotPEtabSBML(condition %in% names(mydata))
        dev.off()
        plottest <- T
        cat(green("Import and plot test for ", fgh, " successful!\n\n\n"))
        
        output <- rbind(output, data.frame(
          modelname = model, import = importtest,
          # fitting_time = format(as.numeric(difftime(fitendtime, fitstarttime, unit = "mins")), digits = 3),
          plot = plottest, chi2 = attr(testobj,"chisquare"), LL = -0.5*testobj$value
          # , bestfit = bestfit, difference = bestfit - testobj$value
        ))
      }
    }
    if (testCases){
      sharedObjects <- paste0("CompiledObjects/",
                              c(paste0("g_",model,".so"),
                                paste0("g_",model,"_deriv.so"),
                                paste0("odemodel_",model,".so"),
                                paste0("odemodel_",model,"_s.so"),
                                paste0("errfn_",model,".so"),
                                paste0("errfn_",model,"_s.so")))
      for (file in sharedObjects) if(file.exists(file)) try(dyn.unload(file))
    }
  }
  if (testCases) {
    simu_output <- output[1]
    output <- cbind(output, chi2_sol = NA, tol_chi2_sol = NA, LL_sol = NA, tol_LL_sol = NA)
    for (model in models) {
      mysolution <- read_yaml(paste0("PEtabTests/", model, "/_", model, "_solution.yaml"))
      output[which(output$modelname == model), "chi2_sol"] <- mysolution$chi2
      output[which(output$modelname == model), "tol_chi2_sol"] <- mysolution$tol_chi2
      output[which(output$modelname == model), "LL_sol"] <- mysolution$llh
      output[which(output$modelname == model), "tol_LL_sol"] <- mysolution$tol_llh
      
      # extract simulation values
      simu_output[which(simu_output$modelname == model), "tol_simus_sol"] <- mysolution$tol_simulations
      mysimulations <- read.csv(paste0("PEtabTests/", model, "/_simulations.tsv"), sep = "\t")
      simu_prediction <- subset(predictions, modelname == model)
      
      # iterate through simulation points
      for (nrow in 1:nrow(mysimulations)) {
        simu_row <- mysimulations[nrow,]
        simu_time <- simu_row$time
        simu_obs <- simu_row$observableId %>% as.character()
        simu_condi <- simu_row$simulationConditionId %>% as.character()
        simu_obspars <- simu_row$observableParameters
        
        if(!is.null(simu_obspars) & length(unique(simu_prediction$pred.condition)) > 1){
          simu_condi <- paste0(simu_condi, "_", simu_obspars)
        }
        pred_row <- subset(simu_prediction, pred.time == simu_time & pred.name == simu_obs & pred.condition == simu_condi)
        # retransform simulation value according to observable transformation
        if(pred_row$obs.transformation == "log10") pred_row$pred.value <- 10**pred_row$pred.value
        if(pred_row$obs.transformation == "log") pred_row$pred.value <- exp(pred_row$pred.value)
        simu_value <- pred_row$pred.value
        
        simu_output[which(simu_output$modelname == model), paste0("simu_", nrow)] <- simu_value
        simu_output[which(simu_output$modelname == model), paste0("simu_", nrow,"_sol")] <- mysimulations$simulation[nrow]
      }
    }
  }
  
  testendtime <- Sys.time()
  mytimediff <- as.numeric(difftime(testendtime, teststarttime, unit = "secs"))
  if (mytimediff > 3600) {
    cat(green(paste0("Test done in ", as.character(format(as.numeric(difftime(testendtime, teststarttime, unit = "hours")), digits = 3)), " hours.\n")))
  } else
    if (mytimediff > 60) {
      cat(green(paste0("Test done in ", as.character(format(as.numeric(difftime(testendtime, teststarttime, unit = "mins")), digits = 3)), " minutes.\n")))
    } else {
      cat(green(paste0("Test done in ", as.character(format(as.numeric(difftime(testendtime, teststarttime, unit = "secs")), digits = 3)), " seconds.\n")))
    }
  
  if (testCases) {
    
    # check simulations
    for (model in models) {
      correctORnot <- NULL
      modelrow <- subset(simu_output, modelname == model)
      modelrow_woNA <- modelrow[colSums(!is.na(modelrow)) > 0]
      for (ncol in seq(3,(ncol(modelrow_woNA)),2)) {
        simuCompare <- abs(modelrow_woNA[[ncol]]-modelrow_woNA[[ncol+1]]) < modelrow_woNA$tol_simus_sol
        correctORnot <- c(correctORnot, simuCompare)
      }
      if(length(unique(correctORnot)) == 1){
        SimuPassed <- unique(correctORnot)
      } else SimuPassed <- FALSE
      simu_output[which(simu_output$modelname == model), "Passed"] <- SimuPassed
    }
    
    output <- cbind(output,
                    X2Passed = (abs(output$chi2 - output$chi2_sol) < output$tol_chi2_sol),
                    LLPassed = (abs(output$LL - output$LL_sol) < output$tol_LL_sol)
    )
  }
  
  if (!testCases) simu_output <- NULL
  return(list(output = output,simu_output = simu_output))
}



#' Import condition.grid from PEtab 
#' 
#' @description This function imports the experimental conditions from the PEtab condition file as a gondition.grid.
#'  
#' @param conditions PEtab condition file as .tsv
#'   
#' @return condition.grid as data frame.
#'   
#' @author Marcus Rosenblatt and Svenja Kemmer
#'   
#' @export
#' 
getConditionsSBML <- function(conditions,data){
  condition.grid_orig <- read.csv(file = conditions, sep = "\t")
  mydata <- read.csv(file = data, sep = "\t")
  
  # handle preequilibration conditions
  myCons <- condition.grid_orig$conditionId
  mypreeqCons <- NULL
  for (con in myCons){if(paste0("preeq_", con)%in%myCons) mypreeqCons <- c(mypreeqCons, con)}
  condition.grid_orig <- filter(condition.grid_orig, !conditionId%in%mypreeqCons)
  condition.grid_orig$conditionId <- sub("preeq_", "", condition.grid_orig$conditionId)
  
  # check which conditions are observed
  condis_obs <- mydata$simulationConditionId %>% unique()
  # check which observables exist
  observables <- mydata$observableId %>% unique()
  
  # replace "" by NA
  if(!is.null(mydata$observableParameters)){
    mydata$observableParameters <- mydata$observableParameters %>% as.character()
    mydata <- mydata %>% mutate(observableParameters = ifelse(observableParameters == "",NA,observableParameters))
  }
  if(!is.null(mydata$noiseParameters)){
    mydata$noiseParameters <- mydata$noiseParameters %>% as.character()
    mydata <- mydata %>% mutate(noiseParameters = ifelse(noiseParameters == "",NA,noiseParameters))
  }
  # generate columns for observableParameters
  if(!is.numeric(mydata$observableParameters) & !is.null(mydata$observableParameters)){
    condition.grid_obs <- data.frame(conditionId = condis_obs)
    for (obs in observables){
      data_obs <- subset(mydata, observableId == obs)
      for (condition in condis_obs){
        if(condition %in% data_obs$simulationConditionId){
          row_pars <- NULL
          obs_par <- subset(data_obs, simulationConditionId == condition)$observableParameters %>% unique() %>% as.character()
          if(length(obs_par)==1){
            if(!is.na(obs_par)){
              # one or more observable parameters?
              if(str_detect(obs_par,";")){
                myobspars <- strsplit(obs_par,";")[[1]]
                for(i in 1:length(myobspars)) {
                  row_pars <- c(row_pars, myobspars[i])
                }
              } else row_pars <- c(row_pars, obs_par)
            }
            if(!is.null(row_pars)) for (par in 1:length(row_pars)) {
              col_name <- paste0("observableParameter",par,"_",obs)
              condition.grid_obs[which(condition.grid_obs$conditionId==condition),col_name] <- row_pars[par]
            }
          } else {
            col_name <- paste0("observableParameter1_",obs)
            add <- NULL
            for(j in 2:length(obs_par)){
              add <- rbind(add, subset(condition.grid_obs, conditionId==condition))
            }
            condition.grid_obs <- rbind(condition.grid_obs, add)
            condition.grid_obs[which(condition.grid_obs$conditionId==condition),col_name] <- obs_par
            condition.grid_obs$conditionId <- as.character(condition.grid_obs$conditionId)
            condition.grid_obs$conditionId[which(condition.grid_obs$conditionId==condition)] <- paste0(condition,"_", obs_par)
            
            condition.grid_orig <- rbind(condition.grid_orig, add)
            condition.grid_orig$conditionId <- as.character(condition.grid_orig$conditionId)
            condition.grid_orig$conditionId[which(condition.grid_orig$conditionId==condition)] <- paste0(condition,"_", obs_par)
          }
        }
      } 
    }
    mycondition.grid <- suppressWarnings(inner_join(condition.grid_orig,condition.grid_obs, by = "conditionId"))
    # avoid warning if not all conditions are observed
  } else mycondition.grid <- condition.grid_orig
  
  # generate columns for noiseParameters
  if(!is.numeric(mydata$noiseParameters) & !is.null(mydata$noiseParameters)) 
  {
    if(exists("mycondition.grid")) {condition.grid_orig <- mycondition.grid}
    condition.grid_noise <- data.frame(conditionId = condis_obs)
    for (obs in observables) 
    {
      data_obs <- subset(mydata, observableId == obs)
      for (condition in condis_obs) 
      {
        if(condition %in% data_obs$simulationConditionId){
          row_pars <- NULL
          noise_par <- subset(data_obs, simulationConditionId == condition)$noiseParameters %>% unique() %>% as.character()
          if(!is.na(noise_par)){
            # one or more observable parameters?
            if(str_detect(noise_par,";")){
              myobspars <- strsplit(noise_par,";")[[1]]
              for(i in 1:length(myobspars)) {
                row_pars <- c(row_pars, myobspars[i])
              }
            } else row_pars <- c(row_pars, noise_par)
          }
          if(!is.null(row_pars)) for (par in 1:length(row_pars)) {
            col_name <- paste0("noiseParameter",par,"_",obs)
            condition.grid_noise[which(condition.grid_noise$conditionId==condition),col_name] <- row_pars[par]
          }
        }
      } 
      
    }
    mycondition.grid <- suppressWarnings(inner_join(condition.grid_orig,condition.grid_noise, by = "conditionId"))
    # avoid warning if not all conditions are observed
  }
  
  if(!exists("mycondition.grid")) mycondition.grid <- condition.grid_orig
  rownames(mycondition.grid) <- mycondition.grid$conditionId
  # mycondition.grid$conditionId <- NULL ## we need this column in cases with just one condition!
  
  # check if all conditions are observed
  if(nrow(mycondition.grid) < nrow(condition.grid_orig)) print("There exist non-observed conditions!")
  
  for(i in 1:nrow(mycondition.grid)){
    for(j in 1:ncol(mycondition.grid)){
      if(is.na(mycondition.grid[i,j])) mycondition.grid[i,j] <- "1"
    }
  }
  
  return(list(condition_grid=mycondition.grid, preeqCons=mypreeqCons))
}



#' Import initials from SBML. 
#' 
#' @description This function imports initial values or equations describing the same from SBML and writes them in a named vector. 
#'  
#' @param model SBML file as .xml
#'   
#' @return Named vector of initials.
#'   
#' @author Marcus Rosenblatt, Svenja Kemmer and Frank Bergmann
#'   
#' @export
#' 
getInitialsSBML <- function(model, conditions){
  
  condition.grid <- read.csv(file = conditions, sep = "\t")
  model = readSBML(model)$getModel()
  initials <- NULL
  species <- NULL
  
  # now we can go through all species
  for ( i in 0:(model$getNumSpecies()-1) ){
    current <- model$getSpecies(i)
    
    # now the species can have several cases that determine
    # their initial value 
    
    # it could be that the species is fully determined by an assignment rule
    # (that apply at all times), so we have to check rules first
    rule <- model$getRule(current$getId())
    if (!is.null(rule))
    {
      # ok there is a rule for this species so lets figure out what its type
      # is as that determines whether it applies at t0
      rule_type <- rule$getTypeCode()
      type_name <- libSBML::SBMLTypeCode_toString(rule_type, 'core')
      if (type_name == "AssignmentRule")
      {
        # the initial value is determined by the formula
        math <- rule$getMath()
        if (!is.null(math))
        {
          formula <- libSBML::formulaToL3String(math)
          # print(paste('Species: ', current$getId(), ' is determined at all times by formula: ', formula))
          initials <- c(initials,formula)
          species <- c(species,current$getId())
          
          # no need to look at other values so continue for another one
          next
        }
      }
      
      if (type_name == "RateRule")
      {
        math <- rule$getMath()
        if (!is.null(math))
        {
          formula <- libSBML::formulaToL3String(math)
          # print(paste('Species: ', current$getId(), ' has an ode rule with formula: ', formula))
          initials <- c(initials,formula)
          species <- c(species,current$getId())
          
          # even though there is an ODE attached to the species, its initial value is needed
        }
      }
    }
    
    
    # it could have an initial assignment
    ia <- model$getInitialAssignment(current$getId())
    if (!is.null(ia))
    {
      math <- ia$getMath()
      if (!is.null(math))
      {
        formula <- libSBML::formulaToL3String(math)
        # print(paste("Species: ", current$getId(), " has an initial assignment with formula: ", formula))
        initials <- c(initials,formula)
        species <- c(species,current$getId())
        
        # as soon as you have that formula, no initial concentration / amount applies
        # so we don't have to look at anything else for this species
        next
      }
    }
    
    
    # it could have an initial amount
    if (current$isSetInitialAmount())
    {
      # print (paste("Species: ", current$getId(), "has initial amount: ", current$getInitialAmount()))
      initials <- c(initials,current$getInitialAmount())
      species <- c(species,current$getId())
    }
    
    # it could have an initial concentration
    if (current$isSetInitialConcentration())
    {
      # print (paste("Species: ", current$getId(), "has initial concentration: ", current$getInitialConcentration()))
      initials <- c(initials,current$getInitialConcentration())
      species <- c(species,current$getId())
    }
  }
  
  names(initials) <- species
  
  
  ## extract compartments
  
  # check if compartments exist
  if(model$getNumCompartments()>0)
  {
    
    #initialize vectors
    comp_name <- NULL
    comp_size <- NULL
    
    for ( i in 0:(model$getNumCompartments()-1) )
    {
      # get compartment name and size
      which <- model$getCompartment(i)$getId()
      if(which %in% names(condition.grid)){
        size <- condition.grid[which] %>% as.numeric()
      } else size <- model$getCompartment(i)$getSize()
      
      comp_name <- c(comp_name,which)
      comp_size <- c(comp_size,size)
      
    }
    compartments <- comp_size
    names(compartments) <- comp_name
  }
  
  return(list(initials = initials, compartments = compartments))
}


#' Import observables from PEtab. 
#' 
#' @description This function imports observables from the PEtab observable file.
#'  
#' @param observables PEtab observable file as .tsv
#'   
#' @return Eqnvec of observables.
#'   
#' @author Marcus Rosenblatt and Svenja Kemmer
#'   
#' @export
#' 
getObservablesSBML <- function(observables){
  ## Load observables
  myobs <- read.csv(file = observables, sep = "\t") %>% as.data.frame()
  obsNames <- myobs$observableId %>% as.character()
  
  # # rename observables with _obs
  # obsNames <- paste0(obsNames,"_obs")
  
  obsFormula <- myobs$observableFormula %>% as.character()
  obsFormula[which(myobs$observableTransformation=="log")] <- paste0("log(", obsFormula[which(myobs$observableTransformation=="log")], ")")
  obsFormula[which(myobs$observableTransformation=="log10")] <- paste0("log10(", obsFormula[which(myobs$observableTransformation=="log10")], ")")
  names(obsFormula) <- obsNames
  observables <- obsFormula %>% as.eqnvec()
  
  # collect observable transformations as attribute
  if(!is.null(myobs$observableTransformation)){
    obsscales <- myobs$observableTransformation %>% as.character()
  } else obsscales <- rep("lin", length(obsNames))
  
  attr(observables,"obsscales") <- obsscales
  return(observables)
}


#' Import Parameters from PEtab 
#' 
#' @description This function imports fixed and fitted parameters from the PEtab parameter file as named vectors.
#'  
#' @param parameters PEtab parameter file as .tsv
#'   
#' @return constraints and pouter as list of named vectros.
#'   
#' @author Marcus Rosenblatt and Svenja Kemmer
#'   
#' @export
#' 
getParametersSBML <- function(parameters, model){
  mypars <- read.csv(file = parameters, sep = "\t")
  fixed <- mypars %>% filter(estimate == 0)
  constraints <- NULL
  if(nrow(fixed)>0){
    for(i in 1:length(fixed$parameterScale)) {
      parscale <- fixed$parameterScale[i]
      par <- fixed$parameterId[i] %>% as.character()
      value <- fixed$nominalValue[i]
      if(parscale == "lin") constraints <- c(constraints, value)
      if(parscale == "log10") constraints <- c(constraints, log10(value))
      if(parscale == "log") constraints <- c(constraints, log(value))
      else paste("This type of parameterScale is not supported.")
      names(constraints)[i] <- par
    } 
    parscales <- fixed$parameterScale %>% as.character()
    pars <- fixed$parameterId %>% as.character()
    names(parscales) <- pars
    attr(constraints,"parscale") <- parscales
  }
  estimated <- mypars %>% filter(estimate == 1)
  pouter <- NULL
  parlower <- NULL
  parupper <- NULL
  if(nrow(estimated)>0){
    for(i in 1:length(estimated$parameterScale)) {
      parscale <- estimated$parameterScale[i]
      par <- estimated$parameterId[i] %>% as.character()
      value <- estimated$nominalValue[i]
      lowervalue <- estimated$lowerBound[i]
      uppervalue <- estimated$upperBound[i]
      if(parscale == "lin"){
        pouter <- c(pouter, value)
        parlower <- c(parlower, lowervalue)
        parupper <- c(parupper, uppervalue)
      } else if(parscale == "log10"){
        pouter <- c(pouter, log10(value))
        parlower <- c(parlower, log10(lowervalue))
        parupper <- c(parupper, log10(uppervalue))
      } else if(parscale == "log"){
        pouter <- c(pouter, log(value))
        parlower <- c(parlower, log(lowervalue))
        parupper <- c(parupper, log(uppervalue))
      } else paste("This type of parameterScale is not supported.")
      names(pouter)[i] <- par
      names(parlower)[i] <- par
      names(parupper)[i] <- par
    } 
    parscales <- estimated$parameterScale %>% as.character()
    pars <- estimated$parameterId %>% as.character()
    names(parscales) <- pars
    attr(pouter,"parscale") <- parscales
    attr(pouter,"lowerBound") <- parlower
    attr(pouter,"upperBound") <- parupper
  }
  
  # check if additional parameters exist in SBML file
  model = readSBML(model)$getModel()
  n_pars <- model$getNumParameters()
  SBMLfixedpars <- NULL
  count <- 1
  for (i in 0:(n_pars-1)) {
    mypar <- model$getParameter(i)$getId()
    if(!mypar %in% names(pouter) & !mypar %in% names(constraints)){
      value <- model$getParameter(i)$getValue()
      SBMLfixedpars <- c(SBMLfixedpars, value)
      names(SBMLfixedpars)[count] <- mypar
      count <- count + 1
    }
  }
  
  return(list(constraints=constraints, pouter=pouter, SBMLfixedpars = SBMLfixedpars))
}


#' Import reactions from SBML. 
#' 
#' @description This function imports reactions from SBML. Reactions are written to an eqnlist object. 
#' Assignment rules for input functions are substituted in the rate column. Time is introduced as an additionel state t.
#'  
#' @param model SBML file as .xml
#'   
#' @return Eqnlist of reactions.
#'   
#' @author Marcus Rosenblatt and Svenja Kemmer
#'   
#' @export
#' 
getReactionsSBML <- function(model, conditions){
  m = readSBML(model)$getModel()
  
  # Initialization
  reactions <- NULL
  events <- NULL
  compartments <- NULL
  
  # import compartments
  N_species <- m$getNumSpecies()
  compartments <- do.call(c, lapply(0:(N_species-1), function(i){ m$getSpecies(i)$getCompartment()}))
  names_compartments <- do.call(c, lapply(0:(N_species-1), function(i){ m$getSpecies(i)$getId() }))
  names(compartments) <- names_compartments
  
  # if(unique(compartments)[1]=="default") compartments <- NULL
  
  # import reactions and adjust by means of compartments
  N_reactions <- m$getNumReactions()
  for (reaction in 0:(N_reactions-1)){
    Reactantstring <- ""
    Productstring <- ""
    eq <- m$getReaction(reaction)
    Reactantnr <- eq$getNumReactants(reaction)
    if(Reactantnr > 0) Reactantstring <- paste0( eq$getReactant(0)$getStoichiometry(), "*", eq$getReactant(0)$getSpecies())
    if(Reactantnr > 1) for (s in 1:(Reactantnr-1)) {
      Reactantstring <- paste0(Reactantstring, " + ",
                               paste0(eq$getReactant(s)$getStoichiometry(), "*", eq$getReactant(s)$getSpecies()))
    }
    Productnr <- eq$getNumProducts(reaction)
    if(Productnr > 0) Productstring <- paste0( eq$getProduct(0)$getStoichiometry(), "*", eq$getProduct(0)$getSpecies())
    if(Productnr > 1) for (s in 1:(Productnr-1)) {
      Productstring <- paste0(Productstring, " + ",
                              paste0(eq$getProduct(s)$getStoichiometry(), "*", eq$getProduct(s)$getSpecies()))
    }
    formula <- eq$getKineticLaw()$getFormula()
    if(str_detect(formula, "Function"))
      rate <- formula # to be double checked  # works for Borghans now
    else 
      rate <- gsub("pow", "", gsub(", ", "**", formula))
    #rate <- replaceOperation("pow", "**", eq$getKineticLaw()$getFormula())
    if(!is.null(compartments)){
      if(Reduce("|", str_detect(rate, unique(compartments)))){
        pos <- which(strsplit(rate, "")[[1]]=="*")[1]
        rate <- substr(rate,pos+1,length(strsplit(rate, "")[[1]]=="*"))
      }
    }
    reactions <- reactions %>% addReaction(Reactantstring, Productstring, rate)
  }
  reactions$rates <- gsub(" ","",reactions$rates)
  # import functions
  N_fundefs <- m$getNumFunctionDefinitions()
  if (N_fundefs > 0){
    for (fun in 0:(N_fundefs-1)){
      mymath <- m$getFunctionDefinition(fun)$getMath()
      # print(fun)
      string <- function_def_to_string(m$getFunctionDefinition(fun)) %>% gsub(" ","",.)
      # print(string)
      first <- strsplit(string, "=")[[1]][1]
      second <- strsplit(string, "=")[[1]][2]
      # print(first)
      # print(second)
      first <- gsub("\\(", "\\\\\\(", first)
      second <- gsub("\\(", "\\\\\\(", second)
      reactions$rates <- gsub(first,
                              second, reactions$rates)  # substitute m$getRule(0)$getVariable() by m$getRule(0)$getFormula()
      #print(formulaToL3String(mymath$getChild(mymath$getNumChildren()-1)))
    }
  }
  
  # import inputs
  N_rules <- m$getNumRules()
  if (N_rules > 0){
    for (rule in 0:(N_rules-1)){
      # substitute m$getRule(0)$getVariable() by m$getRule(0)$getFormula()
      reactions$rates <- replaceSymbols(m$getRule(rule)$getVariable(),
                                        paste0("(",m$getRule(rule)$getFormula(), ")"), reactions$rates)  
    }
  }
  
  reactions$rates <- gsub(" ","",reactions$rates)
  
  # replace function based inputs by events (done in reactions)
  for(fun in c("piecewise")){
    for(reaction in reactions$rates){
      if(str_detect(reaction, fun)){
        split <- str_split(reaction, fun)[[1]][2]
        count_bracket <- 0
        done <- F
        for(z in 1:nchar(split)){
          if(substr(split, z, z)=="(") count_bracket <- count_bracket+1
          if(substr(split, z, z)==")") count_bracket <- count_bracket-1
          if(count_bracket==0 & !done) {done <- T; pos <- z}
        }
        #pos <- which(strsplit(split, "")[[1]]==")")[2]
        event <- paste0(fun, substr(split, 1, pos))
        events <- c(events, event)
      }
    }
  }
  events <- unique(events)
  if(!is.null(events)) for(i in 1:length(events)){
    replace <- gsub("\\(", "\\\\\\(", events[i])
    replace <- gsub("\\*", "\\\\\\*", replace)
    replace <- gsub("\\+", "\\\\\\+", replace)
    #replace <- gsub("\\)", "\\\\\\)", replace)
    reactions$rates <- gsub(replace, paste0("event", i), reactions$rates)
    reactions <- reactions %>% addReaction("", paste0("event", i), "0")
  }
  
  # replace mathematical expressions 
  reactions$rates <- replaceSymbols(c("t", "TIME", "T"), "time", reactions$rates)
  
  TransformEvents <- function(events){
    if(!is.null(events)){
      do.call(rbind, lapply(1:length(events), function(i){
        myevent <- events[i]
        if(str_detect(myevent, "piecewise") & (str_detect(myevent, "leq") | str_detect(myevent, "lt"))){
          expr1 <- strsplit(myevent, ",")[[1]][2]
          expr1 <- gsub(paste0(strsplit(expr1, "\\(")[[1]][1],"\\("), "", expr1)
          expr2 <- strsplit(strsplit(myevent, ",")[[1]][3], ")")[[1]][1]
          if(expr1=="time") timepoint <- expr2 else 
            if(str_detect(expr1, "time-")) timepoint <- gsub("time-", "", expr1) else cat("Warning: Event not yet supported.")
          first <- strsplit(strsplit(myevent, "\\(")[[1]][2], ",")[[1]][1]
          second <- strsplit(strsplit(myevent, ",")[[1]][4], ")")[[1]][1]
          if(!is.na(suppressWarnings(as.numeric(timepoint)))) timepoint <- as.numeric(timepoint) # avoid warning if variable is not numeric
          if(!is.na(suppressWarnings(as.numeric(first)))) first <- as.numeric(first)
          if(!is.na(suppressWarnings(as.numeric(second)))) second <- as.numeric(second)
          return(data.frame(var=paste0("event",i), time=c(0,timepoint), value=c(first, second), method="replace"))
        } else {cat("Warning: Event not yet supported"); return(myevent)}
      }))
    } else return(NULL)
  }
  events <- TransformEvents(events)
  
  
  ## check for preequilibration conditions and handle them via events
  preeqEvents <- NULL
  myconditions <- read.csv(file = conditions, sep = "\t")
  myCons <- myconditions$conditionId
  mypreeqCons <- NULL
  attrib <- NULL
  for (con in myCons){if(paste0("preeq_", con)%in%myCons) mypreeqCons <- c(mypreeqCons, con)}
  if(!is.null(mypreeqCons)){
    for (con in mypreeqCons){
      mycongrid <- filter(myconditions, conditionId==con | conditionId==paste0("preeq_", con))
      if(ncol(mycongrid)>1){
        for(i in 2:ncol(mycongrid)){
          preeqEvents <- addEvent(preeqEvents, var=names(mycongrid)[i], time=0, value=mycongrid[[which(mycongrid$conditionId==con),i]], method="replace")
          attrib <- c(attrib, mycongrid[[which(mycongrid$conditionId==paste0("preeq_",con)),i]])
        }
      }
    }
  }
  mystates <- reactions$states
  reactions_orig <- reactions
  attr(preeqEvents, "initials") <- attrib
  if(!is.null(preeqEvents)) for(i in 1:nrow(preeqEvents)){
    events <- rbind(events, preeqEvents[i,])
    reactions <- reactions %>% addReaction("", preeqEvents[[i,"var"]], "0")
  }
  
  
  # for(i in 1:length(reactions$rates)){
  #   reaction <- reactions$rates[i]
  #   if(str_detect(reaction, "pow")){
  #     reaction_new <- gsub("pow", "", reaction)
  #     # reaction_new <- gsub(", ", "**", reaction_new)
  #     reaction_new <- gsub(",", "**", reaction_new)
  #     reactions$rates[i] <- reaction_new
  #   } else reactions$rates[i] <- reaction
  # }
  
  mydata <- as.data.frame(reactions)
  reactions <- as.eqnlist(mydata, compartments)
  
  return(list(reactions=reactions, events=events, reactions_orig=reactions_orig, preeqEvents=preeqEvents, mystates=mystates))
}


## function from Frank
function_def_to_string <- function(fun)
{
  if (is.null(fun)) return;
  id <- fun$getId()
  
  math <- fun$getMath()
  if (is.null(math)) return;
  
  # the function will be of the form lambda(a, b, formula)
  # so the last child of the AST_Node is the actual math everything in front of it 
  # the arguments
  #
  # So to generate the desired output function_definition_id(args) = math
  # we use: 
  
  
  num_children <- math$getNumChildren()
  
  result <- paste(id, '(', sep="")
  for (i in 1:(num_children-1) ) # this leaves out the last one
  {
    result <- paste(result, libSBML::formulaToL3String(math$getChild(i-1)), sep="")
    
    if (i < (num_children-1))
      result <- paste(result, ', ', sep="")
    
  }
  
  result <- paste(result, ') = ', libSBML::formulaToL3String(math$getChild(num_children - 1)), sep="")
  
  return (result)
  
}



