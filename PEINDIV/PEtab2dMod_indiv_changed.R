#' Rename parscales to the names needed in the base trafo
#'
#' @param parscales setNames(PETABPars$parameterScale, PETABPars$parameterId)
#' @param est.grid data.table
#'
#' @return
#' @export
#'
#' @examples
updateParscalesToBaseTrafo <- function(parscales, est.grid) {
  # Get name mapping between est.grid pars and outer pars
  parseg_outer <- getEstGridParameterMapping(est.grid)
  parsouter_eg <- setNames(names(parseg_outer), parseg_outer)
  # Get scales for outer pars
  parscales_eg <- copy(parscales)
  names(parscales_eg) <- parsouter_eg[names(parscales)]
  
  # Determine if there are any duplicated outer pars with different scales.
  dupes <- names(parscales_eg)[duplicated(names(parscales_eg))]
  dupes <- unique(dupes)
  for (d in dupes) if (length(unique(parscales_eg[d])) > 1) 
    stop("The following parameter refers to the same structural model parameter, but has different ",
         "scales in different conditions. This is not allowed. \n",
         "Parameter: ", d , "\n",
         "Outer pars: ", paste0(names(unique(parscales_eg[d])), collapse = ", "))
  
  # If all went fine, remove the duplicates and return updated parscales
  parscales_eg[!duplicated(names(parscales_eg))]
  parscales_eg
}



#' #' Determine symbolic trafos in a vector
#' #'
#' #' @param trafo_string 
#' #'
#' #' @return vector TRUE/FALSE
#' #' @export
#' #' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' #' @md
#' #'
#' #' @examples
#' #' trafo_string <- c(STAT5A = "207.6 * ratio", STAT5B = "207.6 - 207.6 * ratio", 
#' #' pApB = "0", pApA = "0", pBpB = "0", nucpApA = "0", nucpApB = "0", 
#' #' nucpBpB = "0")
#' #' is_symbolic_trafo(trafo_string)
#' is_symbolic_trafo <- function(trafo_string){
#'   is_numeric <- !is.na(as.numeric(trafo_string))
#'   is_character <- vapply(trafo_string, function(ts) identical(getSymbols(ts),names(ts)), TRUE)
#'   is_symbolic <- !is_numeric & !is_character
#' }


#' Title
#'
#' @param trafo_string 
#'
#' @return
#' @export
#'
#' @examples
#' trafo_string <- c(STAT5A = "207.6 * ratio", STAT5B = "207.6 - 207.6 * ratio", 
#'                   pApB = "alpha", pApA = "beta", pBpB = "pBpB", nucpApA = "0.1", nucpApB = "0", 
#'                   nucpBpB = "0")
#' getTrafoType(trafo_string)
getTrafoType <- function(trafo_string) {
  vapply(names(trafo_string), function(nm) {
    ts <- trafo_string[nm]
    pd <- getParseData(parse(text = ts))
    if (nrow(pd) > 2) return("TRAFO")
    if (pd[1,"token"] == "SYMBOL") {
      if (pd[1,"text"] == nm) return("SYMBOL")
      else return("TRAFO")
      }
    if (pd[1,"token"] == "NUM_CONST") return("NUMBER")
    stop("Unkown trafo type: ", ts)
  }, FUN.VALUE = "TYPE")
}

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
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de) building on the original function of Marcus and Svenja
#' @md
#'
#' @export
#' 
#' # setwd(rstudioapi::getActiveProject())
#' # devtools::load_all()
#' # f <- list.files("BenchmarkModels")
#' # modelname = f[1]
#' # path2model = "BenchmarkModels/"
#' # testCases = FALSE
#' # path2TestCases = "PEtabTests/"
#' # compile = TRUE
#' # SBML_file = NULL
#' # observable_file = NULL
#' # condition_file = NULL
#' # data_file = NULL
#' # parameter_file = NULL

# >>>> comment out <<<<<<<<<<< ----
setwd(rstudioapi::getActiveProject())
devtools::load_all()
f <- list.files("PEtabTests/")
i <- 2

modelname = f[i]
path2model = "BenchmarkModels/"
testCases = TRUE
path2TestCases = "PEtabTests/"
compile = TRUE
SBML_file = NULL
observable_file = NULL
condition_file = NULL
data_file = NULL
parameter_file = NULL
# >>>> comment out <<<<<<<<<<< ----

importPEtabSBML_indiv <- function(modelname = "Boehm_JProteomeRes2014",
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
  ## Read previously imported files --------------------
  if(is.null(modelname)) modelname <- "mymodel"
  
  mywd <- getwd()
  dir.create(paste0(mywd,"/CompiledObjects/"), showWarnings = FALSE)
  
  setwd(paste0(mywd,"/CompiledObjects/"))
  if(compile == FALSE & file.exists(paste0(modelname,".rds"))){
    petab <- readRDS(paste0(modelname,".rds"))
    loadDLL(petab$obj_data)
    return(petab)
  }
  setwd(mywd)
  
  ## load required packages
  require(libSBML) # => Not very nice, better explicitly import the required functions
  
  ## Define path to SBML and PEtab files --------------------
  if(testCases == FALSE){
    if(is.null(SBML_file))       SBML_file       <- paste0(path2model, modelname, "/model_", modelname, ".xml")
    if(is.null(observable_file)) observable_file <- paste0(path2model, modelname, "/observables_", modelname, ".tsv")
    if(is.null(condition_file))  condition_file  <- paste0(path2model, modelname, "/experimentalCondition_", modelname, ".tsv")
    if(is.null(data_file))       data_file       <- paste0(path2model, modelname, "/measurementData_", modelname, ".tsv")
    if(is.null(parameter_file))  parameter_file  <- paste0(path2model, modelname, "/parameters_", modelname, ".tsv")
  } else{
    SBML_file       <- paste0(path2TestCases, modelname, "/_model.xml")
    observable_file <- paste0(path2TestCases, modelname, "/_observables.tsv")
    condition_file  <- paste0(path2TestCases, modelname, "/_conditions.tsv")
    data_file       <- paste0(path2TestCases, modelname, "/_measurements.tsv")
    parameter_file  <- paste0(path2TestCases, modelname, "/_parameters.tsv")
  }
  if(!file.exists(SBML_file)){       cat(paste0("The file ",mywd,SBML_file, " does not exist. Please check spelling or provide the file name via the SBML_file argument.")); return(NULL)}
  if(!file.exists(observable_file)){ cat(paste0("The file ",mywd,observable_file, " does not exist. Please check spelling or provide the file name via the observable_file argument.")); return(NULL)}
  if(!file.exists(condition_file)){  cat(paste0("The file ",mywd,condition_file, " does not exist. Please check spelling or provide the file name via the condition_file argument.")); return(NULL)}
  if(!file.exists(data_file)){       cat(paste0("The file ",mywd,data_file, " does not exist. Please check spelling or provide the file name via the data_file argument.")); return(NULL)}
  if(!file.exists(parameter_file)){  cat(paste0("The file ",mywd,parameter_file, " does not exist. Please check spelling or provide the file name via the parameter_file argument.")); return(NULL)}
  
  # .. Read grids to see original fiile structure => Delete ----
  # tsv_observable <- fread(observable_file)
  # tsv_condition  <- fread(condition_file)
  # tsv_data       <- fread(data_file)
  # tsv_parameter  <- fread(parameter_file)
  
  ## Model Definition - Equations --------------------
  mylist           <- getReactionsSBML(SBML_file, condition_file)
  myreactions      <- mylist$reactions
  myreactions_orig <- mylist$reactions_orig
  myevents         <- mylist$events
  mypreeqEvents    <- mylist$preeqEvents
  mystates         <- mylist$mystates
  myobservables    <- getObservablesSBML(observable_file)
  
  ## Get Data ------------
  mydataSBML <- getDataPEtabSBML(data_file, observable_file)
  mydata     <- mydataSBML$data
  myerrors   <- mydataSBML$errors
  myerr <- NULL
  
  ## Define constraints, initials, parameters and compartments --------------
  myparameters   <- getParametersSBML(parameter_file, SBML_file)
  # [ ] Check constraints
  myconstraints  <- myparameters$constraints
  SBMLfixedpars  <- myparameters$SBMLfixedpars
  myfit_values   <- myparameters$pouter
  myinitialsSBML <- getInitialsSBML(SBML_file, condition_file)
  mycompartments <- myinitialsSBML$compartments
  myinitials     <- myinitialsSBML$initials
  # set remaining event initials to 0
  inits_events <- setdiff(unique(myevents$var), unique(mypreeqEvents$var))
  inits_events <- setNames(rep(0, length(inits_events)), inits_events)
  pars_est <- setNames(nm = names(myfit_values))
  
  
  
  ## Parameter transformations -----------
  # .. Generate condition.grid -----
  grid <- getConditionsSBML(conditions = condition_file, data = data_file)
  mypreeqCons      <- grid$preeqCons
  mycondition.grid <- grid$condition_grid
  attr(mydata, "condition.grid") <- mycondition.grid
  # .. Build fix.grid, est.grid and trafo -----
  
  # Copy condition.grid, take unique identifying column only
  cg <- data.table::data.table(mycondition.grid)
  cg <- cg[,!"conditionName"]
  data.table::setnames(cg, "conditionId", "condition")
  # Initialize fix.grid and est.grid
  # Determine which columns contain values and/or parameter names
  is_string  <- vapply(cg[,-1], function(x) any(is.na(as.numeric(x))), FUN.VALUE = TRUE)
  is_string  <- which(is_string) + 1
  is_numeric <- vapply(cg[,-1], function(x) any(!is.na(as.numeric(x))), FUN.VALUE = TRUE)
  is_numeric <- which(is_numeric) + 1
  # For mixed columns (string & numeric), need NA in the respective places
  fix.grid <- data.table::copy(cg)
  fix.grid <- fix.grid[,.SD,.SD = c(1, is_numeric)]
  fix.grid[,(names(fix.grid)[-1]) := lapply(.SD, as.numeric), .SDcols = -1]
  fix.grid[,`:=`(ID = 1:.N)]
  est.grid <- data.table::copy(cg)
  est.grid <- est.grid[,.SD,.SD = c(1, is_string)]
  est.grid[,(names(est.grid)[-1]) := lapply(.SD, function(x) {replace(x, !is.na(as.numeric(x)), NA)}), .SDcols = -1]
  est.grid[,`:=`(ID = 1:.N)]
  gridlist <- list(est.grid = est.grid, fix.grid = fix.grid)
  
  trafo <- setNames(nm = unique(c(getParameters(myreactions), getSymbols(myobservables), getSymbols(myerrors), getSymbols(as.character(myevents$value)))))
  trafo <- trafo[trafo != "time"]
  
  # .. Fill values into grid and trafo -----
  parameterlist <- list(
    list(par = SBMLfixedpars, overwrite = FALSE),
    list(par = mycompartments, overwrite = FALSE),
    list(par = myinitials, overwrite = FALSE),
    list(par = myconstraints, overwrite = FALSE),
    list(par = inits_events, overwrite = FALSE),
    list(par = pars_est, overwrite = FALSE)
  )
  for (pl in parameterlist) {
    par <- pl$par
    trafoType <- getTrafoType(par)
    is_symbolic_trafo <- trafoType %in% c("TRAFO")
    par_grid <- par[!is_symbolic_trafo]
    par_traf <- par[ is_symbolic_trafo]
    if (any(!is_symbolic_trafo)) gridlist <- add_pars_to_grids(pars = par_grid, gridlist = gridlist, FLAGoverwrite = pl$overwrite)
    if (any( is_symbolic_trafo)) trafo <- repar("x~y", trafo, x = names(par_traf), y = par_traf)
  }
  
  # [ ] Pre-Equi Events
  # branch trafo for different conditions
  # # set preequilibration event initials to corresponding values
  # if(!is.null(mypreeqEvents)){
  #   mypreeqEvents2replace <- filter(mypreeqEvents, !var%in%mystates)
  #   mytrafoL <- repar("x~y", mytrafoL , x = unique(mypreeqEvents2replace$var), y = attr(mypreeqEvents2replace, "initials"))
  # }
  # [ ] Need example for preeqEvents
  
  # .. Parscales -----
  # adjust symbolic trafo
  parscales <- attr(myfit_values,"parscale")
  parscales <- updateParscalesToBaseTrafo(parscales, gridlist$est.grid)
  # if (length(nm <- setdiff(names(parscales), getSymbols(trafo)))) 
  #   stop("undefined parameters in trafo: ", paste0(nm, collapse = ", "))
  trafo <- repar("x ~ 10**(x)", trafo = trafo, x = names(which(parscales=="log10")))
  trafo <- repar("x ~ exp(x)" , trafo = trafo, x = names(which(parscales=="log")))
  
  # Adjust fix.grid
  fg <- gridlist$fix.grid
  for (nm in intersect(names(fg), names(parscales))) {
    scale <- parscales[nm]
    if (scale == "log10") fg[[nm]] <- log10(fg[[nm]])
    if (scale == "log") fg[[nm]] <- log(fg[[nm]])
    fg[[nm]][!is.finite(fg[[nm]])] <- -1000
  }
  gridlist$fix.grid <- fg
  
  
  # -------------------------------------------------------------------------#
  # Model Compilation ----
  # -------------------------------------------------------------------------#
  setwd(paste0(mywd,"/CompiledObjects/"))
  myg <- Y(myobservables, myreactions, compile=TRUE, modelname=paste0("g_",modelname))
  
  myodemodel <- odemodel(myreactions, forcings = NULL, events = myevents, fixed=NULL,
                         estimate = getParametersToEstimate(est.grid = gridlist$est.grid,
                                                            trafo = trafo,
                                                            reactions = myreactions),
                         modelname = paste0("odemodel_", modelname),
                         jacobian = "inz.lsodes", compile = TRUE)
  
  tolerances <- 1e-7
  myx <- Xs(myodemodel,
            optionsOde = list(method = "lsoda", rtol = tolerances, atol = tolerances, maxsteps = 5000),
            optionsSens = list(method = "lsodes", lrw=200000, rtol = tolerances, atol = tolerances))
  
  if(length(getSymbols(myerrors))){
    setwd(paste0(mywd,"/CompiledObjects/"))
    myerr <- Y(myerrors, f = c(as.eqnvec(myreactions), myobservables), states = names(myobservables), attach.input = FALSE, compile = TRUE, modelname = paste0("errfn_", modelname))
    setwd(mywd)
  }
  
  # [ ] Pre-equilibration
  mypSS <- Id()
  
  myp <- P(trafo, compile = TRUE, modelname = paste0("P_", modelname))
  
  setwd(mywd)
  
  # .. Collect lists -----
  symbolicEquations <- list(
    reactions = myreactions,
    observables = myobservables,
    errors  = myerrors,
    trafo = trafo)
  fns <- list(
    g = myg,
    x = myx,
    p1 = mypSS, # [ ] Pre-Equilibration
    p0 = myp
  )
  
  # .. Generate high-level fns -----
  prd <- PRD_indiv(prd0 = Reduce("*", fns), est.grid = gridlist$est.grid, fix.grid = gridlist$fix.grid)
  obj_data <- normL2_indiv(mydata, Reduce("*", fns), errmodel = myerr,
                           est.grid = gridlist$est.grid, fix.grid = gridlist$fix.grid,
                           times = seq(0,max(as.data.frame(mydata)$time), len=501))
  # .. Collect final list -----
  petab <- list(
    symbolicEquations = symbolicEquations,
    odemodel = myodemodel,
    # [ ] complete data specification could be lumped: data, gridlist, myerr
    data = mydata,
    gridlist = gridlist,
    e = myerr,
    fns = fns,
    prd = prd,
    obj_data = obj_data,
    pars = myfit_values
  )
  
  # .. Save everything -----
  setwd(paste0(mywd,"/CompiledObjects/"))
  saveRDS(petab, paste0(modelname, ".rds"))
  setwd(mywd)
  
  petab
  
}

# >>>> from here: comment out <<<<<<<<<<< ----

# -------------------------------------------------------------------------#
# Testing ----
# -------------------------------------------------------------------------#

# setwd(rstudioapi::getActiveProject())
# devtools::load_all()
# f <- list.files("BenchmarkModels")
# 
# 
# petab <- importPEtabSBML(modelname = f[3],
#                          path2model = "BenchmarkModels/",
#                          testCases = FALSE,
#                          path2TestCases = "PEtabTests/",
#                          compile = TRUE,
#                          SBML_file = NULL,
#                          observable_file = NULL,
#                          condition_file = NULL,
#                          data_file = NULL,
#                          parameter_file = NULL)
# 
# p <- petab$fns$p0
# x <- petab$fns$x
# times <- seq(0,max(as.data.frame(petab$data)$time), len=501)
# pred <- petab$prd(times, petab$pars, FLAGbrowserN = 1)
# plotCombined(pred, petab$data)


# -------------------------------------------------------------------------#
# Test models ----
# -------------------------------------------------------------------------#
setwd(rstudioapi::getActiveProject())
devtools::load_all()
f <- list.files("PEtabTests/")
i <- 2
# ..  -----
# debugonce(importPEtabSBML_indiv)
petab <- importPEtabSBML_indiv(modelname = f[i],
                         path2model = "BenchmarkModels/",
                         testCases = TRUE,
                         path2TestCases = "PEtabTests/",
                         compile = TRUE,
                         SBML_file = NULL,
                         observable_file = NULL,
                         condition_file = NULL,
                         data_file = NULL,
                         parameter_file = NULL)

p <- petab$fns$p0
x <- petab$fns$x
times <- seq(0,max(as.data.frame(petab$data)$time), len=501)
pred <- petab$prd(times, petab$pars, FLAGbrowserN = 1)
plotCombined(pred, petab$data)
i <- i+1

# -------------------------------------------------------------------------#
#  ----
# -------------------------------------------------------------------------#


# prd(times, myfit_values, FLAGbrowser = 1)
# prd(times, myfit_values, FLAGbrowser = 2)

# p <- P_indiv(myp, est.grid = gridlist$est.grid, fix.grid = gridlist$fix.grid)
# wup <- p(myfit_values)
# wup
# myfit_values

# myp(myfit_values)
# rp <- tempfile()
# Rprof(rp)
# obj_data(myfit_values) 
# Rprof(NULL)
# summaryRprof(rp)
# pv <- profvis::profvis(prof_input = rp); htmlwidgets::saveWidget(pv, paste0(rp, ".html")); browseURL(paste0(rp, ".html"))

# debugonce(obj_data)
lapply(1:10,function(i)petab$obj_data(petab$pars))

parallel::mclapply(1:12, function(i) obj_data(myfit_values), mc.cores = 4)



# obj_data(myfit_values, FLAGbrowser = T)

b1 <- rbenchmark::benchmark(petab$obj_data(petab$pars), replications = 20)

b1.2 <- rbenchmark::benchmark(mclapply(1:36, function(i) petab$obj_data(petab$pars), mc.cores= 12, mc.preschedule = TRUE), 
                              replications = 3)



importPEtabSBML(modelname, path2model)
# debugonce(obj)
# obj(pouter)
# 
# rp <- tempfile()
# Rprof(rp)
# obj(pouter)
# Rprof(NULL)
# summaryRprof(rp)
# pv <- profvis::profvis(prof_input = rp); htmlwidgets::saveWidget(pv, paste0(rp, ".html")); browseURL(paste0(rp, ".html"))

b2 <- rbenchmark::benchmark(obj(pouter), replications = 20)
b2.2 <- rbenchmark::benchmark(mclapply(1:36, function(i) obj(pouter), mc.cores= 12), 
                              replications = 3)

writeLines(capture.output(print(list(
  b1,
  b1.2,
  b2,
  b2.2
))), "~/wup.txt")


# -------------------------------------------------------------------------#
# p in R ----
# -------------------------------------------------------------------------#
p <- P(c("a"  ="exp(b)"), compile = TRUE, modelname = "p")

p(c(b = -1000))
p(c(b = -1000)) %>% getDerivs()

# Exit ----
