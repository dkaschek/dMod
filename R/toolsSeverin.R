#' Run an any R function on a remote HPC system with SLURM
#' 
#' @description Generates R and bash scrips, transfers them to remote via ssh.
#' 'sshpass' needs to be installed on your local machine to circumvent password entry.
#' @details \code{distribute_computing()} generates R and bash scrips designed to run
#' on a HPC system managed by the SLURM batch manager. The current workspace 
#' together with the scripts are exported and transferred to remote via SSH. If
#' ssh-key authentication is not possible, the ssh-password can be passed and
#' is used by sshpass (which has to be installed on the local machine).
#' 
#' The code to be executed remotely is passed to the \code{...} argument, its final
#' output is stored in \code{cluster_result}, which is loaded in the local
#' workspace by the \code{get()} function.
#' 
#' It is possible to either run repetitions of the same program realization (by 
#' use of the \code{no_rep} parameter), or to pass a list of parameter arrays
#' via \code{var_values}. The parameters to be changed for every run *must* be
#' named \code{var_i} where i corresponds to the i-th array in the \code{var_values}
#' parameter.
#' 
#' @param ... R code to be remotely executed, parameters to be changed for runs
#' must be named \code{var_i} see "Details".
#' @param jobname Name in characters for the run must be unique, existing will be overwritten.
#' Must not contain the string "Minus".
#' @param partition Define the partition used for the SLURM manager. Default is
#' "single". Should only be changed if explicitly wanted.
#' @param cores Number of cores per node to be used. If a value is set to 16 < 
#' value < 25, the number of possible nodes isseverely limited. 16 or less is
#' possible on all nodes, more on none.
#' @param nodes Nodes per task. Default is 1, should not be changed since 
#' \code{distributed_computing()} is set up to be regulated by number of repititions.
#' @param walltime Estimated runtime in the format \code{hh:mm:ss}, default is 1h.
#' Jobs will be canceled after the defined time.
#' @param ssh_passwd To be set when the sshpass should be used to automatically 
#' authenticate on remote via passphrase. This is an obvious security nightmare...
#' @param machine SSH address in the form \code{user@@remote_location}.
#' @param var_values List of parameter arrays. The number of arrays (i.e. list 
#' entries) must correspond to the number of parameters in the function passed 
#' to \code{...}. These parameters must be named \code{var_i} where the i must 
#' replaced by the indes of the corresponding array in \code{var_values}. The 
#' length of the arrays define the number of nodes used with the j-th node use
#' the j-th entry of the arrays for the corresponding \code{var_i}. If \code{no_rep}
#' is used, \code{var_valus} must be set to \code{NULL}.
#' @param no_rep Number of repetitions. Usage of this parameter and \code{var_values}
#' are mutual exclusive. When used the function passed to \code{...} is executed
#' \code{no_rep} times simultaneously using one node per realization. 
#' @param recover Logical parameter, if set to \code{TRUE} nothing is calculated, the
#' functions \code{check()}, \code{get()} and \code{purge()} can be used on the
#' results generated previously under the same \code{jobname}.
#' @param purge_local Logical, if set to \code{TRUE} the \code{purge()} function also removes
#' local files.
#' @param compile Logical, if set to \code{TRUE} the source files are transferred
#' to remote and compiled there. If set to \code{FALSE}, the local shared objects
#' are transferred and used instead.
#' @param custom_folders named vector with exact three entries named 'compiled',
#' 'output' and 'tmp'. The values are strings with relative paths from the current
#' working directory to the respective directory of the compiled files, the temporary
#' folder from which files will be copied to the cluster and the output folder in
#' which the calculated result from the cluster will be saved.
#' The default is \code{NULL}, then everything is done from the current working directory.
#' If only a subset of the folders should be changed, all other need to be set to
#' \code{./}.
#' @param resetSeeds logical, if set to \code{TRUE} (default) the parameter vector
#' with random seeds \code{.Random.seeds} from the transferred work space is deleted
#' on remote. This ensures that each node has uses a different set of  (pseudo) random
#' numbers. Set to FALSE at own risk.
#' @param returnAll logical if set to \code{TRUE} (default) all results are returned, if set
#' to \code{FALSE} only the \code{*result.RData} files are returned.
#'  
#' @return List of functions \code{check()}, \code{get()} and \code{purge()}. 
#' \code{check()} checks, if the result is ready. 
#' \code{get()} copies all files from the remote working directory to local and 
#' loads all present results (even if not all nodes where done) in the current 
#' active workspace as the object \code{cluster_result} which is a list with the
#' results of each node as entries.
#' \code{purge()} deletes the temporary folder on remote and if \code{purge_local}
#' is set to \code{TRUE}.
#' @examples
#' \dontrun{
#' out_distributed_computing <- distributed_computing(
#' {
#'   mstrust(
#'     objfun=objective_function,
#'     center=outer_pars,
#'     studyname = "study",
#'     rinit = 1,
#'     rmax = 10,
#'     fits = 48,
#'     cores = 16,
#'     iterlim = 700,
#'     sd = 4
#'   )
#' },
#' jobname = "my_name",
#' partition = "single",
#' cores = 16,
#' nodes = 1,
#' walltime = "02:00:00",
#' ssh_passwd = "password",
#' machine = "cluster",
#' var_values = NULL,
#' no_rep = 20,
#' recover = F,
#' compile = F
#' )
#' out_distributed_computing$check()
#' out_distributed_computing$get()
#' out_distributed_computing$purge()
#' result <- cluster_result
#' print(result)
#' 
#' 
#' # calculate profiles
#' var_list <- profile_pars_per_node(best_fit, 4)
#' profile_jobname <- paste0(fit_filename,"_profiles_opt")
#' method <- "optimize"
#' profiles_distributed_computing <- distributed_computing(
#'   {
#'     profile(
#'       obj = obj,
#'       pars =  best_fit,
#'       whichPar = (as.numeric(var_1):as.numeric(var_2)),
#'       limits = c(-5, 5),
#'       cores = 16,
#'       method = method,
#'       stepControl = list(
#'         stepsize = 1e-6,
#'         min = 1e-4, 
#'         max = Inf, 
#'         atol = 1e-2,
#'         rtol = 1e-2, 
#'         limit = 100
#'       ),
#'       optControl = list(iterlim = 20)
#'     )
#'   },
#'   jobname = profile_jobname,
#'   partition = "single",
#'   cores = 16,
#'   nodes = 1,
#'   walltime = "02:00:00",
#'   ssh_passwd = "password",
#'   machine = "cluster",
#'   var_values = var_list,
#'   no_rep = NULL,
#'   recover = F,
#'   compile = F
#' )
#' profiles_distributed_computing$check()
#' profiles_distributed_computing$get()
#' profiles_distributed_computing$purge()
#' profiles  <- NULL
#' for (i in cluster_result) {
#'   profiles <- rbind(profiles, i)
#' }
#' }
#' 
#' @export
distributed_computing <- function(
  ...,
  jobname,
  partition = "single",
  cores = 16,
  nodes = 1,
  walltime = "01:00:00",
  ssh_passwd = NULL,
  machine = "cluster",
  var_values = NULL,
  no_rep = NULL,
  recover = T,
  purge_local = F,
  compile = F,
  custom_folders = NULL,
  resetSeeds = TRUE,
  returnAll = TRUE
  # called_function = "func(a = var_1, b = var_1, name = jobname,id = 01)",
) {
  original_wd <- getwd()
  if (is.null(custom_folders)) {
    output_folder_abs <- "./"
  } else if(!is.null(custom_folders) & !all(length(custom_folders) == 3 & sort(names(custom_folders)) == c("compiled", "output", "tmp"))) {
    warning("'custom_folders' must be named vector with exact three elements:\n
            'compiled', 'output', 'tmp', containing relative paths to the resp folders\n
            input is wrong, ignored.\n")
  } else {
    compiled_folder <- custom_folders["compiled"]
    output_folder <- custom_folders["output"]
    tmp_folder <- custom_folders["tmp"]
    
    system(paste0("cp ", compiled_folder, "* ", tmp_folder))
    
    setwd(output_folder)
    output_folder_abs <- getwd()
    
    setwd(original_wd)
    setwd(tmp_folder)
  }
  
  on.exit(setwd(original_wd))
  
  # - definitions - #
  
  # relative path to the working directory, will now allways be used
  
  wd_path <- paste0("./",jobname, "_folder/")
  data_path <- paste0(getwd(),"/",jobname, "_folder/")
  
  # number of repetitions
  if(!is.null(no_rep) & is.null(var_values)) {
    num_nodes <- no_rep - 1
  } else if(is.null(no_rep) & !is.null(var_values)) {
    num_nodes <- length(var_values[[1]]) - 1
  } else {
    stop("I dont know what you want how often done. Please set either 'no_rep' or pass 'var_values' (_not_ both!)")
  }
  
  # define the ssh command depending on 'sshpass' being used
  if(is.null(ssh_passwd)){
    ssh_command <- "ssh "
    scp_command <- "scp "
  } else {
    ssh_command <- paste0("sshpass -p ", ssh_passwd, " ssh ")
    scp_command <- paste0("sshpass -p ", ssh_passwd, " scp ")
  }
  
  # - output functions - #
  # Structure of the output 
  out <- structure(vector("list", 3), names = c("check", "get", "purge"))
  
  # check function
  out[[1]] <- function() {
    
    result_length <- length(
      suppressWarnings(
        system(
          paste0(ssh_command, machine, " 'ls ", jobname, "_folder/ | egrep *result.RData'"),
          intern = TRUE)
      )
    )
    
    if (result_length == num_nodes +1) {
      cat("Result is ready!\n")
      return(TRUE)
    }
    else if (result_length  < num_nodes +1) {
      cat("Result from", result_length, "out of", (num_nodes +1), "nodes are ready.")
      return(FALSE)
    }
    setwd(original_wd)
    
  }
  
  # get function
  out[[2]] <- function () {
    # copy all files back
    if (returnAll == T) {
      system(
        paste0(
          "mkdir -p ", output_folder_abs, "/", jobname,"_folder/results/; ",
          ssh_command, "-n ", machine, # go to remote
          " 'tar -C ", jobname, "_folder", " -czf - ./'", # compress all files on remote)
          " | ", # pipe to local
          "",
          "tar -C ", output_folder_abs, "/", jobname,"_folder/results/ -xzf -"
        )
      )
    } else {
      # copy only result files back
      system(
        paste0(
          "mkdir -p ", output_folder_abs, "/", jobname,"_folder/results/; ",
          ssh_command, "-n ", machine, " '",
          "find ", jobname, "_folder -type f -name \"*result.RData\" -exec tar -czf - {} +'",
          " | tar --strip-components=1 -xz -C ", shQuote(paste0(output_folder_abs, "/", jobname,"_folder/results/"))
        )
      )
    }

    
    # get list of all currently available output files
    # setwd(paste0(jobname,"_folder/results"))
    result_list <- structure(vector(mode = "list", length = num_nodes+1))
    result_files <- list.files(path=paste0(output_folder_abs,"/",jobname,"_folder/results/"),pattern = glob2rx("*result.RData"))
    # setwd("../../")
    
    # result_files <- Sys.glob(file.path(paste0(wd_path, "/results/*RData")))
    
    for (i in seq(1, length(result_files))) {
      cluster_result <- NULL
      check <- try(load(file = paste0(output_folder_abs,"/",jobname,"_folder/results/",result_files[i])), silent = TRUE) 
      if (!inherits("try-error", check)) result_list[[i]] <- cluster_result
    }
    
    if (length(result_files) != num_nodes +1) {
      cat("\n\tNot all results ready\n")
    }
    setwd(original_wd)
    # results_cluster <- Sys.glob(paste0(wd_path, "*.RData")) %>% map_dfr(load)
    .GlobalEnv$cluster_result <- result_list
  }
  
  
  
  # purge function
  out[[3]] <- function (purge_local = FALSE) {
    # remove files remote
    system(
      paste0(
        ssh_command, machine, " rm -rf ", jobname, "_folder"
      )
    )
    # also remove local files if want so
    if (purge_local) {
      system(
        paste0("rm -rf ", output_folder_abs, "/", jobname,"_folder")
      )
    }
    setwd(original_wd)
  }
  
  # if recover == T, stop here
  if(recover) {
    setwd(original_wd)
    return(out)
  } 
  
  
  
  
  
  
  
  
  
  # - calculation - #
  # create wd for this run
  system(
    paste0("rm -rf ",jobname,"_folder")
  )
  system(
    paste0("mkdir ",jobname,"_folder/")
  )
  
  
  
  # export current workspace
  save.image(file = paste0(wd_path,jobname, "_workspace.RData")) 
  
  
  # WRITE R
  
  
  # generate list of currently loaded packages
  package_list <- sapply(strsplit(search(), "package:", fixed = TRUE), function(v) v[2])
  package_list <- package_list[!is.na(package_list)]
  package_list <- paste(paste0("try(library(", package_list, "))"), collapse = "\n")
  if (compile) {
    load_so <- paste0("dyn.load('",jobname,"_shared_object.so')")
  } else (
    load_so <- ""
  )
  
  
  
  # generate parameter lists
  if (!is.null(var_values)) {
    var_list <- paste(
      lapply(
        seq(1,length(var_values)),
        function(i) {
          if( class(var_values[[i]]) == "character") {
            paste0("var_values_", i, "=c('", paste(var_values[[i]], collapse="','"),"')")
          } else {
            paste0("var_values_", i, "=c(", paste(var_values[[i]], collapse=","),")")
          }
          
          
        }
      ),
      collapse = "\n"
    )
    # cat(variable_list)
    
    # List of all names of parameters that will be changes between runs
    var_names <- paste(lapply(seq(1,length(var_values)), function(i) paste0("var_",i)))
    
    # Variables per run
    var_per_run <- paste(
      lapply(
        seq(1, length(var_values)),
        function(i) {
          paste0("var_", i, "=var_values_",i,"[(as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID')) + 1)]")
        }
      ),
      collapse = "\n"
    )
  } else {
    var_list <- ""
    var_per_run <- ""
  }
  
  # cat(var_per_run)
  
  
  # define fixed pars
  fixedpars <- paste(
    "node_ID = Sys.getenv('SLURM_ARRAY_TASK_ID')",
    "job_ID = Sys.getenv('SLURM_JOB_ID')",
    paste0("jobname = ", "'",jobname,"'"),
    sep = "\n"
  )
  
  
  # WRITE R
  expr <- as.expression(substitute(...))
  # expr <- as.expression(substitute(called_function))
  cat(
    paste(
      "#!/usr/bin/env Rscript",
      "",
      "# Load packages",
      package_list,
      "try(library(tidyverse))",
      "",
      "# Load environment",
      paste0("load('",jobname,"_workspace.RData')"),
      "",
      if (resetSeeds == TRUE & exists(".Random.seed")) {
        paste0("# remove random seeds\nrm(.Random.seed)\nset.seed(as.numeric(Sys.getenv('SLURM_JOB_ID')))")
        
      },
      "",
      "# load shared object if precompiled",
      load_so,
      "",
      "files <- list.files(pattern = '.so$')",
      "for (f in files) dyn.load(f)",
      "",
      "# List of variablevalues",
      var_list,
      "",
      "# Define variable values per run",
      var_per_run,
      "",
      "# Fixed parameters",
      fixedpars,
      "",
      "",
      "",
      "# Paste function call",
      paste0("cluster_result <- try(", as.character(expr),")"),
      sep = "\n",
      "save(cluster_result, file = paste0(jobname,'_', node_ID, '_result.RData'))" #'_',job_ID, 
    ),
    file = paste0(wd_path, jobname,".R")
  )
  
  
  
  
  # WRITE BASH
  cat(
    paste(
      "#!/bin/bash",
      "",
      "# Job name",
      paste0("#SBATCH --job-name=",jobname),
      "# Define format of output, deactivated",
      paste0("#SBATCH --output=",jobname,"_%j-%a.out"),
      "# Define format of errorfile, deactivated",
      paste0("#SBATCH --error=",jobname,"_%j-%a.err"),
      "# Define partition",
      paste0("#SBATCH --partition=", partition),
      "# Define number of nodes per task",
      paste0("#SBATCH --nodes=", nodes),
      "# Define number of cores per node",
      paste0("#SBATCH --ntasks-per-node=",cores),
      "# Define walltime",
      paste0("#SBATCH --time=",walltime),
      "# Define of repetition",
      paste0("#SBATCH -a 0-", num_nodes),
      "",
      "",
      "# Load R modules",
      "module load math/R",
      # paste0("export OPENBLAS_NUM_THREADS=",cores),
      paste0("export OMP_NUM_THREADS=","1"), # paste0("export OMP_NUM_THREADS=",cores),
      paste0("export MKL_NUM_THREADS=", "1"), # paste0("export MKL_NUM_THREADS=",cores),
      "",
      "# Run R script",
      paste0("Rscript ", jobname, ".R"),
      sep = "\n" 
    ),
    file = paste0(wd_path,jobname,".sh")
  )
  
  
  if (compile) {
    compile_files <- Sys.glob(paste0("*.cpp"))
    compile_files <- append(compile_files, Sys.glob(paste0("*.c")))
    compile_files <- paste(compile_files, collapse = " ")
    
    tar_locale <- paste0(
      "tar -jcf - ", compile_files, " ",wd_path, "*"
    )
    tar_remote <- paste0(
      "tar -C ./ -jxf - ; mv -t ./", jobname, "_folder ",compile_files,"; "
    )
    
    # list of files to compile
    sourcefiles <- paste(
      paste0( 
        # jobname, "_folder/",
        c(list.files(pattern = glob2rx("*.c")), list.files(pattern = glob2rx("*.cpp")))
      ), 
      collapse = " "
    )
    
    compile_remote <- paste0(
      # " module load math/R; R CMD SHLIB  -o ", jobname, "_shared_object.so ", sourcefiles," ; "
      " module load math/R; R CMD SHLIB ", sourcefiles, " -o ", jobname, "_shared_object.so; "
    )
    # compile_remote <- paste0(
    #   " module load math/R; R CMD cat(getwd()) "
    # )
  } else {
    compile_files <- Sys.glob(paste0("*.so"))
    compile_files <- append(compile_files, Sys.glob(paste0("*.o")))
    compile_files <- paste(compile_files, collapse = " ")
    
    tar_locale <- paste0(
      "tar -jcf - ", compile_files, " ",wd_path, "*"
    )
    tar_remote <- paste0(
      "tar -C ./ -jxf - ; mv -t ./", jobname, "_folder ",compile_files,"; "
    )
    
    # list of files to compile
    sourcefiles <- paste(
      paste0( 
        # jobname, "_folder/",
        c(list.files(pattern = glob2rx("*.c")), list.files(pattern = glob2rx("*.cpp")))
      ), 
      collapse = " "
    )
    compile_remote <- ""
  }
  ##
  # transfer and run files
  system(
    paste0(
      tar_locale, # compress all files in the local working dir
      " | ", ssh_command, machine, # pipe to ssh session on remote
      " 'if [ -d ", jobname, "_folder ]; then rm -Rf ", jobname,"_folder; fi ;", # remove folder if it exists
      " mkdir -p ", jobname,"_folder; ", # create new wd on remote
      tar_remote, # uncompress files in wd on remote, if necessary move files
      "cd ", jobname, "_folder; ", # change in said wd
      compile_remote, # compile files if said so, if not nothing happen
      "sbatch ", jobname, ".sh'" # start bash script
    )
  )
  
  setwd(original_wd)
  return(out)
}



#' Generate parameter list for distributed profile calculation
#' 
#' @description Generates list of \code{WhichPar} entries to facillitate distribute
#' profile calculation.
#' @details Lists to split the parameters for which the profiles are calculated
#' on the different nodes.
#' 
#' @param parameters list of parameters 
#' @param fits_per_node numerical, number of parameters that will be send to each node.
#' 
#' @return List with two arrays: \code{from} contains the number of the starting
#' parameter, while \code{to} stores the respective upper end of the parameter list
#' per node.
#' @examples
#' \dontrun{
#' parameter_list <- setNames(1:10, letters[1:10])
#' var_list <- profile_pars_per_node(parameter_list, 4)
#' }
#' 
#' @export
profile_pars_per_node <- function(parameters, fits_per_node) {
  # get the number of parameters
  n_pars <- length(parameters)
  
  # Get number of fits per node
  fits_per_node <- fits_per_node
  
  # determine the number of nodes necessary
  no_nodes <- 1:ceiling(n_pars/fits_per_node)
  
  # generate the lists which parameters are send to wich node
  pars_from <- fits_per_node
  pars_to_vec <- fits_per_node
  while (pars_from < (n_pars)) {
    pars_from <- pars_from + fits_per_node
    pars_to_vec <- c(pars_to_vec, pars_from)
  }
  pars_to_vec[length(pars_to_vec)] <- n_pars
  
  pars_from_vec <- c(1, pars_to_vec+1)
  pars_from_vec <- head(pars_from_vec, -1)
  
  out <- list(from=pars_from_vec, to=pars_to_vec)
  
  return(out)
}



## Use Julia to calculate steady states -----------------------------------------

#' Install the julia setup  
#' 
#' @description Installs Julia and the necessary julia packages
#' 
#' 
#' @param installJulia boolean, default \code{false}. If set to true, juliaup and via this then Julia is installed. 
#' @param installJuliaPackages boolean, default \code{true}. If set to true, the necessary packages are installed.
#' 
#' @return nothing
#' 
#' @export
installJuliaForSteadyStates <- function(installJulia = FALSE, installJuliaPackages = TRUE) {
  
  tryCatch(
    {
      system("git clone git@github.com:SeverinBang/JuliaSteadyStates.git ~/.JuliaSteadyStates/")
    },
    finally = {
      cat("github.com:SeverinBang/JuliaSteadyStates.git could not be cloned (again), check if ~/.JuliaSteadyStates already exists, if not write Severin your github username to be added to the repository")
    }
  )
  
  # install Julia
  if (installJulia) {
    system("sh -i ~/.JuliaSteadyStates/installJuliaUp.sh -y")
    system("juliaup add release")
  }
  
  if (installJuliaPackages) {
    # install packages
    system("julia -e 'using Pkg; Pkg.add(\"CSV\")'")
    system("julia -e 'using Pkg; Pkg.add(\"DataFrames\")'")
    system("julia -e 'using Pkg; Pkg.add(\"Symbolics\")'")
    system("julia -e 'using Pkg; Pkg.add(\"Catalyst\")'")
    system("julia -e 'using Pkg; Pkg.add(\"SymbolicUtils\")'")
    system("julia -e 'using Pkg; Pkg.add(\"Graphs\")'")
  }
  
}


#' Calculated the steady states of a given model  
#' 
#' @description Uses julia to calculate the steady state transformations
#' 
#' @param el the equation list
#' @param forcings vector of strings, default \code{c("","")}. The names of the forcings which will be set to zero.
#' @param neglect vector of strings, default \code{c("","")}. The names of the variables which will be neglected as fluxParameters and therefore will not be solved for.
#' @param verboseLevel integer, default \code{1}. The level of verbosity of the output, right now only 1 (all) and 0 (no) is implementes.
#' 
#' @return named vector with the steady state transformations. The names are the inner, the values are the outer parameters
#' 
#' @export
steadyStateToolJulia <- function(
    el,
    forcings = NULL,
    neglect = NULL,
    verboseLevel = 1
) {
  # prepare things:
  myWD <- getwd()
  dModEqnFileName = "EquationsForSteadyStates"
  
  dMod::write.eqnlist(el, file = paste0(dModEqnFileName, ".csv"))
  
  inputPath <- file.path(myWD, paste0(dModEqnFileName, ".csv"))
  fileName <- "SteadyStatesFromJulia"
  
  if (is.null(forcings)) {
    forcings <- c("","")
  }
  
  if (is.null(neglect)) {
    neglect <- c("","")
  }
  
  # load julia
  if (!requireNamespace("JuliaCall", quietly = TRUE)) {
    warning("The 'JuliaCall' package must be installed.")
    return(NULL)
  }
  if (dir.exists(file.path(Sys.getenv("HOME"),".juliaup/bin"))) {
    JuliaCall::julia_setup(JULIA_HOME = file.path(Sys.getenv("HOME"),".juliaup/bin"))
  } else if (file.exists("/usr/bin/julia")) {
    JuliaCall::julia_setup(JULIA_HOME = "/usr/bin/")
  } else {
    stop("No Julia installation found, please use juliaup to install julia.")
    return(NULL)
  }
  JuliaCall::julia_setup(JULIA_HOME = file.path(Sys.getenv("HOME"),".juliaup/bin"))
  
  
  # call the julia steady state tool:
  julia_source(file.path(Sys.getenv("HOME"),".JuliaSteadyStates/ODESteadyStateTrafo_function.jl"))
  
  julia_call("determineSteadyStateTrafos", inputPath, forcings, neglect, myWD, fileName, verboseLevel = julia_eval(paste0("Int(", verboseLevel, ")")))
  
  # load the results
  steadyStatesFile = read.csv(paste0(myWD,"/",fileName, ".csv" ), dec = ".", sep = ",")
  
  steadyStates = data.table(keys = steadyStatesFile$Keys, values = steadyStatesFile$Values)
  steadyStates = steadyStates[!(keys %in% forcings)]
  
  
  sstates <- steadyStates$values
  
  sstates <- str_replace_all(sstates, "_init", "")
  
  names(sstates) <- steadyStates$keys
  
  
  return(sstates)
}




