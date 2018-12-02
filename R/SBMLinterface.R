#' Import an SMBL model
#'
#' Requires AMICI https://github.com/ICB-DCM/AMICI/ and dependencies to be installed on your system
#' Big thanks go to Daniel Weindl!
#' 
#' @param modelpath Path to the sbml file
#' @param amicipath Path to your amici-python-installation, e.g.: AMICIPATH/python
#'
#' @return list of eqnlist, parameters and inits
#' @export
#' @importFrom rjson fromJSON 
#' @importFrom stringr str_replace_all
import_sbml <- function(modelpath, amicipath = NULL) {
  
  importscript <- system.file("code/sbmlAmiciDmod.py", package = "dMod")
  tmpfile_json <- tempfile()
  
  if (!is.null(amicipath))
    amicipath <- paste0("PYTHONPATH=", amicipath)
  
  run_import_script_call <- paste0('bash -c "', "cd ~/.virtualenvs/amici/bin && source activate &&", amicipath, "&&",  " python ", importscript, " ", modelpath, " ", tmpfile_json, '"')
  system(run_import_script_call)
  json_content <- rjson::fromJSON(file = tmpfile_json)
  
  S <- do.call(cbind, json_content[["S"]])
  S[S==0] <- NA
  
  v <- json_content[["v"]]
  v <- stringr::str_replace_all(v, "\\*\\*", "^")
  
  states <- json_content[["stateNames"]]
  
  pars <- setNames(json_content[["p"]], json_content[["parameterNames"]])
  x0 <- setNames(json_content[["x0"]], json_content[["stateNames"]])
  
  reactions <- eqnlist(smatrix = S, states = states, rates = v)
  
  observables <- json_content[["observables"]]
  
  print(names(json_content))
  
  out <- list(reactions = reactions, pars = pars, inits = x0, observables = observables)
  return(out)
}
