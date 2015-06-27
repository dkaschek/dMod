#' Select parameter values with lowest Chi^2 among profiles.
#' 
#' @param prf A profiles as returned from \link[dMod]{profile}.
#' @param context If TRUE, the chi^2 and other context of the profile is
#'                returned. This parameter defaults to FALSE in which case the
#'                output can be used as outer parameters directly.
#'                
#' @return Parameter values with lowest chi^2 w/ or w/o profile context. If all
#'         profiles are invalid, NULL is returned.
#' 
#' @details 
#' On profiling a model, parameter values yielding a lower chi^2 than the one of
#' the fit providing the optimal values for the profiles might be encountered.
#' This function extracts the set of parameter values possessing the lowest
#' chi^2 among all profiles. Profiles which are invalid due to e.g., integration
#' problem on calculating them, are ignored. If all profiles are invalid, NULL
#' is returned.
#' 
#' @export

plSelectMin <- function(prf, context = FALSE) {
  
  # Remove invalid profiles.
  apt <- sapply(prf, class) == "matrix"
  if (!any(apt)) {
    return(NULL)
  } else {
    prf <- prf[apt]
  }
  
  # Minium chi^2 parameter values sets per profile.
  chi2MinAllProfiles <- sapply(prf, function(species) {
    return(species[which.min(species[, "value"]), ])        
  })
  
  # Minimum chi^2 parameter values set accross all profiles.
  chi2MinBest <- chi2MinAllProfiles[, which.min(chi2MinAllProfiles[1, ])]
  
  if (context) {
    return(chi2MinBest)   
  } else {
    return(chi2MinBest[-(1:4)])  
  }
  
}
