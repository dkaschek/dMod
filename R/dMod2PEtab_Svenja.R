# -------------------------------------------------------------------------#
# PEtab export stuff ----
# -------------------------------------------------------------------------#
#
# [PURPOSE]
#
#
# [AUTHOR]
# Svenja Kemmer
#
# [Date]
# 2021-03-02 
#
rm(list = ls(all.names = TRUE))
.currentwd <- getwd()
.exportFolder <- "PEtabExports"

# -------------------------------------------------------------------------#
# General definitions ----
# -------------------------------------------------------------------------#

modelname <- "TestCase"

#' Export observables to PEtab. 
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


