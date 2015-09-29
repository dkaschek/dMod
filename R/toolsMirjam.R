#' Print without attributes
#' 
#' @param x The object to be printed out
#' @param list_attributes prints a list of attribute names if TRUE (default=TRUE)
#' @details To suppress the printout of attributes like "deriv". 
#' @export
print0 <- function(x, list_attributes=TRUE ) {
  attributes_all <- names(attributes(x))
  attributes_rm <- attributes_all[!(attributes_all %in% c("dim","names","dimnames","row.names","col.names"))]
  attributes(x)[attributes_rm] <- NULL
  print.default(x)
  if(list_attributes)
    cat("Attributes:",attributes_all)
}

#' Find common elements in lists
#' @description generalisation of \link{intersect} to a list of vectors.
#' @param list contains a set of vectors.
#' @param byNames if set TRUE common names are checked instead of common values (default: TRUE) 
#' @examples testList <-list(c(a=1,b=5,c=3,e=5), c(d=3,b=1,q=1,c=5,i=2)) 
#' intersectList(testList,FALSE) 
#' intersectList(testList,TRUE)
intersectList <- function(list, byNames=TRUE){
  inter <- list[[1]]
  if(byNames)
    inter <- names(inter)
  lapply(list[-1], function(l){
    if(byNames)
      l <- names(l)
    inter <<- intersect(inter,l)
  })
  return(inter)
}

#' Find fluxes with low contributions
#' 
#' @param out list of integrated fluxes, named "state.fluxEquation". States are not allowed to contain a "."-symbol.
#' @details  This can be used if the fluxes were integrated (by using \link{Xf} to integrate the fluxes),
#' otherwise see \link{getZeroFluxes}.
getZeroFluxesInt <- function(out, rtol = .05, atol = 0) {
  

  
  states <- sapply(strsplit(names(out)[-1], ".", fixed=TRUE), function(v) v[1])
  fluxes <- sapply(strsplit(names(out)[-1], ".", fixed=TRUE), function(v) paste(v[-1], collapse=""))
  unique.states <- unique(states)
  unique.fluxes <- unique(fluxes)
  
  out.groups <- lapply(unique.states, function(s) {
    
    selected <- which(states == s)
    out.selected <- out[selected+1]
    names(out.selected) <- fluxes[selected]
    
    abssum.extended <- rep(0, length(unique.fluxes))
    names(abssum.extended) <- unique.fluxes    
    abssum.extended[names(out.selected)] <- out.selected
    
    # Normalize with respect to the L1 norm of the state derivative (sum of all fluxes)
    state.dot <- sum(out.selected)
    if(state.dot < 1e-10){
      state.dot <- 1e-10
    }
    abssum.normed.extended <- abssum.extended/state.dot
    
    return(list(abssum.extended, abssum.normed.extended))
    
  })
  #rownames(out.groups) <- unique.states
  out.groups.abs <- do.call(rbind, lapply(out.groups, function(g) g[[1]]))
  rownames(out.groups.abs) <- unique.states
  out.groups.rel <- do.call(rbind, lapply(out.groups, function(g) g[[2]]))
  rownames(out.groups.rel) <- unique.states
  
  zero.fluxes.abs <- unique.fluxes[apply(out.groups.abs, 2, function(v) all(v < atol))]
  zero.fluxes.rel <- unique.fluxes[apply(out.groups.rel, 2, function(v) all(v < rtol))]
  zero.fluxes <- union(zero.fluxes.abs, zero.fluxes.rel)
  non.zero.fluxes <- unique.fluxes[!unique.fluxes%in%zero.fluxes]
  
  #return(list(out.groups, zero.fluxes, non.zero.fluxes))
  return(list(fluxes.abs = out.groups.abs, 
              fluxes.rel = out.groups.rel, 
              fluxes.zero = zero.fluxes, 
              fluxes.nonzero = non.zero.fluxes))
  
}
