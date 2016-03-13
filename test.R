Psimu <- function(parameters = NULL, condition = NULL) {
  

  p2p <- function(pars, fixed=NULL, deriv = TRUE) {
    
    prediction <- prediction[[condition]]
    pars <- as.parvec(pars)
    
    # Get parameter vector
    overlap <- intersect(colnames(prediction), names(pars))
    if (length(overlap) == 0) return(pars)
    
    last.row <- as.numeric(prediction[nrow(prediction), overlap])
    
    out <- unclass(pars)
    out[overlap] <- last.row
    
    jac.out <- attr(pars, "deriv")
    
    print(out)
    
    # Get derivatives
    myderiv <- NULL
    dstates <- attr(prediction, "deriv")
    
    if (deriv && !is.null(dstates)) {
      
      states <- overlap
      pars.names <- names(pars)
      
      names.last.row.deriv <- outer(states, pars.names, function(x, y) paste(x, y, sep = "."))
      last.row.deriv <- structure(rep(0, length(names.last.row.deriv)),
                                  names = names.last.row.deriv)
      
      print(last.row.deriv)
      
      overlap <- intersect(names.last.row.deriv, colnames(dstates))
      
      last.row.deriv[overlap] <- as.numeric(dstates[nrow(dstates), overlap])
      
      
      jac.matrix <- t(matrix(last.row.deriv, nrow = length(pars), ncol = length(states)))
      rownames(jac.matrix) <- states
      colnames(jac.matrix) <- pars.names
      
      jac.matrix <- jac.matrix %*% submatrix(jac.out, rows = colnames(jac.matrix))
      jac.out[rownames(jac.matrix), ] <- jac.matrix
      
      print(jac.matrix)
      
    }
    
    as.parvec(out, deriv = jac.out)
    
    
  }
  
  parfn(p2p, parameters, condition)
  
}