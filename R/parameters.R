#' Parameter transformation
#' 
#' @param P named character vector. Names correspond to the parameter names to be fed into
#' the model. The values of P are equations that express the parameters in terms of other
#' parameters
#' @param parameters character vector. Optional. If given, the generated parameter
#' transformation returns values for each element in \code{parameters}. If elements of
#' \code{parameters} are not in \code{names(P)} the identity transformation is assumed.
#' @return a function \code{p(P)} representing the parameter transformation. The result of p(P) 
#' contains an attribute "deriv" which is the jacobian of the parameter transformation
#' evaluated at \code{P}.
# P <- function(P, parameters=NULL) {
#   
#   # get outer parameters
#   fparse <- getParseData(parse(text=P))
#   symbols <- unique(fparse$text[fparse$token == "SYMBOL"])
#   if(is.null(parameters)) {
#     parameters <- symbols 
#   } else {
#     identity <- parameters[which(!parameters%in%symbols)]
#     names(identity) <- identity
#     P <- c(P, identity)
#   }
#   
#   # expresion list for parameter and jacobian evaluation
#   expressionList <- lapply(P, function(myrel) parse(text=as.character(myrel)))
#   listJac <- unlist(lapply(parameters, function(var) {
#     unlist(lapply(expressionList, function(myexp) paste(deparse(D(myexp, as.character(var))), collapse="")))
#   }))
#   
#   jacNames <- expand.grid.alt(names(P), parameters)
#   jacNames <- paste(jacNames[,1], jacNames[,2], sep=".")
#   
#   dP <- listJac; names(dP) <- jacNames
# 
#   PEval <- funC.algebraic(P)
#   dPEval <- funC.algebraic(dP)
#   
#   # the parameter transformation function to be returned
#   p2p <- function(p, fixed=NULL, derivs = TRUE) {
#     
#     # replace fixed parameters by values given in fixed
#     if(!is.null(fixed)) {
#       is.fixed <- which(names(p)%in%names(fixed)) 
#       if(length(is.fixed)>0) p <- p[-is.fixed]
#       p <- c(p, fixed)
#     }
#     
#     # check for parameters which are not defined in parameters
#     emptypars <- names(p)[!names(p)%in%parameters]
#     
#     x <- as.list(p)
#     values <- PEval(x)[1,]
#     jacValues <- dPEval(x)
#       
#     
#     jacobian <- matrix(jacValues, ncol=length(parameters), nrow=length(P))
#     if(!is.null(fixed)) {
#       jacobian <- matrix(jacobian[,-which(parameters %in% names(fixed))], nrow=length(P))
#       colnames(jacobian) <- parameters[!parameters%in%names(fixed)]
#       rownames(jacobian) <- names(P)
#     } else {
#       colnames(jacobian) <- parameters
#       rownames(jacobian) <- names(P)
#     }
#     
#     # Append zeros for emptypars
#     if(length(emptypars) > 0) {
#       empty <- matrix(0, nrow=dim(jacobian)[1], ncol=length(emptypars))
#       colnames(empty) <- emptypars
#       jacobian <- cbind(jacobian, empty)
#     }    
#     if(derivs) attr(values, "deriv") <- jacobian
#     
#     
#     
#     return(values)
#   }
#   
#   
#   return(p2p)
#   
# }


#' Parameter transformation
#' 
#' @param P named character vector. Names correspond to the parameter names to be fed into
#' the model. The values of P are equations that express the parameters in terms of other
#' parameters
#' @param parameters character vector. Optional. If given, the generated parameter
#' transformation returns values for each element in \code{parameters}. If elements of
#' \code{parameters} are not in \code{names(P)} the identity transformation is assumed.
#' @return a function \code{p(P)} representing the parameter transformation. The result of p(P) 
#' contains an attribute "deriv" which is the jacobian of the parameter transformation
#' evaluated at \code{P}.
P <- function(P, parameters=NULL) {
  
  # get outer parameters
  fparse <- getParseData(parse(text=P))
  symbols <- unique(fparse$text[fparse$token == "SYMBOL"])
  if(is.null(parameters)) {
    parameters <- symbols 
  } else {
    identity <- parameters[which(!parameters%in%symbols)]
    names(identity) <- identity
    P <- c(P, identity)
  }
  
  # expresion list for parameter and jacobian evaluation
  expressionList <- lapply(P, function(myrel) parse(text=as.character(myrel)))
  listExpression <- parse(text = paste("list(", paste(P, collapse=", "), ")"))
  listJac <- unlist(lapply(parameters, function(var) {
    unlist(lapply(expressionList, function(myexp) paste(deparse(D(myexp, as.character(var))), collapse="")))
  }))
  listJac <- parse(text = paste("list(", paste(listJac, collapse=", "), ")"))
  jacNames <- expand.grid(names(P), parameters)
  
  # the parameter transformation function to be returned
  p2p <- function(p, fixed=NULL, derivs = TRUE) {
    
    # replace fixed parameters by values given in fixed
    if(!is.null(fixed)) {
      is.fixed <- which(names(p)%in%names(fixed)) 
      if(length(is.fixed)>0) p <- p[-is.fixed]
      p <- c(p, fixed)
    }
    
    # check for parameters which are not defined in parameters
    emptypars <- names(p)[!names(p)%in%parameters]
    
    # compute transformation output
    out <- with(as.list(p), { 
      
      values <- unlist(eval(listExpression))
      names(values) <- names(P)
      
      jacValues <- unlist(eval(listJac))
      jacobian <- matrix(jacValues, ncol=length(parameters), nrow=length(P))
      if(!is.null(fixed)) {
        jacobian <- matrix(jacobian[,-which(parameters %in% names(fixed))], nrow=length(P))
        colnames(jacobian) <- parameters[!parameters%in%names(fixed)]
        rownames(jacobian) <- names(P)
      } else {
        colnames(jacobian) <- parameters
        rownames(jacobian) <- names(P)
      }
      
      # Append zeros for emptypars
      if(length(emptypars) > 0) {
        empty <- matrix(0, nrow=dim(jacobian)[1], ncol=length(emptypars))
        colnames(empty) <- emptypars
        jacobian <- cbind(jacobian, empty)
      }    
      if(derivs) attr(values, "deriv") <- jacobian
      
      return(values)
      
    })
    
    
    return(out)
  }
  
  
  return(p2p)
  
}

