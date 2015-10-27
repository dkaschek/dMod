
## Class "eqnlist" and its constructor ------------------------------------------


#' Coerce to an equation list
#' @description Translates a reaction network, e.g. defined by a data.frame, into an equation list object.
#' @param data an R object, usually a \code{data.frame}.
#' @param volumes named character, volume parameters for states. Names must be a subset of the states.
#' Values can be either characters, e.g. "V1", or numeric values for the volume. If \code{volumes} is not
#' \code{NULL}, missing entries are treated as 1.
#' @param ... additional arguments to be passed to or from methods.
#' @details If \code{data} is a \code{data.frame}, it must contain columns "Description" (character), 
#' "Rate" (character), and one column per ODE state with the state names. 
#' The state columns correspond to the stoichiometric matrix.
#' @return Object of class \link{eqnlist}
#' @export
as.eqnlist <- function(data, volumes, ...) {
  UseMethod("as.eqnlist", data)
}

#' @export
as.eqnlist.data.frame <- function(data, volumes = NULL) {
  description <- as.character(data$Description)
  rates <- as.character(data$Rate)
  states <- setdiff(colnames(data), c("Description", "Rate"))
  smatrix <- as.matrix(data[, states]); colnames(smatrix) <- states
  
  eqnlist(smatrix, states, rates, volumes, description)
  
}


#' @export
is.eqnlist <- function(x) {
  
  #Empty list
  if (is.null(x$smatrix)) {
    if (length(x$states) == 0 &&
        length(x$rates) == 0 &&
        is.null(x$volumes) &&
        length(x$description) == 0
    ) {
      return(TRUE)
    } else {
      return(FALSE)
    }
  } else {
    #Non-empty list
    if (inherits(x, "eqnlist") &&
        all(names(x) == c("smatrix", "states", "rates", "volumes", "description")) &&
        all(names(x$smatrix) == names(x$states)) &&
        dim(x$smatrix)[1] == length(x$rates) &&
        dim(x$smatrix)[2] == length(x$states) &&
        is.matrix(x$smatrix)
    ) {
      return(TRUE)
    } else {
      return(FALSE)
    }
  }
}


## Class "eqnlist" and its methods ------------------------------------------

#' Generate a table of reactions (data.frame) from an equation list
#' 
#' @param eqnlist object of class \link{eqnlist}
#' @return \code{data.frame} with educts, products, rate and description. The first
#' column is a check if the reactions comply with reaction kinetics.
#' @export
getReactions <- function(eqnlist) {
  
  # Extract information from eqnlist
  S <- eqnlist$smatrix
  rates <- eqnlist$rates
  description <- eqnlist$description
  variables <- eqnlist$states
  
  # Determine lhs and rhs of reactions
  if(is.null(S)) return()
  
  reactions <- apply(S, 1, function(v) {
    
    numbers <- v[which(!is.na(v))]
    educts <- -numbers[numbers < 0]
    educts[which(educts == 1)] <- " "
    products <- numbers[numbers > 0]
    products[which(products == 1)] <- " "
    educts <- paste(paste(educts, names(educts), paste = ""), collapse="+")
    products <- paste(paste(products, names(products), paste = ""), collapse="+")
    
    reaction <- paste(educts, "->", products)
    return(c(educts, products))
    
  })
  educts <- reactions[1,]
  products <- reactions[2,]
  
  # Determine number of conserved quantities
  S[is.na(S)] <- 0
  v <- MASS::Null(t(S))
  cq <- ncol(v)
  
  # Check for consistency
  exclMarks.logical <- unlist(lapply(1:length(rates), function(i) {
    
    myrate <- rates[i]
    parsedRate <- getParseData(parse(text=myrate, keep.source = TRUE))
    symbols <- parsedRate$text[parsedRate$token=="SYMBOL"]
    
    educts <- variables[which(S[i,]<0)]
    
    !all(unlist(lapply(educts, function(e) any(e==symbols))))
    
  }))
  exclMarks <- rep(" ", ncol(reactions))
  exclMarks[exclMarks.logical] <- "!"
  

  # Generate data.frame  
  out <- data.frame(exclMarks, educts, "->", products, rates, description, stringsAsFactors = FALSE)
  colnames(out) <- c("Check", "Educt",  "->",  "Product", "Rate", "Description")
  rownames(out) <- 1:nrow(out)
  
  # Attributes
  attr(out, "cq") <- cq
  
  return(out)
  
}

#' @export
addReaction <- function(x, ...) {
  UseMethod("addReaction", x)
}
#' Add reaction to reaction table
#' 
#' @param from character with the left hand side of the reaction, e.g. "2*A + B"
#' @param to character with the right hand side of the reaction, e.g. "C + 2*D"
#' @param rate named character. The rate associated with the reaction. The name is employed as a description
#' of the reaction.
#' @param f equation list, see \link{eqnlist}
#' @return An object of class \link{eqnlist}.
#' @examples 
#' \dontrun{
#' f <- addReaction("2*A+B", "C + 2*D", "k1*B*A^2", NULL)
#' f <- addReaction("C + A", "B + A", "k2*C*A", f)
#' }
#' @export
addReaction.eqnlist <- function(eqnlist, from, to, rate) {
  
  volumes <- eqnlist$volumes
  
  # Analyze the reaction character expressions
  educts <- getSymbols(from)
  eductCoef <- 0
  if(length(educts) > 0) eductCoef <- sapply(educts, function(e) sum(getCoefficients(from, e)))
  products <- getSymbols(to)
  productCoef <- 0
  if(length(products) > 0) productCoef <- sapply(products, function(p) sum(getCoefficients(to, p)))
  
  
  # States
  states <- unique(c(educts, products))
  
  # Description
  description <- names(rate)
  if(is.null(description)) description <- ""
  
  # Stoichiometric matrix
  smatrix <- matrix(NA, nrow = 1, ncol=length(states)); colnames(smatrix) <- states
  if(length(educts)>0) smatrix[,educts] <- -eductCoef
  if(length(products)>0) {
    filled <- !is.na(smatrix[,products])
    smatrix[,products[filled]] <- smatrix[,products[filled]] + productCoef[filled]
    smatrix[,products[!filled]] <- productCoef[!filled]  
  }
  
  
  smatrix[smatrix == "0"] <- NA
  
  
  # data.frame
  mydata <- cbind(data.frame(Description = description, Rate = as.character(rate)), as.data.frame(smatrix))
  row.names(mydata) <- NULL
  
  
  if(!is.null(eqnlist)) {
    mydata0 <- as.data.frame(eqnlist)
    mydata <- combine(mydata0, mydata)
  }
  
  
  as.eqnlist(mydata, volumes = volumes)
  
}


#' Generate list of fluxes from equation list
#' 
#' @param eqnlist object of class \link{eqnlist}.
#' @return list of named characters, the in- and out-fluxes for each state.
#' @export
getFluxes <- function(eqnlist) {
  
  description <- eqnlist$description
  rate <- eqnlist$rates
  variables <- eqnlist$states
  SMatrix <- eqnlist$smatrix
  volumes <- eqnlist$volumes
  
  if(is.null(SMatrix)) return()
  
  volumes.draft <- structure(rep("1", length(variables)), names = variables)
  volumes.draft[names(volumes)] <- volumes
  volumes <- volumes.draft
  
  
  
  
  # generate equation expressions
  terme <- lapply(1:length(variables), function(j) {
    v <- SMatrix[,j]
    nonZeros <- which(!is.na(v))
    var.description <- description[nonZeros]
    positives <- which(v > 0)
    negatives <- which(v < 0)
    volumes.destin <- volumes.origin <- rep(volumes[j], length(v))
    if(length(positives) > 0) {
      volumes.origin[positives] <- sapply(positives, function(i) {
        candidates <- which(SMatrix[i,] < 0)
        myvolume <- unique(volumes[candidates])
        if(length(myvolume) > 1) 
          stop("species from different compartments meet in one reaction")
        if(length(myvolume) == 0) myvolume <- volumes[j]
        
        return(myvolume)
      })
    }
    volumes.ratios <- paste0("*(", volumes.origin, "/", volumes.destin, ")")
    #print(volumes.ratios)
    volumes.ratios[volumes.destin == volumes.origin] <- ""
    
    numberchar <- as.character(v)
    if(nonZeros[1] %in% positives){
      numberchar[positives] <- paste(c("", rep("+", length(positives)-1)), numberchar[positives], sep = "") 
    } else {
      numberchar[positives] <- paste("+", numberchar[positives], sep = "")
    }
    var.flux <- paste0(numberchar[nonZeros], "*(",  rate[nonZeros], ")", volumes.ratios[nonZeros])
    names(var.flux) <- var.description
    return(var.flux)
  })
  
  fluxes <- terme
  names(fluxes) <- variables
  
  return(fluxes)
  
  
}

#' Coerce equation list into a data frame
#' 
#' @param eqnlist object of class \link{eqnlis}
#' @return a \code{data.frame} with columns "Description" (character), 
#' "Rate" (character), and one column per ODE state with the state names. 
#' The state columns correspond to the stoichiometric matrix.
#' @export
as.data.frame.eqnlist <- function(eqnlist) {
  
  if(is.null(eqnlist$smatrix)) return()
  
  data <- data.frame(Description = eqnlist$description,
                     Rate = eqnlist$rate,
                     eqnlist$smatrix, 
                     stringsAsFactors = FALSE)
  
  attr(data, "volumes") <- eqnlist$volumes
  
  return(data)
}

#' Write equation list into a csv file
#' 
#' @param eqnlist object of class \link{eqnlist}
#' @param ... Arguments going to \link[utils]{write.csv}
#' 
#' @export
write.eqnlist <- function(eqnlist, ...) {
  
  data <- as.data.frame(eqnlist)
  write.csv(data, ...)
  
}


#' subset of an equation list
#' 
#' @param x the equation list
#' @param ... logical expression for subsetting
#' @details The argument \code{...} can contain "Educt", "Product", "Rate" and "Description".
#' The "%in%" operator is modified to allow searches in Educt and Product (see examples).
#' 
#' @return An object of class \link{eqnlist}
#' @examples
#' reactions <- data.frame(Description = c("Activation", "Deactivation"), 
#'                         Rate = c("act*A", "deact*pA"), A=c(-1,1), pA=c(1, -1) )
#' f <- as.eqnlist(reactions)
#' subset(f, "A" %in% Educt)
#' subset(f, "pA" %in% Product)
#' subset(f, grepl("act", Rate))
#' @export subset.eqnlist
#' @export
subset.eqnlist <- function(eqnlist, ...) {
  
  # Do selection on data.frame
  data <- getReactions(eqnlist)
  if(is.null(data)) return()
  
  data.list <- list(Educt = lapply(data$Educt, getSymbols), 
                    Product = lapply(data$Product, getSymbols),
                    Rate = data$Rate,
                    Description = data$Description,
                    Check = data$Check)
  
  "%in%" <- function(x, table) sapply(table, function(mytable) any(x == mytable))
  select <- which(eval(substitute(...), data.list))
  
  # Translate subsetting on eqnlist entries
  # smatrix
  smatrix <- submatrix(eqnlist$smatrix, rows = select)
  empty <- sapply(1:ncol(smatrix), function(i) all(is.na(smatrix[, i])))
  smatrix <- submatrix(smatrix, cols = !empty)
  
  # states and rates
  states <- colnames(smatrix)
  rates <- eqnlist$rates[select]
  
  # volumes
  volumes <- eqnlist$volumes
  if(!is.null(volumes)) volumes <- volumes[names(volumes) %in% states]
  
  # description
  description <- eqnlist$description[select]
  
  eqnlist(smatrix, states, rates, volumes, description)
  
  
}


#' Print equation list
#' 
#' @param eqnlist object of class \link{eqnlist}
#' @export
print.eqnlist <- function(eqnlist) {
  
  cat("Reaction table:\n")
  reactions <- getReactions(eqnlist)
  print(reactions)
  cat("Number of conserved quantities:", attr(reactions, "cq"), "\n")
  
}

#' Pander on equation list
#' 
#' @param eqnlist object of class \link{eqnlist}
#' @return pander object
#' @export
pander.eqnlist<- function(eqnlist) {
  
  pander::panderOptions("table.alignment.default", "left")
  pander::panderOptions("table.split.table", Inf)
  pander::panderOptions("table.split.cells", Inf)
  out <- getReactions(eqnlist)
  exclude <- "Check"
  out <- out[, setdiff(colnames(out), exclude)]
  out$Rate <- paste0(format.eqnvec(as.character(out$Rate)))
  pander::pander(out)
  
}


## Class "eqnvec" and its constructors --------------------------------------------



#' Coerce to an equation vector
#' 
#' @param x an R object, usually a named character or an \link{eqnlist}.
#' @details If \code{x} is of class \code{eqnlist}, \link{getFluxes} is called and coerced
#' into a vector of equations.
#' @return object of class \link{eqnvec}.
#' @export
as.eqnvec <- function(x) {
  UseMethod("as.eqnvec", x)
}

#' @export
as.eqnvec.character <- function(char) {

  eqnvec(char, names(char))
  
}

#' Transform equation list into vector of equations
#' 
#' @description An equation list stores an ODE in a list format. The function
#' translates this list into the right-hand sides of the ODE.
#' @param eqnlist equation list, see \link{eqnlist}
#' @return An object of class \link{eqnvec}. 
#' @export
as.eqnvec.eqnlist <- function(eqnlist) {
  
  terme <- getFluxes(eqnlist)
  if(is.null(terme)) return()
  terme <- lapply(terme, function(t) paste(t, collapse=" "))
  
  
  terme <- do.call(c, terme)
  
  eqnvec(terme, names(terme))
  
}

#' @export
c.eqnlist <- function(...) {
  
  out <- lapply(list(...), as.data.frame)
  out <- Reduce(combine, out)
  
  as.eqnlist(out)
  
}


#' @export
is.eqnvec <- function(x) {
  if (inherits(x, "eqnvec") &&
      length(x) == length(names(x))
  )
    return(TRUE)
  
  else
    return(FALSE)
}


## Class "eqnvec" and its methods --------------------------------------------


#' Encode equation vector in format with sufficient spaces
#' 
#' @param eqnvec object of class \link{eqnvec}. Alternatively, a named parsable character vector.
#' @return named character
#' @export format.eqnvec
#' @export
format.eqnvec <- function(eqnvec) {
  
  eqns <- sapply(eqnvec, function(eqn) {
    parser.out <- getParseData(parse(text = eqn, keep.source = TRUE))
    parser.out <- subset(parser.out, terminal == TRUE)
    parser.out$text[parser.out$text == "*"] <- "·"
    out <- paste(parser.out$text, collapse = "")
    return(out)
  })
  
  patterns <- c("+", "-", "·", "/")
  for(p in patterns) eqns <- gsub(p, paste0(" ", p, " "), eqns, fixed = TRUE)
  
  return(eqns)
    
  
}

#' Print equation vector
#' 
#' @param object of class \link{eqnvec}.
#' @export
print.eqnvec <- function(eqnvec) {
  
  cat("Equations:\n")
  out <- as.data.frame(unclass(eqnvec))
  colnames(out) <- ""
  print(out)
}



#' Pander on equation vector
#'
#' @param object of class \link{eqnvec}. 
#' @return pander object
#' @export
pander.eqnvec <- function(eqnvec) {
  
  pander::panderOptions("table.alignment.default", "left")
  pander::panderOptions("table.split.table", Inf)
  pander::panderOptions("table.split.cells", Inf)
  out <- as.data.frame(unclass(eqnvec), stringsAsFactors = FALSE)
  colnames(out) <- "" #  as.character(substitute(eqnvec))
  out[, 1] <- format.eqnvec(out[, 1])
  pander::pander(out)
  
}


#' Symbolic time derivative of equation vector given an equation list
#' 
#' @param observable named character vector. Names correspond to observable names, the chars 
#' correspond to the observation function
#' @param eqnlist equation list
#' @details Observables are translated into an ODE
#' @return An object of class \link{eqnvec}
#' @examples 
#' reactions <- data.frame(Description = c("Activation", "Deactivation"), 
#'                         Rate = c("act*A", "deact*pA"), A=c(-1,1), pA=c(1, -1))
#' f <- generateEquations(reactions)
#' myobs <- c(tA = "s1*(pA + A)", dA = "s2*(pA-A)")
#' f <- addObservable(myobs, f)
#' @export
dot <- function(observable, eqnlist) {
  UseMethod("dot", observable)
}
#' @export
#' @rdname dot
dot.eqnvec <- function(observable, eqnlist) {
  
 
  # Analyze the observable character expression
  symbols <- getSymbols(observable)
  states <- intersect(symbols, eqnlist$states)
  derivatives <- lapply(observable, function(obs) {
    out <- lapply(as.list(states), function(x) paste(deparse(D(parse(text=obs), x), width.cutoff = 500),collapse=""))
    names(out) <- states
    return(out)
  })
  
  # Generate equations from eqnist
  f <- as.eqnvec(eqnlist)
  
  newodes <- sapply(derivatives, function(der) {
    
    out <- sapply(names(der), function(n) {
      d <- der[n]
      if(d != "0") paste( paste("(", d, ")", sep="") , paste("(", f[names(d)], ")",sep=""), sep="*") else return("0")
    })
    out <- paste(out, collapse = "+")
    
    return(out)
    
  })
  
  as.eqnvec(newodes)
}


#'@export
c.eqnvec <- function(...) {
 
  out <- lapply(list(...), unclass)
  out <- do.call(c, out)
  if (any(duplicated(names(out)))) {
    stop("Names must be unique")
  }
  
  as.eqnvec(out)
}

#' @export
"[.eqnvec" <- function(x, ...) {
  out <- unclass(x)[...]
  class(out) <- "eqnvec"
  return(out)
}


#' Evaluation of algebraic expressions defined by characters
#' 
#' @param x Name character vector, the algebraic expressions
#' @param compile Logical. The function is either translated into a C file to be compiled or is
#' evaluated in raw R.
#' @return A prediction function \code{f(mylist, attach.input = FALSE)} where \code{mylist} is a list of numeric 
#' vectors that can
#' be coerced into a matrix. The names correspond to the symbols used in the algebraic expressions. 
#' The argument \code{attach.input} determines whether \code{mylist} is attached to the output.
#' The function \code{f} returns a matrix.
#' @examples 
#' \dontrun{
#' myfun <- funC0(c(y = "a*x^4 + b*x^2 + c"))
#' out <- myfun(list(a = -1, b = 2, c = 3, x = seq(-2, 2, .1)), attach.input = TRUE)
#' plot(out[, "x"], out[, "y"])
#' }
#' 
#' @export
funC0 <- function(x, compile = FALSE, modelname = NULL) {
    
  # Get symbols to be substituted by x[] and y[]
  outnames <- names(x)
  innames <- getSymbols(x)
  
  x.new <- paste0(x, collapse = ", ")
  x.new <- paste0("list(", x.new, ")")
  x.expr <- parse(text = x.new)
  
  ## Compiled version based on inline package
  ## Non-compiled version based on with() and eval()
  if(compile) {
    
    # Do the replacement to obtain C syntax
    x <- replaceOperation("^", "pow", x)
    x <- replaceSymbols(innames, paste0("x[", (1:length(innames))-1, "+i* *k]"), x)
    names(x) <- paste0("y[", (1:length(outnames)) - 1, "+i* *l]")
    
    # Paste into equation
    x <- x[x != "0"]
    expr <- paste(names(x), "=", x, ";")
    
    # Put equation into C function
    if(is.null(modelname)) {
      funcname <- paste0("funC0_", paste(sample(c(0:9, letters), 8, replace = TRUE), collapse = ""))
    } else {
      funcname <- modelname
    }
    body <- paste(
      "#include <R.h>\n", 
      "#include <math.h>\n", 
      "void", funcname, "( double * x, double * y, int * n, int * k, int * l ) {\n",
      "for(int i = 0; i< *n; i++) {\n",
      paste(expr, collapse="\n"),
      "\n}\n}"
    )
    
    filename <- paste(funcname, "c", sep = ".")
    sink(file = filename)
    cat(body)
    sink()
    system(paste0(R.home(component="bin"), "/R CMD SHLIB ", filename))
    .so <- .Platform$dynlib.ext
    dyn.load(paste0(funcname, .so))
    
    # Generate the C function by the inline package
    #myCfun <- inline::cfunction(sig=c(x = "double", y = "double", n = "integer", k = "integer", l = "integer"),
    #                            body=body,
    #                            language="C",
    #                            convention=".C"
    #)
    
    
    # Generate output function
    myRfun <- function(x, attach.input = FALSE) {
      
      # Translate the list into matrix and then into vector
      M <- do.call(rbind, x[innames])
      if(length(M) == 0) M <- matrix(0)
      x <- as.double(as.vector(M))
      
      # Get integers for the array sizes
      n <- as.integer(dim(M)[2])
      k <- as.integer(length(innames))
      if(length(k) == 0) k <- as.integer(0)
      l <- as.integer(length(outnames))
      
      
      # Initialize output vector
      y <- double(l*n)
      
      # Evaluate C function and write into matrix
      loadDLL(func = funcname, cfunction = funcname)
      out <- matrix(.C(funcname, x = x, y = y, n = n, k = k, l = l)$y, nrow=length(outnames), ncol=n)
      rownames(out) <- outnames
      
      rownames(M) <- innames
      if(attach.input)
        out <- rbind(M, out)
        
      
      return(t(out))    
      
    }
    
    
    
  } else {
    
    # Generate output function
    myRfun <- function(x, attach.input = FALSE) {
      
      # Translate the list into matrix and then into vector
      M <- do.call(rbind, x[innames])
      if(length(M) == 0) M <- matrix(0)
      
      out.list <- with(x, eval(x.expr))
      out.matrix <- do.call(cbind, out.list)
      colnames(out.matrix) <- outnames
      rownames(out.matrix) <- NULL
      
      rownames(M) <- innames
      if(attach.input)
        out.matrix <- cbind(t(M), out.matrix)
      
      return(out.matrix)
      
    }
    
  }
  
  
  attr(myRfun, "equations") <- x
  
  return(myRfun)
  
}


