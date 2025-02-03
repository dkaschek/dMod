#' Optimize an Objective Function using Various Optimization Methods using optimx
#'
#' This function extends the existing optimization capabilities by integrating the `optimx` package,
#' allowing the use of various optimization algorithms such as "L-BFGS-B", "BFGS", "Nelder-Mead", etc.
#'
#' @param objfn A function that computes and returns a list with components `value`, `gradient`, and `hessian`.
#'              This represents the objective function to be minimized.
#' @param parinit A numeric vector of initial parameter values.
#' @param method A character string specifying the optimization method. Defaults to "L-BFGS-B".
#'               Available methods include those supported by the `optimx` package, such as "BFGS",
#'               "Nelder-Mead", "L-BFGS-B", etc.
#' @param lower A numeric vector of lower bounds for the parameters (used only by methods that support
#'              box constraints, e.g., "L-BFGS-B"). Defaults to `-Inf`.
#' @param upper A numeric vector of upper bounds for the parameters (used only by methods that support
#'              box constraints, e.g., "L-BFGS-B"). Defaults to `Inf`.
#' @param control A list of control parameters to pass to the optimization algorithm.
#' @param ... Additional arguments to pass to the objective function.
#' 
#' @return A list containing:
#'         - `value`: The value of the objective function at the optimum.
#'         - `gradient`: The gradient at the optimum.
#'         - `hessian`: The Hessian at the optimum.
#'         - `argument`: The optimized parameters.
#'         - `converged`: Logical indicating if the optimizer converged.
#'         - `iterations`: The number of function evaluations.
#' @import optimx
#' @export
optimize <- function(objfn, parinit, method = "L-BFGS-B", lower = -Inf, upper = Inf, control = list(), ...) {
  
  # Sanitize the initial parameters.
  sanePars <- sanitizePars(parinit, list(...)$fixed)
  parinit <- sanePars$pars
  
  # Ensure lower/upper bounds are vectors with proper names.
  if (length(lower) == 1) {
    lower <- rep(lower, length(parinit))
    names(lower) <- names(parinit)
  } else if (is.null(names(lower))) {
    names(lower) <- names(parinit)
  }
  
  if (length(upper) == 1) {
    upper <- rep(upper, length(parinit))
    names(upper) <- names(parinit)
  } else if (is.null(names(upper))) {
    names(upper) <- names(parinit)
  }
  
  # Initialize cache to store previously evaluated parameter results.
  cache <- new.env(hash = TRUE, parent = emptyenv())
  
  # Helper function to generate a unique key for each parameter set.
  generate_key <- function(par) paste0(format(par, digits = 8), collapse = ",")
  
  # Define the function wrappers required by optimx with error handling and caching.
  fn <- function(par, ...) {
    names(par) <- names(parinit)
    key <- generate_key(par)
    
    if (exists(key, envir = cache)) {
      result <- get(key, envir = cache)
    } else {
      result <- try(objfn(par, ...), silent = TRUE)
      assign(key, result, envir = cache)
    }
    
    if (inherits(result, "try-error") || !is.finite(result$value)) {
      return(1e10)  # Large penalty for solver failures
    }
    result$value
  }
  
  gr <- function(par, ...) {
    names(par) <- names(parinit)
    key <- generate_key(par)
    
    if (exists(key, envir = cache)) {
      result <- get(key, envir = cache)
    } else {
      result <- try(objfn(par, ...), silent = TRUE)
      assign(key, result, envir = cache)
    }
    
    if (inherits(result, "try-error") || !all(is.finite(result$gradient))) {
      return(rep(0, length(par)))  # Return zero gradient on failure
    }
    unname(result$gradient)
  }
  
  hess <- function(par, ...) {
    names(par) <- names(parinit)
    key <- generate_key(par)
    
    if (exists(key, envir = cache)) {
      result <- get(key, envir = cache)
    } else {
      result <- try(objfn(par, ...), silent = TRUE)
      assign(key, result, envir = cache)
    }
    
    if (inherits(result, "try-error") || !all(is.finite(result$hessian))) {
      return(diag(length(par)))  # Return identity matrix if Hessian fails
    }
    unname(result$hessian)
  }
  
  # Perform the optimization.
  optim_result <- optimx::optimr(
    par = as.numeric(parinit),
    fn = fn,
    gr = gr,
    hess = hess,
    method = method,
    lower = unname(lower),
    upper = unname(upper),
    control = control,
    ...
  )
  
  # Extract the optimized parameters.
  final_par <- structure(optim_result$par, names = names(parinit))
  
  # Evaluate the objective function at the optimum.
  final_result <- objfn(final_par, ...)
  
  # Attach optimization metadata.
  final_result$argument <- final_par
  final_result$converged <- !as.logical(optim_result$convergence)
  final_result$iterations <- optim_result$counts["function"]
  
  return(final_result)
}
