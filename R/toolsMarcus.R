#' Compute a differentiable box prior
#' 
#' @param p parameter vector
#' @param mu vector of means of boxes
#' @param sigma half box width
#' @param k shape of box; if 0 a quadratic prior is obtained, the higher k the more box shape, gradient at border of the box (-sigma, sigma) is equal to sigma*k
#' @return list with entries: value (numeric, the weighted residual sum of squares), 
#' gradient (numeric, gradient) and 
#' hessian (matrix of type numeric).
constraintExp2 <- function(p, mu, sigma = 1, k = 0.05, kmin=1e-5) {
  par <- names(mu)
  t <- p[par]
  s <- sigma
  k <- sapply(k, function(ki){
    if(ki < kmin){
      kmin
    } else ki
  })
  
  gr <- rep(0, length(p)); names(gr) <- names(p)
  hs <- matrix(0, length(p), length(p), dimnames = list(names(p), names(p)))
  
  val <- sum(0.5*(exp(k*((t-mu)/s)^2)-1)/(exp(k)-1))
  gr <- (k*(t-mu)/(s^2)*exp(k*((t-mu)/s)^2)/(exp(k)-1))
  diag(hs)[par] <- k/(s*s)*exp(k*((t-mu)/s)^2)/(exp(k)-1)*(1+2*k*(t-mu)/(s^2))
  
  dP <- attr(p, "deriv")
  if(!is.null(dP)) {
    gr <- as.vector(gr%*%dP); names(gr) <- colnames(dP)
    hs <- t(dP)%*%hs%*%dP; colnames(hs) <- colnames(dP); rownames(hs) <- colnames(dP)
  }
  
  out <- list(value=val,gradient=gr,hessian=hs)
  class(out) <- c("obj", "list")
  
  return(out)
  
}