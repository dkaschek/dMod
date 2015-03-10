gauss <- function(t, mu, s) {
  zero <- 1e-10
  n <- 50/(s*sqrt(2*pi))
  out <- n * exp(-0.5*((t-mu)/s)^2)
  out[out<zero] <- zero
  #out <- rep(0, length(t))
  return(out)
}
gaussC <- function(t, mu, s) {
  n <- 1/(s*sqrt(2*pi))
  val <- n * exp(-0.5*((t-mu)/s)^2)
  #out[out<zero] <- zero
  gr <- n * exp(-0.5*((t-mu)/s)^2)*(-(t-mu)/s)
  hs <- n * exp(-0.5*((t-mu)/s)^2)*(-(t-mu)/s)^2 +n * exp(-0.5*((t-mu)/s)^2)*(-1/s)
  return(list(value=val,grad=gr,hessian=hs))
}


constraint <- function(pp, par, mu=6, s=0.1) {
  t <- pp[par]
  
  res <- rep(0, length(pp)); names(res) <- names(pp)
  res[par] <- sqrt(0.5)*(t-mu)/s
  J <- matrix(0, nrow=length(res), ncol=length(res)); rownames(J) <- colnames(J) <- names(pp)
  diag(J)[par] <- sqrt(0.5)/s
  
  val <- sum((0.5*((t-mu)/s)^2))
  gr <- rep(0, length(pp)); names(gr) <- names(pp)
  gr[par] <- ((t-mu)/(s^2))
  
  hs <- matrix(0, length(pp), length(pp))
  colnames(hs) <- names(pp); rownames(hs) <- names(pp)
  diag(hs)[par] <- 1/(s*s)
  
  dP <- attr(pp, "deriv")
  if(!is.null(dP)) {
    gr <- as.vector(gr%*%dP); names(gr) <- colnames(dP)
    hs <- t(dP)%*%hs%*%dP; colnames(hs) <- colnames(dP); rownames(hs) <- colnames(dP)
  }
  
  out <- list(value=val,gradient=gr,hessian=hs, res = res, J = J)
  class(out) <- c("obj", "list")
  
  return(out)

}

repWithNames <- function(x, names){
  repnum <- rep(x,length(names))
  names(repnum) <- names
  return(repnum)
}


getSigmaFromProfileList <- function(proflist, alpha=0.32) {
  
  
  data <- do.call(rbind, lapply(names(proflist), function(n) {
    
    
    
    values <- proflist[[n]][,1]
    zerovalue <- proflist[[n]]["out",1]
    parvalues <- (proflist[[n]][,-2])[,n]
    deltavalues <- values - zerovalue
    
    data.frame(name = n, delta = deltavalues, par = parvalues)
    
  }))
  
  #return(data)
  
  

  threshold <- qchisq(1-alpha, 1)
  
  sigma <- lapply(unique(data$name), function(n) {
    
    subdata <-   subset(data, name==n)
    mid <- which.min(subdata$delta)
    subdataLeft <- subdata[1:mid,]
    subdataRight <- subdata[mid:(dim(subdata)[1]),]
    left <- subdataLeft$par[which.min(abs(threshold - subdataLeft$delta))]
    right <- subdataRight$par[which.min(abs(threshold - subdataRight$delta))]
    
    return( (right-left)/2)
    
    
  })
  
  names(sigma) <- unique(data$name)
  
  return(sigma)
  
}
