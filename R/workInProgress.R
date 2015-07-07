similarity <- function(fitlist, x, times, fixed, deriv = FALSE) {
  
  combinations <- combn(1:nrow(fitlist), 2)
  out <- do.call(rbind, lapply(1:ncol(combinations), function(j) {
    
    i1 <- combinations[1, j]
    i2 <- combinations[2, j]
    
    par1 <- unlist(fitlist[i1, -1])
    par2 <- unlist(fitlist[i2, -1])
    
    pred1 <- wide2long(x(times, par1, fixed = fixed, deriv = deriv))
    pred2 <- wide2long(x(times, par2, fixed = fixed, deriv = deriv))
    
    rss <- sum(((pred1$value - pred2$value)/(pred1$value + pred2$value + 1))^2)

    data.frame(i1 = c(i1, i2), i2 = c(i2, i1), rss = rss)    
    
  }))
  
  return(out)
  
}


