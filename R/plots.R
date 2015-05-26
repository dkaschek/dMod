plotPrediction <- function(out, ...) {

  require(ggplot2)
  if(any(class(out)=="list")) out <- wide2long(out)
  
  out <- subset(out, ...)
  
  ggplot(out, aes(x=time, y=value, group=condition, color=condition)) + facet_wrap(~name, scales="free") + geom_line()
  
}

plotPredictionCont <- function(out, ...) {
  
  require(ggplot2)
  
  mynames <- c("time", "name", "value", "sigma", "condition")
  
  
  forclist <- lapply(out, function(o) attr(o, "forc"))
  out <- cbind(wide2long.list(out), sigma=NA)
  out <- subset(out, ...)
  targets <- as.character(unique(out$name))
  
  print(targets)
  
  forc <- lbind(lapply(forclist, function(fo) {
    names(fo) <- substr(names(fo), 1, nchar(names(fo))-1)
    print(names(fo))
    data <- do.call(rbind, lapply(targets[targets%in%names(fo)], function(t) {
      print(t)
      data.frame(time = fo[[t]][,1], name = t, value = fo[[t]][,2], sigma=1/sqrt(fo[[paste0("weight", t)]][,2]))
    }))
    
    return(data)
        
  }))
  
  ggplot(rbind(out[, mynames], forc[, mynames]), 
         aes(x = time, 
             y = value, ymin = value - sigma, ymax = value + sigma, 
             group = condition, color = condition, fill=condition)) + 
    facet_wrap(~name, scales = "free") + 
    geom_line(data = out) + 
    geom_line(data = forc, lty=2) + 
    geom_ribbon(data = forc, alpha=0.3, lty=0)
  
  
}


plotCombined <- function (prediction, data, ...) 
{
  require(ggplot2)
  mynames <- c("time", "name", "value", "sigma", "condition")
  if (any(class(prediction) == "list")) 
    prediction <- cbind(wide2long(prediction), sigma = NA)
  if (any(class(data) == "list")) 
    data <- lbind(data)
  
  
  prediction <- subset(prediction, ...)
  data <- subset(data, ...)
  
  
  ggplot(rbind(prediction[, mynames], data[, mynames]), aes(x = time, 
                                                            y = value, ymin = value - sigma, ymax = value + sigma, 
                                                            group = condition, color = condition)) + facet_wrap(~name, 
                                                                                                                scales = "free") + geom_line(data = prediction) + geom_point(data = data) + 
    geom_errorbar(data = data, width = 0)
}

plotData <- function (mydata, ...) 
{
  require(ggplot2)
  if (any(class(mydata) == "list")) 
    mydata <- lbind(mydata)
  mydata <- subset(mydata, ...)
  ggplot(mydata, aes(x = time, y = value, ymin = value - sigma, 
                     ymax = value + sigma, group = condition, color = condition)) + 
    facet_wrap(~name, scales = "free") + geom_point() + geom_errorbar(width = 0)
}


# plotData <- function(mydata) {
#   require(ggplot2)
#   if(any(class(mydata)=="list")) mydata <- lbind(mydata)
#   
#   ggplot(mydata, aes(x=time, y=value, ymin=value-sigma, ymax=value+sigma, group=condition, color=condition)) +
#     facet_wrap(~name, scales="free") +
#     geom_point() +
#     geom_errorbar(width=0)
#   
# }


# plotCombined <- function(prediction, data) {
#   require(ggplot2)
#   mynames <- c("time", "name", "value", "sigma", "condition")
#   if(any(class(prediction)=="list")) prediction <- cbind(wide2long(prediction), sigma=NA)
#   if(any(class(data)=="list")) data <- lbind(data)
#   
#   out <- rbind(prediction[,mynames], data[,mynames])
#   
#   ggplot(out, aes(x=time, y=value, ymin=value-sigma, ymax = value+sigma, group=condition, color=condition)) +
#     facet_wrap(~name, scales="free") +
#     geom_line(data=prediction) +
#     geom_point(data=data) + geom_errorbar(data=data, width=0)
#   
#   
# }


plotArray <- function(x, times, fitlist, data = NULL, fixed=NULL, binwidth=10) {
  
  require("ggplot2")
  require("wq")
  
  n <- dim(fitlist)[1]
  chisquare <- as.data.frame(log10(fitlist[,1]))
  colnames(chisquare) <- "logchisquare"
  
  
  out <- do.call(rbind, lapply(1:n, function(i) {
    myout <- x(times, c(fitlist[i,-1], fixed), deriv=FALSE)
    myout <- wide2long(myout)
    myout <- cbind(myout, logchisquare = log10(fitlist[i,1]))
    return(myout)
  }))
  
  
  P1 <- ggplot(out, aes(x=time, y=value, group=logchisquare, color=logchisquare)) +
    facet_wrap(~name, scales="free") 
  if(!is.null(data)) {
    data <- cbind(data, logchisquare = NA)
    P1 <- P1 + 
      geom_point(aes(x=time, y=value), data = data, color="gray") + 
      geom_errorbar(aes(x=time, y=value, ymin=value-sigma, ymax=value+sigma), data=data, color="gray", width=0)
  }
  P1 <- P1 +
    geom_line(alpha=.3) +
    scale_color_gradientn(colours = rainbow(7)) +
    theme(legend.position="none")
  
  
  #   P3 <- ggplot(out, aes(x=logchisquare, group=logchisquare, fill=logchisquare)) + 
  #     geom_histogram(binwidth=binwidth) + 
  #     scale_fill_gradientn(colours = rainbow(7)) +
  #     theme(legend.position="none") +
  #     xlab("") +
  #     coord_flip()
  
  
  
  P2 <- ggplot(chisquare, aes(x = 1:length(logchisquare), y=logchisquare, color=logchisquare)) +
    geom_point() +
    scale_color_gradientn(colours = rainbow(7)) + 
    theme(legend.position="none") +
    xlab("sorted index") +
    ylab(expression(log[10](chi^2)))
  
  layOut(list(P2, 1, 1),
         list(P1, 1, 2:5))
  
  
}


plotObjective <- function(out) {
  require("ggplot2")
  require("wq")
  
  value <- out$value
  
  gradient <- out$gradient
  npar <- length(gradient)
  names <- factor(paste(names(gradient), 1:npar, sep=", "), levels = paste(names(gradient), 1:npar, sep=", "))
  gradient.data <- data.frame(name = names, value = gradient)
  
  
  
  hessian <- out$hessian
  hessian.data <- data.frame(x = as.factor(1:npar), y=as.factor(rep(1:npar, each=npar)), hessian = as.vector(hessian))
  
  
  P1 <- ggplot(gradient.data, aes(x=name, y= value)) + 
    geom_bar(stat="identity") + 
    coord_flip() + ylab("gradient value") + xlab("parameter name, i") + 
    ggtitle(paste("obj. value:", value)) 
  
  
  P2 <- ggplot(hessian.data, aes(x=x, y=y, z=hessian, fill=hessian)) + 
    geom_tile(color="gray") + scale_fill_gradient2() + xlab("i") + ylab("j") + 
    theme(legend.position=c(0, 1), legend.justification=c(0,1))
  
  layOut(list(P2, 1, 1:2),
         list(P1, 1, 3))
  
  
  
}



plotFitList <- function(fitlist) {
  require(ggplot2)
    
  ggplot(fitlist, aes(x=chisquare, y=value)) + 
    facet_wrap(~name, scales="free") + 
    geom_point()  
  
}

plotProfile <- function(..., maxvalue = 5) {
  
  require(ggplot2)
  
  arglist <- list(...)
  
  data <- do.call(rbind, lapply(1:length(arglist), function(i) {
    proflist <- arglist[[i]]
    subdata <- do.call(rbind, lapply(names(proflist), function(n) {
      
      values <- proflist[[n]][,1]
      zerovalue <- proflist[[n]]["out",1]
      parvalues <- proflist[[n]][,n]
      deltavalues <- values - zerovalue
      
      subset(data.frame(name = n, delta = deltavalues, par = parvalues, proflist = i), delta <= maxvalue)
      
    }))
    return(subdata)
  }))
  
  data$proflist <- as.factor(data$proflist)
  
  
  threshold <- c(1, 2.7, 3.84)
  
  p <- ggplot(data, aes(x=par, y=delta, group=proflist, color=proflist)) + facet_wrap(~name, scales="free_x") + 
    geom_line() + geom_point(aes=aes(size=1), alpha=1/3) +
    geom_hline(yintercept=threshold, lty=2, color="gray") + 
    ylab(expression(paste("CL /", Delta*chi^2))) +
    scale_y_continuous(breaks=c(1, 2.7, 3.84), labels = c("68% / 1   ", "90% / 2.71", "95% / 3.84")) +
    xlab("parameter value")
  
  
  return(p)
  
}

plotPaths <- function(..., whichPar = NULL, sort = FALSE) {
  
  require(ggplot2)
  
  arglist <- list(...)
  
  data <- do.call(rbind, lapply(1:length(arglist), function(i) {
    # choose a proflist
    proflist <- arglist[[i]]
    if(is.null(whichPar)) whichPar <- names(proflist)
    subdata <- do.call(rbind, lapply(whichPar, function(n) {
      # chose a profile
      paths <- proflist[[n]][,-(1:4)]
      values <- proflist[[n]][,1]
      combinations <- expand.grid.alt(whichPar, colnames(paths))
      if(sort) combinations <- apply(combinations, 1, sort) else combinations <- apply(combinations, 1, identity)
      combinations <- combinations[,-which(combinations[1,] == combinations[2,])]
      combinations <- combinations[,!duplicated(paste(combinations[1,], combinations[2,]))]
      
      
        
      
      path.data <- do.call(rbind, lapply(1:dim(combinations)[2], function(j) {
        data.frame(chisquare = values, 
                   name = n,
                   proflist = i,
                   combination = paste(combinations[,j], collapse = " - "),
                   x = paths[,combinations[1,j]],
                   y = paths[,combinations[2,j]])
      }))
      
      return(path.data)
      
    }))
    
    return(subdata)
    
  }))
  
  data$proflist <- as.factor(data$proflist)
  
  
  p <- ggplot(data, aes(x=x, y=y, group=interaction(name, proflist), color=name, lty=proflist)) + 
    facet_wrap(~combination, scales="free") + 
    geom_path() + geom_point(aes=aes(size=1), alpha=1/3) +
    xlab("parameter 1") + ylab("parameter 2") +
    scale_linetype_discrete(name = "profile\nlist") +
    scale_color_discrete(name = "profiled\nparameter")
  
  return(p)
  
}


plotFluxes <- function(out, fluxEquations, pars) {
  
  require(scales)
  
  if(any(class(out)=="list")) out <- wide2long(out)
  
  nFluxes <- length(fluxEquations)
  if(is.null(names(fluxEquations))) names(fluxEquations) <- paste0("reaction", 1:nFluxes)
  fluxEquations <- c(fluxEquations, sum = paste(fluxEquations, collapse="+"))
  
  # Evaluate fluxes
  fluxes <- with(c(as.list(out), as.list(pars)), {
    flux <- do.call(rbind, lapply(1:(nFluxes+1), function(i) {
      ev <- eval(parse(text=fluxEquations[i]))
      nout <- data.frame(time = out$time, 
                         name = names(fluxEquations)[i], 
                         value = ev, 
                         condition = out$condition)
    }))
    return(flux)
  })
  
  fluxes1 <- subset(fluxes, name != "sum")
  fluxes2 <- subset(fluxes, name == "sum")
  
  
  P <- ggplot(out, aes(x=time, y=value, group=name, fill=name)) + facet_wrap(~condition) +
    geom_density(stat="identity", position="stack", alpha=0.3, color="darkgrey", size=0.4, data=fluxes1) +
    geom_line(aes(x=time, y=value, group=NULL, fill=NULL), color="black", data=fluxes2)
  
  
  return(P)
  
}

