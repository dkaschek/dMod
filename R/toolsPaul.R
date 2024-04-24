## Flussempfindlichkeitsanalyse

#' Scale dependent flux sensitivity Analysis
#' @description This function performs the scale dependent flux sensitivity analysis for a simulated system. It prints a pdf-file to the working directory for each condition of the system.
#' @param frame Reactionlist of the model
#' @param x Xs function of the model
#' @param obstates The states under observation whose changes are analysed
#' @param bestfit The parameter values which are used for the simulation
#' @param time Array of observed timepoints
#' @param p P function of the system
#' @param bound Absolute value which a flux sensitive component has to reach at any point in time to be shown in the outputplots
#' @param figurename Name of the saved pdf
#' @param timeunit Unit of the time at the x-axis of the plots
#' @param time,width Dimension of the plots
#' @param linewidth Linewidth of the coloured fluxes in the plots
#' @example inst/examples/scaledependent.R
#' @export

scaledependent <- function(frame, x, obstates, bestfit, time, p, bound=0.001, figurename= "scaledependent", timeunit = "h", width = 16, height = 9, linewidth = 2){
  library(ggplot2)
  library(dMod)
  library(deSolve)
  library(cowplot)
  library(cmocean)
  library(latex2exp)

  Flusse <- getFluxes(frame)
  states <- names(Flusse)

  #Flusse zu Verabeitung in Form bringen
  realfluxes <- c()
  for (fluss in Flusse){
    realfluxes <- c(realfluxes, c(paste0(fluss, collapse = "")))
  }
  realfluxes <- structure(realfluxes, names = names(Flusse))

  vector_collection <- NULL
  for (state in obstates){
    #Zuweisung fur observable-Funktion
    myvector <- structure(c(realfluxes[[state]]), names = c(paste0("Fluss", state)))
    vector_collection <- c(vector_collection, myvector)

    #Generierung der abgeleiteten Flusse
    myvec <- c()
    myvecnames <- c()
    for (i in seq(1, length(states))){
      name <- paste0("Fluss", state, states[i])
      deriv <- paste0(deparse(D(parse(text=paste0(Flusse[[state]], collapse = "")), states[i])), collapse = "")
      derivflux <- paste("(", deriv, ")/", states[i])
      assign(name, derivflux)
      #Vorbereitung zur Normierung
      myvec <- c(myvec, derivflux)
      myvecnames <- c(myvecnames, paste0("Fluss", state, states[i]))
    }
    myvec <- structure(myvec, names = myvecnames)
    vector_collection <- c(vector_collection, myvec)
    #Normierung der Flusse
    mysum <- NULL
    for (i in myvec) {
      if(is.null(mysum)){
        mysum <- paste0("((", i,")^2)^0.5")
      } else mysum <- paste0(mysum, "+ ((", i,")^2)^0.5")
    }
    sum <- structure(mysum, names = c(paste0("sum", state)))
    vector_collection<- c(vector_collection, sum)
  }
  observables <- as.eqnvec(vector_collection)
  g <- Y(observables, frame)

  fit <- (g*x*p)(time, bestfit)

  conditions <- getConditions(p)

  #Plotfunktion
  lapply(conditions, function(condition){
    mydata <- as.data.frame(fit[[condition]])

    for (state in obstates){
      for (fluss in states){
        mydata[paste0("Fluss", state, fluss, "_norm")] <- mydata[paste0("Fluss", state, fluss)]/mydata[paste0("sum", state)]
    }}

    for (name in obstates){
      myplots <- lapply(1:length(states), function(i){ # hier war mclapply
        mystate <- states[i]
        if (is.nan(mydata[paste0("Fluss", name, mystate, "_norm")][1,1])){
          mydata[paste0("Fluss", name, mystate, "_norm")][1,1] <- 0}
        if(any(mydata[paste0("Fluss", name, mystate, "_norm")][-1, ] >= bound) || any(mydata[paste0("Fluss", name, mystate, "_norm")][-1, ] <= -1*bound)){
          p <- ggplot(mydata, aes(time, eval(parse(text=paste0("Fluss", name))))) +
            geom_line(aes(color = eval(parse(text=paste0("Fluss", name, mystate, "_norm")))), linewidth = linewidth) +
            theme_dMod()+ theme(panel.grid = element_blank(), plot.title = element_text(hjust = 0.5, size = 35))+ theme(axis.text = element_text(size = 25), axis.title = element_text(size = 25)) +
            labs(title = TeX(paste0("$\\chi_{", mystate, "}$")), x = paste0("time [", timeunit, "]"), y = TeX(paste0("$v_{",name,"}$"))) +
            scale_color_gradientn(colours = head(tail(cmocean("balance")(512), -40), -40),
                                  n.breaks= 4, limits = c(-1, 1),
                                  guide = "none")#guide_colorbar(title = NULL, barwidth = 1.5, barheight = 15))
          return(p)
        }
      })

      myplots <- Filter(function(x) !is.null(x), myplots)
      myplots[[length(myplots)]] <- myplots[[length(myplots)]] + scale_color_gradientn(colours = head(tail(cmocean("balance")(512), -40), -40),
                                                                                       n.breaks= 4, limits = c(-1, 1),guide = guide_colorbar(title = NULL, barwidth = 2.5, barheight = 25))

      #size <- ceiling(sqrt(length(myplots))) # fur Dimension des Plotrasters
      #myplots <- c(list(NULL), myplots) # Titelzeile leer lassen
      title <- ggdraw() + draw_label(TeX(paste0("$v_{$", name, "}$   ", condition)), fontface = 'bold', x = 0.45, hjust = 0, size = 40) + theme(plot.margin = margin(0, 0, 0, 7))
      plot <- plot_grid(plotlist = myplots, ncol = length(myplots), nrow = 1, rel_widths = c(rep(1, length(myplots)-1), 1.2))
      assign(paste0("p", name), plot_grid(title, plot, ncol = 1, rel_heights = c(0.1, 1)))
    }

    plotcols <- lapply(obstates, function(name){
      eval(parse(text = paste0("p", name)))
    })

    cairo_pdf(filename = paste0(figurename, condition,".pdf"), onefile = TRUE, width = width, height = height)
    for (plotcol in plotcols){
    print(plotcol)
    }
    dev.off()

  })


  gc()

}


#' Scale invariant flux sensitivity Analysis
#' @description This function performs the scale invariant flux sensitivity analysis for a simulated system. It prints a pdf-file to the working directory for each condition of the system.
#' @param frame Reactionlist of the model
#' @param x Xs function of the Model
#' @param obstates The states under observation whose changes are analysed
#' @param bestfit The parameter values which are used for the simulation
#' @param time Array of observed timepoints
#' @param p P function of the system
#' @param bound Absolute value which a flux sensitive component has to reach at any point in time to be shown in the outputplots
#' @param figurename Name of the saved pdf
#' @param timeunit Unit of the time at the x-axis of the plots
#' @param time,width Dimensions of the plots
#' @param linewidth Linewidth of the coloured fluxes in the plots
#' @example inst/examples/scaleinvariant.R
#' @export

scaleinvariant <- function(frame, x, obstates, bestfit, time, p, bound=0.001, figurename= "scaleinvariant", timeunit = "h", width = 16, height = 9, linewidth = 2){
  library(ggplot2)
  library(dMod)
  library(deSolve)
  library(cowplot)
  library(cmocean)
  library(latex2exp)

  Flusse <- getFluxes(frame)
  states <- names(Flusse)

  #Flusse zu Verabeitung in Form bringen
  realfluxes <- c()
  for (fluss in Flusse){
    realfluxes <- c(realfluxes, c(paste0(fluss, collapse = "")))
  }
  realfluxes <- structure(realfluxes, names = names(Flusse))

  vector_collection <- NULL
  for (state in obstates){
    #Zuweisung fur observable-Funktion
    myvector <- structure(c(realfluxes[[state]]), names = c(paste0("Fluss", state)))
    vector_collection <- c(vector_collection, myvector)

    #Generierung der abgeleiteten Flusse
    myvec <- c()
    myvecnames <- c()
    for (i in seq(1, length(states))){
      name <- paste0("Fluss", state, states[i])
      deriv <- paste0(deparse(D(parse(text=paste0(Flusse[[state]], collapse = "")), states[i])), collapse = "")
      derivflux <- paste("(", deriv, ")")
      assign(name, derivflux)
      #Vorbereitung zur Normierung
      myvec <- c(myvec, derivflux)
      myvecnames <- c(myvecnames, paste0("Fluss", state, states[i]))
    }
    myvec <- structure(myvec, names = myvecnames)
    vector_collection <- c(vector_collection, myvec)
  }
  observables <- as.eqnvec(vector_collection)
  g <- Y(observables, frame)

  fit <- (g*x*p)(time, bestfit)

  conditions <- getConditions(p)

  #Plotfunktion
  lapply(conditions, function(condition){
    mydata <- as.data.frame(fit[[condition]])


    for (state  in states){
      mydata[paste0(state, "Punkt")] <- c(c(0), diff(mydata[[paste0(state)]]))/ c(c(1), diff(mydata[["time"]]))
      mydata[paste0("Fluss", state, "Punkt")] <- c(c(0), diff(mydata[[paste0("Fluss", state)]]))/ c(c(1), diff(mydata[["time"]]))
    }


    for (obstate in obstates){
      mydata[paste0("sum", obstate)] <- rep(0, length(mydata[["time"]]))
      for (state  in states){
        mydata[paste0("Anteil", obstate, state)] <- mydata[paste0("Fluss", obstate, state)]*mydata[paste0(state, "Punkt")]/mydata[paste0("Fluss", obstate)]
        mydata[paste0("sum", obstate)] <- mydata[paste0("sum", obstate)] + sqrt(mydata[paste0("Anteil", obstate, state)]**2)
      }
      #Normierung
      for (state  in states){
        mydata[paste0("Anteil", obstate, state)] <- mydata[paste0("Anteil", obstate, state)]/mydata[paste0("sum", obstate)]
        #NaN aussortieren
        for (index in which(sapply(mydata[paste0("Anteil", obstate, state)], function(x) is.nan(x)))){
          mydata[paste0("Anteil", obstate, state)][index, 1] <- 0}
      }}

    for (name in obstates){
      myplots <- lapply(1:length(states), function(i){ # hier war mclapply
        mystate <- states[i]
        if (is.nan(mydata[paste0("Anteil", name, mystate)][1,1])){
          mydata[paste0("Anteil", name, mystate)][1,1] <- 0}
        if(any(mydata[paste0("Anteil", name, mystate)][-1, ] >= bound) || any(mydata[paste0("Anteil", name, mystate)][-1, ] <= -1*bound)){
          p <- ggplot(mydata, aes(time, eval(parse(text=paste0("Fluss", name))))) +
            geom_line(aes(color = eval(parse(text=paste0("Anteil", name, mystate)))), linewidth = linewidth) +
            theme_dMod()+ theme(panel.grid = element_blank(), plot.title = element_text(hjust = 0.5, size = 35))+ theme(axis.text = element_text(size = 25), axis.title = element_text(size = 25)) +
            labs(title = TeX(paste0("$\\rho_{", mystate, "}$")), x = paste0("time [", timeunit, "]"), y = TeX(paste0("$v_{",name,"}$"))) +
            scale_color_gradientn(colours = head(tail(cmocean("balance")(512), -40), -40),
                                  n.breaks= 4, limits = c(-1, 1),
                                  guide = "none")#guide_colorbar(title = NULL, barwidth = 1.5, barheight = 15))
          return(p)
        }
      })

      myplots <- Filter(function(x) !is.null(x), myplots)
      myplots[[length(myplots)]] <- myplots[[length(myplots)]] + scale_color_gradientn(colours = head(tail(cmocean("balance")(512), -40), -40),
                                                                                       n.breaks= 4, limits = c(-1, 1),guide = guide_colorbar(title = NULL, barwidth = 2.5, barheight = 25))

      #size <- ceiling(sqrt(length(myplots))) # fur Dimension des Plotrasters
      #myplots <- c(list(NULL), myplots) # Titelzeile leer lassen
      title <- ggdraw() + draw_label(TeX(paste0("$v_{$", name, "}$   ", condition)), fontface = 'bold', x = 0.45, hjust = 0, size = 40) + theme(plot.margin = margin(0, 0, 0, 7))
      plot <- plot_grid(plotlist = myplots, ncol = length(myplots), nrow = 1, rel_widths = c(rep(1, length(myplots)-1), 1.2))
      assign(paste0("p", name), plot_grid(title, plot, ncol = 1, rel_heights = c(0.1, 1)))
    }

    plotcols <- lapply(obstates, function(name){
      eval(parse(text = paste0("p", name)))
    })

    cairo_pdf(filename = paste0(figurename, condition,".pdf"), onefile = TRUE, width = width, height = height)
    for (plotcol in plotcols){
      print(plotcol)
    }
    dev.off()

  })


  gc()

}


