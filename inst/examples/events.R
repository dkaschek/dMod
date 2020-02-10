library(dMod)
library(dplyr)
library(ggplot2)
setwd(tempdir())

#' Check function to check sensitivities
#'
#' @param p pars for the prediction function
#' @param whichpar names of pars to check
#' @param cond indexing for condition
#' @param step stepsize for finite difference
checkSensitivities <- function(p, whichpar, cond = 1, step = 0.1) {
  h <- rep(0, length(p))
  h[which(names(p) == whichpar)] <- step
  
  
  M1 <-  y(times, p, deriv = TRUE)[[cond]]
  M2 <-  y(times, p + h, deriv = TRUE)[[cond]]
  M3 <- attr(y(times, p, deriv = TRUE)[[cond]], "deriv")
  print(colnames(M3))
  
  S1 <- cbind(time = M1[match(times, M1[, 1]), 1], (M2[match(times, M2[, 1]),-1, drop = FALSE] - M1[match(times, M1[, 1]),-1, drop = FALSE])/step)
  print(colnames(S1))
  S2 <- cbind(time = M1[match(times, M1[, 1]), 1], M3[match(times, M3[, 1]), grep(paste0(".", whichpar), colnames(M3), fixed = TRUE), drop = FALSE])
  print(colnames(S2))
  colnames(S2) <- colnames(S1)
  
  
  out <- prdlist(numeric = S1, sens = S2)

  return(out)
  
}

## check with root-triggered events
model <- eqnlist() %>%
  addReaction("A", "B", "kon*A") %>%
  addReaction("B", "A", "koff*B") %>%
  addReaction("B", "0", "degrad*B")%>%
  addReaction("0", "kon", "kbase", "Event state") %>%
  addReaction("kon", "0", "degrad*kon", "Event state") %>%
  odemodel(
    events = rbind(
      #data.frame(var = "B", time = "t_thres", value = "1", root = NA, method = "replace", stringsAsFactors = FALSE),
      data.frame(var = "A", time = "t_A_thres", value = "1", root = "A - A_thres", method = "replace", stringsAsFactors = FALSE)
    ),
    estimate = c("A_thres", "A", "t_A_thres", "kon", "koff", "degrad")
  ) 
x <- model %>% Xs(optionsOde = list(method = "lsoda"), optionsSens = list(method = "lsoda", rtol = 1e-10, atol = 1e-10))

innerpars <- getParameters(x)

p <- eqnvec() %>%
  define("x~x", x = innerpars) %>%
  define("x~0", x = "B") %>%
  define("x~1", x = "A") %>%
  insert("x~exp(x)", x = innerpars) %>%
  P()

outerpars <- getParameters(p)

set.seed(33)
pouter <- structure(rnorm(length(outerpars), -1), names = outerpars)
pouter["A_thres"] <- log(0.8)
# pouter["t_thres"] <- log(0.2)
# pouter["val_A"] <- log(1)
times <- seq(0, 3, .01)

pouter %>% (x*p)(times = times, deriv = TRUE) %>% plot()
# pouter %>% (x*p)(times = times) %>% getDerivs() %>% plot()

y <- x*p

#pdf("~/root_events.pdf")

for (i in 1:length(pouter)) {
  out <- checkSensitivities(pouter, names(pouter)[i], 3, .00000001) %>% as.prdlist
  print(plotPrediction(out) + ggtitle(names(pouter)[i]))
  
}

#dev.off()



(GR-EMAX*((Ce+yps)^hill/((Ce+yps)^hill+EC50^hill)))*Gcure
(Fcure-1)*Gcure


## check with root-triggered events
model2 <- eqnvec(
  PL = "(GR - EMAX*exp(-k*time))*Gcure",
  Gcure = "(Fcure - 1)*Gcure",
  Fcure = "0"
) %>%
  odemodel(
    events = rbind(
      #data.frame(var = "B", time = "t_thres", value = "1", root = NA, method = "replace", stringsAsFactors = FALSE),
      data.frame(var = "Fcure", time = "tPLcure", value = "0", root = "PL - PLcure", method = "replace", stringsAsFactors = FALSE)
    ),
    estimate = c("PLcure", "Fcure", "tPLcure", "EMAX")
  ) 
x <- model2 %>% Xs(optionsOde = list(method = "lsoda"), optionsSens = list(method = "lsoda", rtol = 1e-10, atol = 1e-10))

innerpars <- getParameters(x)

p <- eqnvec() %>%
  define("x~x", x = innerpars) %>%
  define("x~5", x = "PL") %>%
  define("x~1", x = c("Fcure", "Gcure")) %>% 
  P()

outerpars <- getParameters(p)

set.seed(33)
pouter <- structure(rnorm(length(outerpars), -1), names = outerpars)
pouter["GR"] <- 0.03
pouter["PLcure"] <- -5
pouter["EMAX"] <- 10.3
pouter["k"] <- 1 
pouter["tPLcure"] <- 0

times <- seq(0, 10, .1)

pouter %>% (x*p)(times = times, deriv = TRUE) %>% plot()
# pouter %>% (x*p)(times = times) %>% getDerivs() %>% plot()

y <- x*p

#pdf("~/root_events.pdf")

for (i in 1:length(pouter)) {
  out <- checkSensitivities(pouter, names(pouter)[i], 1, .000001) %>% as.prdlist
  print(plotPrediction(out) + ggtitle(names(pouter)[i]))
  
}









## Do the check with method = "replace"
model <- eqnlist() %>%
  addReaction("A+A", "B", "kon*A^2") %>%
  addReaction("B", "A", "koff*B") %>%
  addReaction("B", "0", "degrad*B")%>%
  addReaction("0", "kon", "kbase", "Event state") %>%
  addReaction("kon", "0", "degrad*kon", "Event state") %>%
  odemodel(
    events = rbind(
      data.frame(var = "kon", time = c(0, "toff1"), value = c(1, "kmax"), method = "replace", stringsAsFactors = FALSE),
      data.frame(var = "B", time = c(0, "toff2"), value = c(1., "kmax"), method = "replace", stringsAsFactors = FALSE)
    )
  ) 
x <- model %>% Xs()

innerpars <- getParameters(x)

p <- eqnvec() %>%
  define("x~x", x = innerpars) %>%
  define("x~0", x = "B") %>%
  define("x~1", x = "A") %>%
  insert("x~exp(x)", x = innerpars) %>%
  P()

outerpars <- getParameters(p)

set.seed(2)
pouter <- structure(rnorm(length(outerpars), -1), names = outerpars)
times <- seq(0, 7, .01)

pouter %>% (x*p)(times = times) %>% plot()
# pouter %>% (x*p)(times = times) %>% getDerivs() %>% plot()

y <- x*p

for (i in 1:length(pouter)) {
  out <- checkSensitivities(pouter, names(pouter)[i], 3, .000001) %>% as.prdlist
  print(plotPrediction(out) + ggtitle(names(pouter)[i]))
  
}



## Do the check with method = "add"
library(dMod)
library(dplyr)
library(ggplot2)
setwd("/tmp")



model <- eqnlist() %>%
  addReaction("A+A", "B", "kon*A^2") %>%
  addReaction("B", "A", "koff*B") %>%
  addReaction("B", "0", "degrad*B")%>%
  addReaction("0", "kon", "0", "Event state") %>%
  addReaction("kon", "0", "decay*B*kon", "Event state") %>%
  odemodel(
    events = rbind(
      data.frame(var = "kon", time = "te", value = "ve", method = "add"),
      data.frame(var = "kon", time = "tf", value = "ve", method = "add"),
      data.frame(var = "B"  , time = "te", value = "1" , method = "add")
    )
  ) 
x <- model %>% Xs()

innerpars <- getParameters(x)

p <- eqnvec() %>%
  define("x~x", x = innerpars) %>%
  define("x~0", x = c("B")) %>%
  define("x~1", x = c("kon", "A")) %>%
  P()

outerpars <- getParameters(p)

pouter <- structure(runif(length(outerpars)), names = outerpars)
pouter["te"] <- log(2.092)
pouter["tf"] <- log(3.092)
pouter["ve"] <- log(3.5)
times <- seq(0, 4, .01)

pouter %>% (x*p)(times = times) %>% plot()
#pouter %>% (x*p)(times = times) %>% getDerivs() %>% plot()

y <- x*p


for (i in 1:length(pouter)) {
  out <- checkSensitivities(pouter, names(pouter)[i], 1, .001) %>% as.prdlist()
  print(plotPrediction(out) + ggtitle(names(pouter)[i])) 
  
}


## Do the check with method = "multiply"
library(dMod)
library(dplyr)
library(ggplot2)
setwd("/tmp")



model <- eqnlist() %>%
  addReaction("A+A", "B", "kon*A^2") %>%
  addReaction("B", "A", "koff*B") %>%
  addReaction("B", "0", "degrad*B")%>%
  addReaction("0", "kon", "0", "Event state") %>%
  addReaction("kon", "0", "decay*B*kon", "Event state") %>%
  odemodel(
    events = rbind(data.frame(var = "kon", time = "te", value = "ve", method = "multiply"),
                   data.frame(var = "kon", time = "tf", value = "ve", method = "multiply"),
                   data.frame(var = "B"  , time = "te", value = "3.", method = "multiply"))
  ) 
x <- model %>% Xs()

innerpars <- getParameters(x)

p <- eqnvec() %>%
  define("x~x", x = innerpars) %>%
  define("x~0", x = c("B")) %>%
  define("x~1", x = c("kon", "A")) %>%
  P()

outerpars <- getParameters(p)

pouter <- structure(runif(length(outerpars)), names = outerpars)
pouter["te"] <- log(2.092)
pouter["tf"] <- log(3.092)
pouter["ve"] <- log(3.5)
times <- seq(0, 4, .01)

pouter %>% (x*p)(times = times) %>% plot()
#pouter %>% (x*p)(times = times) %>% getDerivs() %>% plot()

y <- x*p


for (i in 1:length(pouter)) {
  out <- checkSensitivities(pouter, names(pouter)[i], 1, .001) %>% as.prdlist()
  print(plotPrediction(out) + ggtitle(names(pouter)[i])) 
  
}



###

ODEs <- c(
  Ad= "-ka*Ad+Fabs*INPUT1",
  Ac = "ka*Ad-CL/Vc*Ac",
  PD = "kin-(1+(EMAX*(Ac/Vc)/(EC50+(Ac/Vc))))*kout*PD",
  INPUT1 = "0"
)

events <- eventlist() %>% 
  addEvent("INPUT1", "ton_INPUT1_1", "xon_INPUT1_1") %>% 
  addEvent("INPUT1", "toff_INPUT1_1", "0")

model <- odemodel(ODEs, events = events, estimate = "ton_INPUT1_1")
attr(model$extended, "events")


x <- model %>% Xs()

innerpars <- getParameters(x)

p <- eqnvec() %>%
  define("x~x", x = innerpars) %>%
  define("x~0", x = c("INPUT1", "Ad", "Ac")) %>%
  define("x~1", x = c("PD")) %>%
  P()

outerpars <- getParameters(p)

pouter <- structure(runif(length(outerpars)), names = outerpars)
pouter["ton_INPUT1_1"] <- 1
pouter["toff_INPUT1_1"] <- 2
pouter["xon_INPUT1_1"] <- 1
times <- seq(0, 4, .01)

pouter %>% (x*p)(times = times) %>% plot()
pouter %>% (x*p)(times = times) %>% getDerivs() %>% plot()

y <- x*p


for (i in 1:length(pouter)) {
  out <- checkSensitivities(pouter, names(pouter)[i], 1, .01) %>% as.prdlist()
  print(plotPrediction(out) + ggtitle(names(pouter)[i])) 
  
}

