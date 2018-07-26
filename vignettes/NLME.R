library(dMod)
library(nlme)
library(dplyr)
library(ggplot2)
setwd("/tmp")

x <- eqnlist() %>% 
  addReaction("A", "B", "kA*A", "Decay of A") %>%
  addReaction("B", "0", "kB*B", "Decay of B") %>%
  odemodel() %>% Xs()


n_individuals <- 10

covtable <- data.frame(
  ID = 1:n_individuals,
  WGT = round(exp(rnorm(n_individuals, log(70), .2)))
)
rownames(covtable) <- do.call(function(...) paste(..., sep = "_"), covtable)

parameters <- getParameters(x)
p <- eqnvec() %>% 
  define("x~x", x = parameters) %>%
  define("B~0") %>%
  insert("x~exp(x)", x = parameters) %>%
  branch(covtable) %>%
  insert("k ~ k + beta_k * log(WGT/70)", WGT = WGT, k = c("kA", "kB")) %>%
  P()

estimate <- getParameters(p)
ini <- structure(rep(-1, length(estimate)), names = estimate)
times <- seq(0, 10, .1)
timesD <- c(0, 1, 2, 4, 7, 10)

(x*p)(times, ini) %>% plot()


# Simualte data
data <- do.call(rbind, lapply(1:n_individuals, function(i) {
  
  ini[1] <- rnorm(1, ini[1], .5)
  ini[2] <- rnorm(1, ini[2], .5)
  
  (x*p)(timesD, ini, deriv = FALSE, conditions = rownames(covtable)[i]) %>% 
    as.data.frame() %>% 
    mutate(ID = covtable$ID[i], WGT = covtable$WGT[i]) %>%
    mutate(sigma = 0.02) %>%
    mutate(value = rnorm(length(value), value, sigma)) %>%
    mutate(ytype = as.numeric(name))
  
}))

data %>% as.datalist(split.by = c("ID", "WGT")) %>% plot()

model <- modelNLME(x*p, covtable)

fit <- nlme::nlme(value~model(time, name, A, kA, beta_kA, kB, beta_kB, ID, WGT), data = data,
                  fixed = A+kA+beta_kA+kB+beta_kB ~ 1, random = A+kA ~ 1 | ID, 
                  weights = varIdent(form = ~ 1 | name),
                  start = ini)


prediction <- fit %>% coef() %>% split(f = 1:10) %>% 
  lapply(function(pars) (x*p)(times, unlist(pars), conditions = rownames(covtable)[as.numeric(rownames(pars))])[[1]]) %>%
  as.prdlist(names = rownames(covtable))

plot(prediction, as.datalist(data, split.by = "condition")) + facet_wrap(~name*condition, scales = "free") 

# Set up function for saemix
model <- modelSAEMIX(x*p, cores = 1)


saemix.data <- saemixData(data, 
                          name.group = "condition", 
                          name.predictors = c("time", "name"), 
                          name.response = "value",
                          name.ytype = "ytype")

ini <- c(A = 0, kA = -2, kB = 0)
saemix.model <- saemixModel(model, 
                            description = "A -> B -> 0", 
                            psi0 = matrix(ini, nrow = 1), 
                            name.modpar = names(ini),
                            transform.par = rep(0, length(ini)),  
                            error.model = c("combined"),
                            error.init = c(.1, .1),
                            verbose = TRUE)

myfit <- saemix(saemix.model, saemix.data)
