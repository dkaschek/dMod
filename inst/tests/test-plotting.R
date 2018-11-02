# context("Plotting")
#   test_that("PlotX functions behave correctly", {
# 
# 
#     #-!Start example code
#     library(dplyr)
#     library(dMod)
#     source("../examples/example_CCD4/setup.R")
#     
#     construct_parframe <- function(pars, n = 20, seed = 12345, samplefun = rnorm) {
#       set.seed(seed)
#       rnd <- samplefun(n*length(pars))
#       mypars <- matrix(rnd, nrow = n)
#       mypars <- `names<-`(as.data.frame(mypars), names(pars))
#       parframe(mypars)
#     }
#     
#     model <- dMod.frame("Paper version", g, x, p, data, NULL, odemodels = list(model0))
#     model <- appendObj(model, parframes = list(construct_parframe(pars)))
#     model <- mutate(model, fits = list(mstrust(obj, parframes, cores = 4)))
#     model <- appendParframes(model)
#     
#     p1 <- plotCombined(model)
#     
#     model <- dMod.frame("Paper version", g, x, p, data, NULL, odemodels = list(model0))
#     model <- appendObj(model)
#     model <- mutate(model,
#                     fixed = list(c(logkb1 = -1.04, logbcar = 1.78, logkb2 = -2.31)))
#     model <- mutate(model,
#                     pars = list(structure(rnorm(length(getParameters(prd))-3),
#                                           names = getParameters(prd)[!getParameters(prd)%in%names(fixed)])),
#                     parframes = list(construct_parframe(pars)))
#     model <- mutate(model, fits = list(mstrust(obj, parframes, cores = 4, fixed = fixed)))
#     model <- appendParframes(model)
#     
#     p2 <- plotCombined(model)
#     
#     unlink(paste0("*.", c("c", "o", "so")))
#     #-!End example code
#     
#     saveRDS(p1, "resources/plotting/p1.rds")
#     p1_saved <- readRDS("resources/plotting/p1.rds")
#     
#     
#     saveRDS(p2, "resources/plotting/p2.rds")
#     p2_saved <- readRDS("resources/plotting/p2.rds")
#     
# 
#     # Define your expectations here
#     expect_identical(as.character(p1), as.character(p1_saved))
#     expect_identical(as.character(p2), as.character(p2_saved))
#   })
