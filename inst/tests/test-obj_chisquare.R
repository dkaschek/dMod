# # context("SteadyStates")
# # test_that("steady_states_are_steady", {
#   
#   #-!Start example code
#   library(dMod)
#   library(dplyr)
#   # library(ggplot2)
#   setwd(tempdir())
#   devtools::load_all("~/Promotion/Promotion/Software/dMod/")
#   # .. Reactions -----
#   f <- NULL
#   f <- addReaction(f, 
#                    from = "Enz + Sub", 
#                    to = "Compl", 
#                    rate = "k1*Enz*Sub",
#                    description = "production of complex")
#   f <- addReaction(f, 
#                    from = "Compl", 
#                    to = "Enz + Sub", 
#                    rate = "k2*Compl",
#                    description = "decay of complex")
#   f <- addReaction(f, 
#                    from = "Compl", 
#                    to = "Enz + Prod", 
#                    rate = "k3*Compl",
#                    description = "production of product")
#   f <- addReaction(f, 
#                    from = "Enz", 
#                    to = ""     , 
#                    rate = "k4*Enz",
#                    description = "enzyme degradation")
#   
#   # ODE model
#   model <- odemodel(f, modelname = "enzymeKinetics")
#   
#   # Prediction function
#   x <- Xs(model)
#   
#   # .. Observables -----
#   observables <- eqnvec(
#     product = "Prod", 
#     substrate = "(Sub + Compl)", 
#     enzyme = "(Enz + Compl)"
#   )
#   
#   # Generate observation functions
#   g <- Y(observables, x, compile = TRUE, modelname = "obsfn", attach.input = FALSE)
#   
#   # .. Error model -----
#   err <- eqnvec(
#     product = "s0_prod", 
#     substrate = "sqrt((substrate * srel_substrate)^2)", 
#     enzyme = "0.2"
#   )
#   
#   # Generate observation functions
#   e <- Y(err, g, compile = TRUE, modelname = "errfn", attach.input = FALSE)
#   
#   # .. Parameters -----
#   innerpars <- getParameters(e,g,x)
#   # Identity transformation
#   trafo <- repar("x~x", x = innerpars)
#   # Set some initial value parameters
#   trafo <- repar("x~0", x = c("Compl", "Prod"), trafo)
#   # Explicit log-transform of all parameters
#   trafo <- repar("x~exp(x)", x = innerpars, trafo)
#   
#   ## Split transformation into two
#   trafo1 <- trafo2<- trafo
#   
#   # Set the degradation rate in the first condition to 0
#   trafo1["k4"] <- "0"
#   
#   # Generate parameter transformation functions
#   p <- NULL
#   p <- p + P(trafo1, condition = "noDegradation")
#   p <- p + P(trafo2, condition = "withDegradation")
#   
#   # Initialize with randomly chosen parameters
#   set.seed(1234)
#   outerpars <- getParameters(p)
#   pouter <- structure(rnorm(length(outerpars), -2, .5), names = outerpars)
#   times <- 0:100
#   
#   # .. Data -----
#   data <- datalist(
#     noDegradation = data.frame(
#       name = c("product", "product", "product", "substrate", "substrate", "substrate"),
#       time = c(0, 25, 100, 0, 25, 100),
#       value = c(0.0025, 0.2012, 0.3080, 0.3372, 0.1662, 0.0166),
#       sigma = 0.02),
#     withDegradation = data.frame(
#       name = c("product", "product", "product", "substrate", "substrate", "substrate", "enzyme", "enzyme", "enzyme"),
#       time = c(0, 25, 100, 0, 25, 100, 0, 25, 100),
#       value = c(-0.0301,  0.1512, 0.2403, 0.3013, 0.1635, 0.0411, 0.4701, 0.2001, 0.0383),
#       sigma = 0.02)
#   )
#   plot((g*x*p)(times, pouter), data)
#   timesD <- sort(unique(unlist(lapply(data, function(d) d$time))))  
#   
#   d <- as.data.frame(data)
#   
#   
#   # -------------------------------------------------------------------------#
#   # Tests ----
#   # -------------------------------------------------------------------------#
#   # .. 1 gxp, data -----
#   obj1 <- normL2(data, g*x*p) 
#   tofix <- NULL
#   par1 <- pouter[setdiff(names(pouter), tofix)]
#   fix1 <- pouter[tofix]
#   val1 <- obj1(par1, fixed = fix1)
#   fit1 <- trust(obj1, par1, fixed = fix1, rinit = 1, rmax = 10)
#   # .. 2 gxp, data, e -----
#   obj2 <- normL2(data, g*x*p, e) 
#   tofix <- NULL
#   par2 <- pouter[setdiff(names(pouter), tofix)]
#   fix2 <- pouter[tofix]
#   val2 <- obj2(par2, fixed = fix2)
#   fit2 <- trust(obj2, par2, fixed = fix2, rinit = 1, rmax = 10)
#   # .. 3 gxp, data sigma NA, e -----
#   dat3 <- d %>% mutate(sigma = NA) %>% as.data.frame() %>% as.datalist()
#   obj3 <- normL2(dat3, g*x*p, e) 
#   tofix <- NULL
#   par3 <- pouter[setdiff(names(pouter), tofix)]
#   fix3 <- pouter[tofix]
#   val3 <- obj3(par3, fixed = fix3)
#   fit3 <- trust(obj3, par3, fixed = fix3, rinit = 1, rmax = 10)
#   # .. 4 gxp, data sigma NA LLOQ, e -----
#   dat4 <- d %>% mutate(sigma = NA, lloq = case_when(name == "substrate" ~ 0.2, name == "enzyme" ~ 0.1, name == "product" ~ -Inf)) %>% as.data.frame() %>% as.datalist()
#   obj4 <- normL2(dat4, g*x*p, e) 
#   tofix <- NULL
#   par4 <- pouter[setdiff(names(pouter), tofix)]
#   fix4 <- pouter[tofix]
#   val4 <- obj4(par4, fixed = fix4)
#   fit4 <- trust(obj4, par4, fixed = fix4, rinit = 1, rmax = 10)
#   # .. 5 constraint -----
#   obj5 <- constraintL2(structure(rep(0, length(pouter)), names = names(pouter)), 10)
#   tofix <- NULL
#   par5 <- pouter[setdiff(names(pouter), tofix)]
#   fix5 <- pouter[tofix]
#   val5 <- obj5(par5, fixed = fix5)
#   # .. 6 sumobjlist -----
#   val1 + val2
#   val1 + val5
#   fit1 + val5
#   # .. 7 as.parframe -----
#   msfit <- mstrust(obj1, par1, fixed = fix1, fits = 10, sd = 1)
#   pf <- as.parframe(msfit)
#   
#   
#   # -------------------------------------------------------------------------#
#   # Define tests ----
#   # -------------------------------------------------------------------------#
#   #-!End example code
#   
#   
# # })