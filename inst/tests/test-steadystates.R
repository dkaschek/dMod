context("SteadyStates")
test_that("steady_states_are_steady", {

  #-!Start example code
  #-! library(dMod)
  #-! setwd(tempdir())

  reactions <- eqnlist()
  reactions <-   addReaction(reactions, "Tca_buffer", "Tca_cyto", "import_Tca*Tca_buffer", "Basolateral uptake")
  reactions <-   addReaction(reactions, "Tca_cyto", "Tca_buffer", "export_Tca_baso*Tca_cyto", "Basolateral efflux")
  reactions <-   addReaction(reactions, "Tca_cyto", "Tca_canalicular", "export_Tca_cana*Tca_cyto", "Canalicular efflux")
  reactions <-   addReaction(reactions, "Tca_canalicular", "Tca_buffer", "transport_Tca*Tca_canalicular", "Transport bile")

  mysteadies <- steadyStates(reactions)
  #-! print(mysteadies)

  x <- Xs(odemodel(reactions))

  parameters <- getParameters(x)
  trafo <- `names<-`(parameters, parameters)
  trafo <- repar("inner~steadyEqn", trafo, inner = names(mysteadies), steadyEqn = mysteadies)

  pSS <- P(trafo, condition = "steady")

  set.seed(2)
  pars <- structure(runif( length(getParameters(pSS)), 0,1), names = getParameters(pSS))

  prediction <- (x*pSS)(seq(0,10, 0.1), pars, deriv = F)
  #-! plot(prediction)
  #-!End example code

  is_steady <- function(prediction) {
    mypred <- wide2long(prediction)
    steady_conds <- sapply(split(mypred, mypred[c("name", "condition")]), function(i) {
      zapsmall(var(i$value)) == 0
    })
    return(all(steady_conds))
  }

  unlink("*.c")
  unlink("*.o")
  unlink("*.so")
  unlink("_model.csv")
  
  # Define your expectations here
  expect_length(mysteadies, 3)
  expect_true(is_steady(prediction))

})
