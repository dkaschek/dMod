context("Print parameter transformation")


# Set up environment
# Construct parameter transformation from topology and catch print.
reactionlist <- read.csv("dataSets/topo.csv")
f <- generateEquations(reactionlist)
observables <- c(obs_sufRR = "log(scl_Fe * Fe + off_Fe)")
forcings <- NULL
innerpars <- getSymbols(c(f, names(f), observables), exclude = c(forcings, "time"))
names(innerpars) <- innerpars
constraints <- resolveRecurrence(c())
trafo <- replaceSymbols(names(observables), observables, innerpars)
trafo <- replaceSymbols(names(constraints), constraints, trafo)
trafo <- replaceSymbols(innerpars[!grepl("^s_", innerpars)], paste0("exp(log_", innerpars[!grepl("^s_", innerpars)], ")"), trafo)
pT <- P(trafo)
pTOutput <- capture.output(pT)

# Load print reference
refFile <- file("dataSets/topo-para_trans_reference.txt")
reference <- readLines(refFile)



test_that("parameters are printed works", {
  expect_that(pTOutput, equals(reference))
})