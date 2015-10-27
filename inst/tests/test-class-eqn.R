context("Class -- eqnvec")

# Tested functions
# - eqnvec
# - as.eqnvec
# - is.eqnvec
#
# Missing functions
# - getReactions
# - addReaction.eqnlist
# - getFluxes
# - as.data.frame.eqnlist
# - write.eqnlist
# - subset.eqnlist
# - print.eqnlist
# - pander.eqnlist
# - format.eqnvec
# - as.eqnvec.eqnlist
# - format.eqnvec
# - print.eqnvec
# - pander.eqnvec
# - dot.eqnvec


# Tests
test_that("creation and validation of equations works", {
  
  m_eqn <- c(a = "3*x+y", b = "2*z")
  m_names <- names(m_eqn)
  m_eqnvec <- eqnvec(equations = m_eqn, names = m_names)
  m_eqnvec2 <- eqnvec(equations = m_eqn, names = c("y", "z"))
  
  expect_true(is.null(eqnvec()))
  expect_true(is.eqnvec(m_eqnvec))
  expect_equal(names(m_eqnvec), m_names)
  expect_equivalent(m_eqnvec[m_names[1]], m_eqn[1])
  expect_equivalent(m_eqnvec[m_names[2]], m_eqn[2])
})

test_that("coerce different data types to eqnvec works", {
  
  m_eqn <- c(a = "3*x+y", b = "2*z")
  m_names <- names(m_eqn)
  m_eqnvec <- eqnvec(equations = m_eqn, names = m_names)
  
  expect_equal(as.eqnvec(m_eqn), m_eqnvec)
})

test_that("concatenate and other operations work", {
  
  m_eqn <- c(a = "3*x+y", b = "2*z")
  m_names <- names(m_eqn)
  m_eqnvec <- eqnvec(equations = m_eqn, names = m_names)
  m_eqnvec2 <- eqnvec(equations = m_eqn, names = c("y", "z"))
  
  expect_error(c(m_eqnvec, m_eqnvec), "Names must be unique")
  expect_true(is.eqnvec(c(m_eqnvec, m_eqnvec2)))
})



context("Class -- eqnlist")

# Tested functions
# - eqnlist

test_that("creation and validation of equation lists works", {

  m_topo <- read.csv("dataSets/topoIsaR1.csv")
  m_sm <- m_topo[, -(1:2)]
  m_rates <- m_topo[, 2]
  m_description <- m_topo[, 1]
  m_eqnlist <- eqnlist(smatrix = m_sm, rates = m_rates)
  m_eqnlistDsc <- eqnlist(smatrix = m_sm, rates = m_rates, description = m_description)
  m_eqnlistCSV <- as.eqnlist(m_topo)
  
  expect_true(is.eqnlist(eqnlist()))
  expect_true(is.eqnlist(m_eqnlist))
  expect_true(is.eqnlist(m_eqnlistDsc))
  expect_true(is.eqnlist(m_eqnlistCSV))
  expect_equal(m_eqnlistCSV, m_eqnlistDsc)
  expect_equal(m_eqnlist$description, as.character(1:length(m_eqnlist$description)))
  expect_equal(m_eqnlistDsc$description, as.character(m_description))
})



