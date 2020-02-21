# context("nll_wrss")
# test_that("IQR-nll", {

#-!Start example code

devtools::load_all("~/PROJTOOLS/dMod")

source("etc/example_CCD4/001-setup.R")
# ..  -----
# debugonce(dMod:::`+.objlist`)
# debugonce(obj)
# debugonce(nll)
obj(ini)

# ..  -----
identical(obj(ini), readRDS(file.path("expectations","001-objResult.rds")))
# .. 1 normal -----
mydata <- data
resResult <- lapply(setNames(nm = conditions), function(cn) {
  res(mydata[[cn]], prediction[[cn]])
})
identical(mydata, readRDS(file.path("expectations","011-data.rds")))
identical(resResult, readRDS(file.path("expectations","012-resResult.rds")))
identical(lapply(resResult,wrss), readRDS(file.path("expectations","013-wrss.rds")))


# .. 2 error model -----
mydata <- data
resResult <- lapply(setNames(nm = conditions), function(cn) {
  err <- e(prediction[[cn]], ini)
  res(mydata[[cn]], prediction[[cn]], err[[1]])
})
identical(mydata, readRDS(file.path("expectations","021-data.rds")))
identical(resResult, readRDS(file.path("expectations","022-resResultErrpars.rds")))
identical(lapply(resResult,nll), readRDS(file.path("expectations","023-wrss.rds")))

# .. 3 error model, delete sigma first -----
mydata <- as.datalist(lapply(data, function(x) {x$sigma = NA;x}))
resResult <- lapply(setNames(nm = conditions), function(cn) {
  err <- e(prediction[[cn]], ini)
  res(mydata[[cn]], prediction[[cn]], err[[1]])
})
identical(mydata, readRDS(file.path("expectations", "031-data.rds")))
identical(resResult, readRDS(file.path("expectations", "032-resResultErrpars.rds")))
identical(lapply(resResult,nll), readRDS(file.path("expectations", "033-wrss.rds")))

# .. 4 Data with LLOQ, sigma from data -----
mydata <- as.datalist(lapply(data, function(x) {x$lloq = median(x$value);x}))
resResult <- lapply(setNames(nm = conditions), function(cn) {
  # err <- e(prediction[[cn]], ini)
  err <- NULL
  res(mydata[[cn]], prediction[[cn]], err[[1]])
})
identical(mydata, readRDS(file.path("expectations","041-data.rds")))
identical(resResult, readRDS(file.path("expectations","042-resResultErrpars.rds")))
identical(lapply(resResult,wrss), readRDS(file.path("expectations","043-wrss.rds")))

# .. 5 Data with LLOQ, sigma from errorModel -----
mydata <- as.datalist(lapply(data, function(x) {x$lloq = median(x$value);x}))
resResult <- lapply(setNames(nm = conditions), function(cn) {
  err <- e(prediction[[cn]], ini)
  res(mydata[[cn]], prediction[[cn]], err[[1]])
})
identical(mydata, readRDS(file.path("expectations","051-data.rds")))
identical(resResult, readRDS(file.path("expectations","052-resResultErrpars.rds")))
identical(lapply(resResult,nll), readRDS(file.path("expectations","053-wrss.rds")))


# .. 6 Data with LLOQ, delete sigma first -----
mydata <- as.datalist(lapply(data, function(x) {x$lloq = median(x$value);x$sigma = NA;x}))
resResult <- lapply(setNames(nm = conditions), function(cn) {
  err <- e(prediction[[cn]], ini)
  res(mydata[[cn]], prediction[[cn]], err[[1]])
})
identical(mydata, readRDS(file.path("expectations","061-data.rds")))
identical(resResult, readRDS(file.path("expectations","062-resResultErrpars.rds")))
identical(lapply(resResult,nll), readRDS(file.path("expectations","063-wrss.rds")))


for (f in list.files(pattern = "\\.(c|o|so)$")) unlink(f)
#-!End example code


# Define your expectations here

#})