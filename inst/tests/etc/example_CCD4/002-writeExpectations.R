# setwd(rstudioapi::getSourceEditorContext()$path)
source("001-setup.R")
# -------------------------------------------------------------------------#
# 1 Write expectations ----
# -------------------------------------------------------------------------#
# .. -1 create dir -----
try(unlink("expectations", recursive = TRUE))
dir.create("expectations")
# ..0  Objective function  evaluation -----
saveRDS(obj(ini), file.path("expectations","001-objResult.rds"))
# .. 1 normal -----
mydata <- data
resResult <- lapply(setNames(nm = conditions), function(cn) {
  res(mydata[[cn]], prediction[[cn]])
})
saveRDS(mydata, file.path("expectations","011-data.rds"))
saveRDS(resResult, file.path("expectations","012-resResult.rds"))
saveRDS(lapply(resResult,wrss),file.path("expectations","013-wrss.rds"))


# .. 2 error model -----
mydata <- data
resResult <- lapply(setNames(nm = conditions), function(cn) {
  err <- e(prediction[[cn]], ini)
  res(mydata[[cn]], prediction[[cn]], err[[1]])
})
saveRDS(mydata, file.path("expectations","021-data.rds"))
saveRDS(resResult, file.path("expectations","022-resResultErrpars.rds"))
saveRDS(lapply(resResult,nll), file.path("expectations","023-wrss.rds"))

# .. 3 error model, delete sigma first -----
mydata <- as.datalist(lapply(data, function(x) {x$sigma = NA;x}))
resResult <- lapply(setNames(nm = conditions), function(cn) {
  err <- e(prediction[[cn]], ini)
  res(mydata[[cn]], prediction[[cn]], err[[1]])
})
saveRDS(mydata, file.path("expectations", "031-data.rds"))
saveRDS(resResult, file.path("expectations", "032-resResultErrpars.rds"))
saveRDS(lapply(resResult,nll), file.path("expectations", "033-wrss.rds"))

# .. 4 Data with LLOQ, sigma from data -----
mydata <- as.datalist(lapply(data, function(x) {x$lloq = median(x$value);x}))
resResult <- lapply(setNames(nm = conditions), function(cn) {
  # err <- e(prediction[[cn]], ini)
  err <- NULL
  res(mydata[[cn]], prediction[[cn]], err[[1]])
})
saveRDS(mydata, file.path("expectations","041-data.rds"))
saveRDS(resResult, file.path("expectations","042-resResultErrpars.rds"))
saveRDS(lapply(resResult,wrss), file.path("expectations","043-wrss.rds"))

# .. 5 Data with LLOQ, sigma from errorModel -----
mydata <- as.datalist(lapply(data, function(x) {x$lloq = median(x$value);x}))
resResult <- lapply(setNames(nm = conditions), function(cn) {
  err <- e(prediction[[cn]], ini)
  res(mydata[[cn]], prediction[[cn]], err[[1]])
})
saveRDS(mydata, file.path("expectations","051-data.rds"))
saveRDS(resResult, file.path("expectations","052-resResultErrpars.rds"))
saveRDS(lapply(resResult,nll), file.path("expectations","053-wrss.rds"))


# .. 6 Data with LLOQ, delete sigma first -----
mydata <- as.datalist(lapply(data, function(x) {x$lloq = median(x$value);x$sigma = NA;x}))
resResult <- lapply(setNames(nm = conditions), function(cn) {
  err <- e(prediction[[cn]], ini)
  res(mydata[[cn]], prediction[[cn]], err[[1]])
})
saveRDS(mydata, file.path("expectations","061-data.rds"))
saveRDS(resResult, file.path("expectations","062-resResultErrpars.rds"))
saveRDS(lapply(resResult,nll), file.path("expectations","063-wrss.rds"))



# -------------------------------------------------------------------------#
# Clean the folder ----
# -------------------------------------------------------------------------#
for (f in list.files(pattern = "\\.(c|o|so)$")) unlink(f)


