# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> ----
# >>>> TO BE REMOVED <<<<<<<< ----
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> ----
# This file was to perform some manual checks that I didn't break anything while 
#   updating the nll and normL2 function when I ported the stuff back from IQRtools. 
# I just leave it here for completeness at the moment.



# 
# devtools::load_all("~/PROJTOOLS/dMod")
# 
# source("etc/example_CCD4/001-setup.R")
# 
# 
# compare_named_numeric <- function(.x,.y) {
#   .x <- tibble(name = names(.x),x = .x)
#   .y = tibble(name = names(.y), y = .y)
#   out <- merge(.x,.y, all = TRUE)
#   out <-  mutate(out, log10Diff = round(log10(abs(y-x))))
#   out <-  mutate(out, log10Diff = case_when(log10Diff > -20 ~ paste0(log10Diff, "    "), TRUE ~ as.character(log10Diff)))
#   return(out)
# }
# compare_hessians <- function(.x,.y) {
#   .x <- setNames(`dim<-`(.x, NULL), outer(rownames(.x), colnames(.x), paste0))
#   .y <- setNames(`dim<-`(.y, NULL), outer(rownames(.y), colnames(.y), paste0))
#   compare_named_numeric(.x,.y)
# }
# 
# 
# # ..  -----
# # debugonce(dMod:::`+.objlist`)
# # debugonce(obj)
# # debugonce(nll)
# obj(ini)
# 
# # ..  -----
# identical(obj(ini), readRDS(file.path("etc/example_CCD4/expectations","001-objResult.rds")))
# 
# # [x] Values are same
# mapply(.x = obj(ini), .y = readRDS(file.path("etc/example_CCD4/expectations","001-objResult.rds")), FUN = function(.x,.y) identical(.x,.y))
# mapply(.x = obj(ini), .y = readRDS(file.path("etc/example_CCD4/expectations","001-objResult.rds")), FUN = function(.x,.y) {conveniencefunctions::compare(.x,.y)})
# 
# # .. 1 normal -----
# mydata <- data
# resResult <- lapply(setNames(nm = conditions), function(cn) {
#   res(mydata[[cn]], prediction[[cn]])
# })
# # [x] Data is same
# identical(mydata, readRDS(file.path("etc/example_CCD4/expectations","011-data.rds")))
# # [x] res are same: No suffixes .x,.y appended!
# mapply(.x = resResult, .y = readRDS(file.path("etc/example_CCD4/expectations","012-resResult.rds")), FUN = function(.x,.y) {
#   merge(.x,.y)
#   # print("new");print(.x);print("old");print(.y)
# }, SIMPLIFY = F)
# 
# # Values are same
# a <- mapply(.x = lapply(resResult,function(x) nll(x, ini, TRUE)), .y = readRDS(file.path("etc/example_CCD4/expectations","013-wrss.rds")),FUN = function(.x,.y) {
#   # cat("\n\nnew\n");print(.x);cat("\n\nold\n");print(.y)
#   
#   cat(crayon::green("values\n"))
#   print(.x$value)
#   print(.y$value)
#   
#   cat(crayon::green("gradient lengths\n"))
#   print(length(.x$gradient))
#   print(length(.y$gradient))
#   
#   cat(crayon::green("gradients\n"))
#   print(compare_named_numeric(.x$gradient, .y$gradient))
#   
#   
#   cat(crayon::green("hessians\n"))
#   print(compare_hessians(.x$hessian, .y$hessian))
#   a <- NULL
# })
# 
# 
# # .. 2 error model -----
# # expect Difference, because sigma is present in data!
# mydata <- data
# resResult <- lapply(setNames(nm = conditions), function(cn) {
#   err <- e(prediction[[cn]], ini)
#   res(mydata[[cn]], prediction[[cn]], err[[1]])
# })
# # [x] Data is same
# identical(mydata, readRDS(file.path("etc/example_CCD4/expectations","021-data.rds")))
# 
# mapply(.x = resResult, .y = readRDS(file.path("etc/example_CCD4/expectations","022-resResultErrpars.rds")), FUN = function(.x,.y) {
#   merge(.x,.y, all = TRUE)
#   # print("new");print(.x);print("old");print(.y)
# }, SIMPLIFY = F)
# 
# a <- mapply(.x = lapply(resResult,function(x) nll(x, ini, TRUE)), .y = readRDS(file.path("etc/example_CCD4/expectations","023-wrss.rds")),FUN = function(.x,.y) {
#   # cat("\n\nnew\n");print(.x);cat("\n\nold\n");print(.y)
#   
#   cat(crayon::green("values\n"))
#   print(.x$value)
#   print(.y$value)
#   
#   cat(crayon::green("gradient lengths\n"))
#   print(length(.x$gradient))
#   print(length(.y$gradient))
#   
#   cat(crayon::green("gradients\n"))
#   print(compare_named_numeric(.x$gradient, .y$gradient))
#   
#   
#   # cat(crayon::green("hessians\n"))
#   # print(compare_hessians(.x$hessian, .y$hessian))
#   a <- NULL
# })
# 
# 
# 
# 
# 
# # .. 3 error model, delete sigma first -----
# # should be the same again
# mydata <- as.datalist(lapply(data, function(x) {x$sigma = NA;x}))
# resResult <- lapply(setNames(nm = conditions), function(cn) {
#   err <- e(prediction[[cn]], ini)
#   res(mydata[[cn]], prediction[[cn]], err[[1]])
# })
# # [x] Data is same
# identical(mydata, readRDS(file.path("etc/example_CCD4/expectations","031-data.rds")))
# 
# mapply(.x = resResult, .y = readRDS(file.path("etc/example_CCD4/expectations","032-resResultErrpars.rds")), FUN = function(.x,.y) {
#   merge(.x,.y, all = TRUE)
#   # print("new");print(.x);print("old");print(.y)
# }, SIMPLIFY = F)
# 
# a <- mapply(.x = lapply(resResult,function(x) nll(x, ini, TRUE)), .y = readRDS(file.path("etc/example_CCD4/expectations","033-wrss.rds")),FUN = function(.x,.y) {
#   # cat("\n\nnew\n");print(.x);cat("\n\nold\n");print(.y)
#   
#   cat(crayon::green("values\n"))
#   print(.x$value)
#   print(.y$value)
#   
#   cat(crayon::green("gradient lengths\n"))
#   print(length(.x$gradient))
#   print(length(.y$gradient))
#   
#   cat(crayon::green("gradients\n"))
#   print(compare_named_numeric(.x$gradient, .y$gradient))
#   
#   
#   # cat(crayon::green("hessians\n"))
#   # print(compare_hessians(.x$hessian, .y$hessian))
#   a <- NULL
# })
# 
# 
# # .. 4 Data with LLOQ, sigma from data -----
# mydata <- as.datalist(lapply(data, function(x) {x$lloq = median(x$value);x}))
# resResult <- lapply(setNames(nm = conditions), function(cn) {
#   # err <- e(prediction[[cn]], ini)
#   err <- NULL
#   res(mydata[[cn]], prediction[[cn]], err[[1]])
# })
# # [x] Data is same
# identical(mydata, readRDS(file.path("etc/example_CCD4/expectations","041-data.rds")))
# 
# mapply(.x = resResult, .y = readRDS(file.path("etc/example_CCD4/expectations","042-resResultErrpars.rds")), FUN = function(.x,.y) {
#   merge(.x,.y, all = TRUE)
#   # print("new");print(.x);print("old");print(.y)
# }, SIMPLIFY = F)
# 
# a <- mapply(.x = lapply(resResult,function(x) nll(x, ini, TRUE)), .y = readRDS(file.path("etc/example_CCD4/expectations","043-wrss.rds")),FUN = function(.x,.y) {
#   # cat("\n\nnew\n");print(.x);cat("\n\nold\n");print(.y)
#   
#   cat(crayon::green("values\n"))
#   print(.x$value)
#   print(.y$value)
#   
#   cat(crayon::green("gradient lengths\n"))
#   print(length(.x$gradient))
#   print(length(.y$gradient))
#   
#   cat(crayon::green("gradients\n"))
#   print(compare_named_numeric(.x$gradient, .y$gradient))
#   
#   
#   # cat(crayon::green("hessians\n"))
#   # print(compare_hessians(.x$hessian, .y$hessian))
#   a <- NULL
# })
# 
# # .. 5 Data with LLOQ, sigma from errorModel -----
# # expect difference
# mydata <- as.datalist(lapply(data, function(x) {x$lloq = median(x$value);x}))
# resResult <- lapply(setNames(nm = conditions), function(cn) {
#   err <- e(prediction[[cn]], ini)
#   res(mydata[[cn]], prediction[[cn]], err[[1]])
# })
# # [x] Data is same
# identical(mydata, readRDS(file.path("etc/example_CCD4/expectations","051-data.rds")))
# 
# mapply(.x = resResult, .y = readRDS(file.path("etc/example_CCD4/expectations","052-resResultErrpars.rds")), FUN = function(.x,.y) {
#   merge(.x,.y, all = TRUE)
#   # print("new");print(.x);print("old");print(.y)
# }, SIMPLIFY = F)
# 
# a <- mapply(.x = lapply(resResult,function(x) nll(x, ini, TRUE)), .y = readRDS(file.path("etc/example_CCD4/expectations","053-wrss.rds")),FUN = function(.x,.y) {
#   # cat("\n\nnew\n");print(.x);cat("\n\nold\n");print(.y)
#   
#   cat(crayon::green("values\n"))
#   print(.x$value)
#   print(.y$value)
#   
#   cat(crayon::green("gradient lengths\n"))
#   print(length(.x$gradient))
#   print(length(.y$gradient))
#   
#   cat(crayon::green("gradients\n"))
#   print(compare_named_numeric(.x$gradient, .y$gradient))
#   
#   
#   # cat(crayon::green("hessians\n"))
#   # print(compare_hessians(.x$hessian, .y$hessian))
#   a <- NULL
# })
# 
# 
# # .. 6 Data with LLOQ, delete sigma first -----
# mydata <- as.datalist(lapply(data, function(x) {x$lloq = median(x$value);x$sigma = NA;x}))
# resResult <- lapply(setNames(nm = conditions), function(cn) {
#   err <- e(prediction[[cn]], ini)
#   res(mydata[[cn]], prediction[[cn]], err[[1]])
# })
# # [x] Data is same
# identical(mydata, readRDS(file.path("etc/example_CCD4/expectations","061-data.rds")))
# 
# mapply(.x = resResult, .y = readRDS(file.path("etc/example_CCD4/expectations","062-resResultErrpars.rds")), FUN = function(.x,.y) {
#   merge(.x,.y, all = TRUE)
#   # print("new");print(.x);print("old");print(.y)
# }, SIMPLIFY = F)
# 
# a <- mapply(.x = lapply(resResult,function(x) nll(x, ini, TRUE)), .y = readRDS(file.path("etc/example_CCD4/expectations","063-wrss.rds")),FUN = function(.x,.y) {
#   # cat("\n\nnew\n");print(.x);cat("\n\nold\n");print(.y)
#   
#   cat(crayon::green("values\n"))
#   print(.x$value)
#   print(.y$value)
#   
#   cat(crayon::green("gradient lengths\n"))
#   print(length(.x$gradient))
#   print(length(.y$gradient))
#   
#   cat(crayon::green("gradients\n"))
#   print(compare_named_numeric(.x$gradient, .y$gradient))
#   
#   
#   # cat(crayon::green("hessians\n"))
#   # print(compare_hessians(.x$hessian, .y$hessian))
#   a <- NULL
# })
# 
# 
# # -------------------------------------------------------------------------#
# # SectionTitle ----
# # -------------------------------------------------------------------------#
# for (f in list.files(pattern = "\\.(c|o|so)$")) unlink(f)
# #-!End example code
# 
# 
# # Define your expectations here
# 
# #})