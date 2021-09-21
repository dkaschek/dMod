# # -------------------------------------------------------------------------#
# # Testing ----
# # -------------------------------------------------------------------------#
try(setwd(dirname(rstudioapi::getSourceEditorContext()$path)))
setwd("..")
# devtools::load_all()
library(dMod)
f <- list.files("BenchmarkModels")

petab <- importPEtabSBML_indiv(modelname = f[1],
                         path2model = "BenchmarkModels/",
                         testCases = FALSE,
                         path2TestCases = "PEtabTests/",
                         compile = TRUE,
                         SBML_file = NULL,
                         observable_file = NULL,
                         condition_file = NULL,
                         data_file = NULL,
                         parameter_file = NULL)

p <- petab$fns$p0
x <- petab$fns$x
times <- seq(0,max(as.data.frame(petab$data)$time), len=501)
pred <- petab$prd(times, petab$pars, FLAGbrowserN = 1)
plotCombined(pred, petab$data)

# # -------------------------------------------------------------------------#
# # Test models ----
# # -------------------------------------------------------------------------#
# try(setwd(dirname(rstudioapi::getSourceEditorContext()$path)))
# setwd("..")
# devtools::load_all()
# f <- list.files("PEtabTests/")
# i <- 1
# # ..  -----
# # debugonce(importPEtabSBML_indiv)
# petab <- importPEtabSBML_indiv(modelname = f[i],
#                          path2model = "BenchmarkModels/",
#                          testCases = TRUE,
#                          path2TestCases = "PEtabTests/",
#                          compile = TRUE,
#                          SBML_file = NULL,
#                          observable_file = NULL,
#                          condition_file = NULL,
#                          data_file = NULL,
#                          parameter_file = NULL)
# 
# p <- petab$fns$p0
# x <- petab$fns$x
# times <- seq(0,max(as.data.frame(petab$data)$time), len=501)
# pred <- petab$prd(times, petab$pars, FLAGbrowserN = 1)
# plotCombined(pred, petab$data)
# i <- i+1
# 
# # -------------------------------------------------------------------------#
# #  ----
# # -------------------------------------------------------------------------#
# 
# 
# # prd(times, myfit_values, FLAGbrowser = 1)
# # prd(times, myfit_values, FLAGbrowser = 2)
# 
# # p <- P_indiv(myp, est.grid = gridlist$est.grid, fix.grid = gridlist$fix.grid)
# # wup <- p(myfit_values)
# # wup
# # myfit_values
# 
rp <- tempfile()
Rprof(rp)
petab$obj_data(petab$pars)
Rprof(NULL)
summaryRprof(rp)
pv <- profvis::profvis(prof_input = rp); htmlwidgets::saveWidget(pv, paste0(rp, ".html")); browseURL(paste0(rp, ".html"))
# 
# # debugonce(obj_data)
# lapply(1:10,function(i)petab$obj_data(petab$pars))
# 
# parallel::mclapply(1:12, function(i) obj_data(myfit_values), mc.cores = 4)
# 
# 
# 
# # obj_data(myfit_values, FLAGbrowser = T)
# 
b1 <- rbenchmark::benchmark(petab$obj_data(petab$pars), replications = 20)

b1.2 <- rbenchmark::benchmark(mclapply(1:36, function(i) petab$obj_data(petab$pars), 
                                       mc.cores= 12, mc.preschedule = TRUE),
                              replications = 3)



importPEtabSBML(f[1])
# debugonce(obj)
# obj(pouter)
#
# rp <- tempfile()
# Rprof(rp)
# obj(pouter)
# Rprof(NULL)
# summaryRprof(rp)
# pv <- profvis::profvis(prof_input = rp); htmlwidgets::saveWidget(pv, paste0(rp, ".html")); browseURL(paste0(rp, ".html"))

b2 <- rbenchmark::benchmark(obj(pouter), replications = 20)
b2.2 <- rbenchmark::benchmark(mclapply(1:36, function(i) obj(pouter), mc.cores= 12),
                              replications = 3)

writeLines(capture.output(print(list(
  b1,
  b1.2,
  b2,
  b2.2
))), "~/wup.txt")


# # Exit ----
