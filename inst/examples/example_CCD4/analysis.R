# publised in https://www.ncbi.nlm.nih.gov/pubmed/27811075
# D2D: model named Bruno
# simple model: sort of A -> B -> C, conditions differ by initial values of A, B, C
# bcar/betacar: symmetric carotene build of three pieces bio-middle-bio
# 1. cut off one end: bcar -> b10 + bio where b10=bio-middle
# 2. cut off second end: b10 -> bio + middle
# -> measure bcar and b10
# bcry/betacry: one modified end ohbio-middle-obio
# zea: both ends modified ohbio-middle-ohbio
# -> modified rates for cutting modified ends
# no error model fitting, observation function without scalings and offsets

source("inst/examples/example_CCD4/setup.R")
do.profs <- FALSE

## Fitting  -------------------------------------------------

# Data times
timesD <- unique(sort(c(unlist(sapply(data, function(d) c(0, d$time))))))
times <- sort(union(timesD,seq(min(timesD), max(timesD), len=101)))

# Fitting
ini <- pouter # random start
cat("test obj with random ini:",obj(ini)$value, "\n")

ini <- as.parvec(readRDS("inst/examples/example_CCD4/mstrust_CCD4.Rds")) 
cat("test obj with published values:",obj(ini)$value, "\n")

center <- ini
print("start fit")
myfit <- trust(obj, ini, rinit=1, rmax=10, iterlim=500,blather=TRUE, conditions = conditions)
print(paste("fit",myfit$converged,"with",myfit$iterations,"iterations to",myfit$value))

# Prediction  
prediction <- (g*x*p)(c(times, seq(130,210, length.out = 10)), myfit$argument, conditions = conditions)

# Plots  
print(plotCombined(prediction,  data) + theme_bw() + xlab("time[min]") + ylab("concentration [nM]"))

paper_names <- c(obcar= "β-carotene", obcry = "β-cryptoxanthin", obio= "β-io", ob10 = "β-10", oohb10 = "3OH β-10", ozea="β-zea")
prediction <- as.data.frame(prediction, data = data)

P <- ggplot(prediction, aes(x = time, y = value, group = name, color = name, ymin = value - sigma, ymax = value + sigma)) +
  geom_line() +
  facet_wrap(~condition) +
  theme_dMod() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) + 
  scale_color_dMod(labels = paper_names) + 
  geom_point(data = as.data.frame(data)) + geom_errorbar(data = as.data.frame(data), width = 0, linetype = 1) + 
  xlab("") + ylab("") 

print(P)

# Profile likelihood
if(do.profs){
  cat("now do profs \n")
  #profiles.exact  <- do.call(rbind,mclapply(names(myfit$argument), function(n) profile(obj, myfit$argument, n, method = "integrate", limits = c(-5, 5), alpha = 1e-6, stepControl = list(min = .1), conditions = conditions),mc.cores=3))
  profiles.exact <- readRDS("inst/examples/example_CCD4/profiles.Rds")
  # plot of all profiles on log-scale
  print(plotProfile(profiles.exact))
  
  # want to plot the profiles on the non-log scale
  npp <- attr(profiles.exact, "parameters")
  snpp <- substr(npp, 4, nchar(npp))
  profiles.exact[,npp] <- exp(profiles.exact[,npp])
  colnames(profiles.exact) <- c(colnames(profiles.exact)[1:6],snpp)
  attr(profiles.exact, "parameters") <- snpp
  profiles.exact$whichPar <- substr(profiles.exact$whichPar, 4, nchar(profiles.exact$whichPar))
  
  # plot the rate constants on non-log scale
  print(plotProfile(profiles.exact, maxvalue = 5, grepl("k", name) & !(grepl("k_zea", name))  & mode == "data") + facet_wrap(~name, ncol=1))
}
