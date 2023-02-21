library(tidyverse)
library(data.table)
library(dMod)


# 1. Data -----------------------------------------------------------------
data <- as.datalist(
  as.data.frame(
    tibble::tribble(
      ~time,  ~name,     ~value,  ~condition,
      0.00000, "B_obs", 1.0430899, "reference",
      10.34483, "B_obs", 1.8257833, "reference",
      20.68966, "B_obs", 1.8324904, "reference",
      31.03448, "B_obs", 1.7451547, "reference",
      41.37931, "B_obs", 1.7168632, "reference",
      51.72414, "B_obs", 1.7008545, "reference",
      62.06897, "B_obs", 1.5732636, "reference",
      72.41379, "B_obs", 1.4172603, "reference",
      82.75862, "B_obs", 1.4186770, "reference",
      93.10345, "B_obs", 1.4503594, "reference",
      103.44828, "B_obs", 1.4354153, "reference",
      113.79310, "B_obs", 1.4107170, "reference",
      124.13793, "B_obs", 1.3334672, "reference",
      134.48276, "B_obs", 1.2615258, "reference",
      144.82759, "B_obs", 1.2544703, "reference",
      155.17241, "B_obs", 1.2297995, "reference",
      165.51724, "B_obs", 1.1575647, "reference",
      175.86207, "B_obs", 1.2684661, "reference",
      186.20690, "B_obs", 1.1534651, "reference",
      196.55172, "B_obs", 1.1687653, "reference",
      206.89655, "B_obs", 1.1864222, "reference",
      217.24138, "B_obs", 1.1649452, "reference",
      227.58621, "B_obs", 1.0603839, "reference",
      237.93103, "B_obs", 1.1279193, "reference",
      248.27586, "B_obs", 1.0641429, "reference",
      258.62069, "B_obs", 0.9844280, "reference",
      268.96552, "B_obs", 1.0220764, "reference",
      279.31034, "B_obs", 1.0393746, "reference",
      289.65517, "B_obs", 1.0199705, "reference",
      300.00000, "B_obs", 0.9994695, "reference",
      0.00000, "B_obs", 1.0046959, "treatment",
      10.34483, "B_obs", 1.6289925, "treatment",
      20.68966, "B_obs", 1.5229543, "treatment",
      31.03448, "B_obs", 1.4797702, "treatment",
      41.37931, "B_obs", 1.3145159, "treatment",
      51.72414, "B_obs", 1.2620031, "treatment",
      62.06897, "B_obs", 1.1857014, "treatment",
      72.41379, "B_obs", 1.0865067, "treatment",
      82.75862, "B_obs", 1.0479692, "treatment",
      93.10345, "B_obs", 1.0255577, "treatment",
      103.44828, "B_obs", 1.0552660, "treatment",
      113.79310, "B_obs", 1.0278639, "treatment",
      124.13793, "B_obs", 1.1006438, "treatment",
      134.48276, "B_obs", 1.0550504, "treatment",
      144.82759, "B_obs", 1.0271107, "treatment",
      155.17241, "B_obs", 1.0205512, "treatment",
      165.51724, "B_obs", 0.9966564, "treatment",
      175.86207, "B_obs", 0.9092575, "treatment",
      186.20690, "B_obs", 0.9252080, "treatment",
      196.55172, "B_obs", 0.9590331, "treatment",
      206.89655, "B_obs", 0.9861776, "treatment",
      217.24138, "B_obs", 0.9451449, "treatment",
      227.58621, "B_obs", 1.0572657, "treatment",
      237.93103, "B_obs", 0.9400520, "treatment",
      248.27586, "B_obs", 1.0110617, "treatment",
      258.62069, "B_obs", 0.8694348, "treatment",
      268.96552, "B_obs", 0.9706745, "treatment",
      279.31034, "B_obs", 0.9269390, "treatment",
      289.65517, "B_obs", 0.9725193, "treatment",
      300.00000, "B_obs", 0.9615111, "treatment"
    )
  ),
  split.by = "condition"
)




# 2. Equation list --------------------------------------------------------
el <- NULL
el <- addReaction(el, from = "A", to = "B", rate = "k_AB * A", description = "A to B")
el <- addReaction(el, from = "B", to = "C", rate = "k_BC * B", description = "B to C")


# 3. Observables ----------------------------------------------------------
observables <- eqnvec(
  B_obs = "scale_B * B + offset_B"
)


# 4. Error model ----------------------------------------------------------
errors <- c(B_obs = "sigma_B_obs")


# 5. Create ODE model -----------------------------------------------------
model <- odemodel(
  f = el
)
x <- Xs(
  odemodel = model
)


# 6. Observation function -------------------------------------------------
g <- Y(
  g = observables,
  f = el
)


# 7. Compile error model --------------------------------------------------
e <- Y(
  g = errors,
  f = c(as.eqnvec(el), observables),
  states = names(observables) 
)


# 8. Transformation function ----------------------------------------------
paramInner <- c(unique(c(getParameters(el), getSymbols(observables), getSymbols(errors))))
names(paramInner) <- paramInner

trafo <- as.eqnvec(paramInner)
conditionGrid <- attr(data, "condition.grid")
trafo <- insert(trafo, "x~y", x = c("B", "C"), y = 0)
trafoList <- branch(trafo, conditionGrid)

conditions <- names(trafoList)

for(i in seq_along(conditions)) {
  trafoList[[i]] <- insert(trafoList[[i]], paste0("x~x_",conditions[i]) , x = "k_BC")
}
trafoList <- insert(trafoList, "x~10^(x)", x = .currentSymbols)

p <- P(
  trafo = trafoList
)


# 9. Compile prediction function ------------------------------------------
prdFunc <- Reduce("*", list(g, x, p))


# 10. Define Objective Function -------------------------------------------
paramOuter <- structure(rep(-1, length(getParameters(p))), names = getParameters(p))

obj <- normL2(
  data = data,
  x = prdFunc,
  errmodel = e
) + constraintL2(paramOuter, sigma = 10)


# 11. Fit -----------------------------------------------------------------
fitOutput <- mstrust(
  objfun=obj, 
  center=dMod::msParframe(paramOuter, n = 20, seed=1),#pouter,#
  studyname="fits", 
  fits = 20, 
  cores = 1
)


# 12. View fit results ----------------------------------------------------
fitResults  <- as.parframe(fitOutput)
bestFit  <- fitResults[1]


# 12.1 Waterfall Plot -----------------------------------------------------
plotValues(fitResults)


# 12.2 Model predictions together with measured data ----------------------
fitPlot <- plotCombined(prdFunc(seq(300), bestFit), data)
ggsave(
  plot = fitPlot,
  filename = "./fitPlot.pdf",
       width = 200,
       height = 150,
       units = "mm",  
  device = "pdf",         
  title = ""
)


# 13. Profiles ------------------------------------------------------------

profiles <- profile(
  obj = obj,
  pars =  bestFit,
  whichPar = names(paramOuter),
  cores = 1,
  method = "optimize",
  optControl = list(iterlim = 50)
)

profilePlot <- plotProfile(profiles, mode == "data")
ggsave(
  plot = profilePlot,
  filename = "./profilePlot.pdf",
  width = 300,
  height = 200,
  units = "mm",  
  device = "pdf",         
  title = ""
)
