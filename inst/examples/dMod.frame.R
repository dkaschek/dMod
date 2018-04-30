\dontrun{

# ---- Example
library(dMod)
library(conveniencefunctions)

## Model definition (text-based, scripting part)
f <- NULL %>%
  addReaction("A", "B", "k1*A", "translation") %>%
  addReaction("B",  "", "k2*B", "degradation") %>%
  as.eqnvec()


events <- eventlist(var = "A", time = 5, value = "A_add", method = "add")

x <- odemodel(f, events = events) %>% Xs

g <- Y(c(Bobs = "s1*B"), x, compile = T, modelname = "obsfn")

conditions <- c("a", "b")

# there is a bug in
# getParameters(g*x)
parameters <- union(getParameters(g), getParameters(x))

trafo <-
  NULL %>%
  define("x~x", x = parameters) %>%
  branch(conditions = conditions) %>%
  insert("x~x_cond", x = "s1", cond = condition) %>%
  insert("x~1", x = "added", conditionMatch = "a") %>%
  insert("x~5", x = "added", conditionMatch = "b") %>%
  insert("x~exp(x)", x = getSymbols(mytrafo[[i]])) %>%
  {.}

p <- P(trafo)

# Parameter transformation for steady states
pSS <- P(f, method = "implicit", condition = NULL, compile = T)




## Process data

# Data
data <- datalist(
  a = data.frame(time = c(0, 2, 7),
                 name = "Bobs",
                 value = c(.3, .3, .3),
                 sigma = c(.03, .03, .03)),
  b = data.frame(time = c(0, 2, 7),
                 name = "Bobs",
                 value = c(.1, .1, .2),
                 sigma = c(.01, .01, .02))
)


# construct dMod.frame
myframe1 <- dMod.frame("no steady states", g, x, p, data)
print(myframe1)


# Augment by derived objects: prd, obj_data, obj, times, pars
set.seed(4)
myframe2 <- appendObj(myframe1)
print(myframe2)


# Plot the model with random pars
plotCombined(myframe2)

# Fit with prior
myframe3 <- myframe2 %>%
  mutate(constr = list(constraintL2(mu = 0*pars, sigma = 5)),
         obj = list(obj_data + constr))

myframe4 <- myframe3 %>%
  mutate(fits = list(mstrust(obj, pars, studyname = "Fits", fits = 20, cores = 4, blather = T)))


# Inspect Fits
myframe5 <- myframe4 %>%
  appendParframes

# Visualize
# Little Bug: If you want to use the dots that go to subset(), you need to specify the others
plotCombined(myframe5, 1, 1, str_detect(name, "B"))

plotPars(myframe5)

plotValues(myframe5)

plotValues(myframe5, 1, tol =0.0000000001, value >1)



# Profiles
myframe6 <- myframe5 %>%
  mutate(profiles = list(profile(obj, as.parvec(parframes), whichPar = "k1")))

myframe6$profiles %>% plotProfile()


# Validation

myframe7 <- myframe6 %>%
  mutate(vali = list(datapointL2("A", 2, "mypoint", .1, condition = "a")),
         obj_vali = list(obj_data + constr + vali),
         par_vali = list(c(dMod:::sanitizePars(as.parvec(parframes))$pars, "mypoint" = 0.1 )),
         fits_vali = list(mstrust(obj_vali, par_vali)),
         profile_vali = list(profile(obj_vali, fits_vali %>% as.parframe %>% as.parvec, "mypoint")))

myframe7$profile_vali %>% plotProfile()


# stage for commit
git_add_dMod.frame(myframe7)



}
