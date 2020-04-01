
\dontrun{
library(dMod)
# library(dplyr) # devtools::install_github("dlill/conveniencefunctions")
library(conveniencefunctions)
  
setwd(tempdir())  

## Model definition (text-based, scripting part)
reactions <- NULL %>%
  addReaction("A", "B", "k1*A", "translation") %>%
  addReaction("B",  "", "k2*B", "degradation") 

f <- reactions %>%
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


myframe <- dMod.frame("no steady states", g, x, p, data) %>%
  appendObj() %>%
  mutate(constr = list(constraintL2(mu = 0*pars, sigma = 5)),
         obj = list(obj_data + constr),
         fits = list(mstrust(obj, pars, studyname = "Fits", fits = 20, cores = 4, blather = T))) %>%
  appendParframes() %>%
  mutate(profiles = list(profile(obj, as.parvec(parframes), whichPar = "k1"))) %>%
  mutate(vali = list(datapointL2("A", 2, "mypoint", .1, condition = "a")),
         obj_vali = list(obj_data + constr + vali),
         par_vali = list(c(dMod:::sanitizePars(as.parvec(parframes))$pars, "mypoint" = 0.1 )),
         fits_vali = list(mstrust(obj_vali, par_vali)),
         profile_vali = list(profile(obj_vali, fits_vali %>% as.parframe %>% as.parvec, "mypoint"))) %>% 
  {.}

myframe <- myframe %>% tibble::add_column(reactions = list(reactions), fixed = list(NULL))


saveShiny_dMod.frame(myframe, projectname = "dModFrameTest")

}
