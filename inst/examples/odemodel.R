\dontrun{

## Generate a compiled ODE model from an equation vector
## The model will not return sensitivities for "switch"
## Files will be generated in your working directory!

f <- eqnvec(A = "-k*A + switch*F")
model <- odemodel(f, forcings = "F", fixed = "switch")
print(model)

## Generate the same model from an equation list
f <- addReaction(NULL, from = "", to = "A", rate = "switch*F", description = "production")
f <- addReaction(f   , from = "A", to = "", rate = "k*A", description = "degradation")
print(f)

model <- odemodel(f, forcings = "F", fixed = "switch")
print(model)


# create forcings
forc1 <- data.frame(name = "F", time = seq(0,5, 0.1), value = sin(seq(0,5,0.1)))
forc2 <- data.frame(name = "F", time = seq(0,5, 0.1), value = exp(-seq(0,5,0.1)))
forc3 <- data.frame(name = "F", time= 0,              value = 0.1)


x <- Xs(model, forc1, condition = "forc1") + 
  Xs(model, forc2, condition = "forc2") + 
  Xs(model, forc3, condition = "forc3")

g <- Y(c(out1 = "F * A", out2 = "F"), x)

times <-  seq(0,5, 0.001)
pars <- setNames(runif(length(getParameters(x))), getParameters(x))

pred <- (g*x)(times, pars)  
plot(pred)


}
