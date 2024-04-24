#Scale dependent flux sensitivity Analysis for the SIR-Model with two different values for beta

frame <- NULL
frame <-addReaction(frame, from = "S", to = "I", rate = "beta*S*I")
frame <-addReaction(frame, from = "I", to = "R", rate = "gamma*I")
x <-Xs(odemodel(frame))

obstates = c("S", "I", "R")

par <- getParameters(x)
p <- NULL
for (i in c(0.0026, 0.0027)){
  trafo <- repar("x~x", x=par)
  trafo <- repar("x~y", x = "beta", y= i, trafo)
  parameter <- structure(trafo, names = par)
  p <- p+P(parameter, condition= paste0(i))}


time <- seq(0, 14, length.out=1000)
bestfit <- c(S = 7.637603e+02, I = 6.184001e-01, R = 8.251977e-16, gamma = 4.582934e-01)

scaledependent(frame, x, obstates, bestfit, time, p)