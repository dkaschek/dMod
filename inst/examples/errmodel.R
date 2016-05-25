
library(dMod)
library(deSolve)

f <- NULL
f <- addReaction(f, "A", "B", "k1*A", "Production of B")
f <- addReaction(f, "B", "C", "k2*B", "Production of C")

model <- odemodel(f)
p <- P(
  trafo = eqnvec(
    A = "1", B = "0", C = "0", k1 = "exp(k1)", k2 = "exp(k2)", # dynamic parameters
    off_B = "exp(off_B)", # observation parameters
    sigma_rel = "exp(sigma_rel)" # error parameters
  ),
  condition = "C1"
)
x <- Xs(model)
g <- Y(eqnvec(B_obs = "B + off_B"), f)
err <- Y(eqnvec(B_obs = "sigma_rel*(B + off_B)"), f, attach.input = FALSE)

## Simulate data
ptrue <- c(k1 = -1, k2 = -2, off_B = -1, sigma_rel = log(0.1))
times <- 0:20
prediction <- (g*x*p)(times, ptrue)
datasheet <- subset(wide2long(prediction), name == "B_obs")
datasheet$sigma <- 0.1*datasheet$value
datasheet$value <- datasheet$value + rnorm(length(datasheet$value), 0, datasheet$sigma)
data <- as.datalist(datasheet)

## Fit data with error model
myfit <- trust(normL2(data, g*x*p, err), ptrue + rnorm(length(ptrue), 0, 1), rinit = 1, rmax = 10)

## Plotting
out <- ggdata_fn(g*x*p, err, data, times = seq(0, 20, len = 100), pars = myfit$argument)
ggplot(out$prediction, aes(x = time, y = value, ymin = value-sigma, ymax = value+sigma, 
                           group = condition, color = condition, fill = condition)) +
  facet_wrap(~name, scales = "free") +
  geom_line() + geom_ribbon(alpha = .2, lty = 0) + 
  geom_point(data = out$data) + geom_errorbar(data = out$data, width = 0) +
  theme_dMod() + scale_color_dMod() + scale_fill_dMod()
