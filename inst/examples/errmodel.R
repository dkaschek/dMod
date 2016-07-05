
library(dMod)
library(deSolve)

f <- NULL
f <- addReaction(f, "A", "B", "k1*A", "Production of B")
f <- addReaction(f, "B", "C", "k2*B", "Production of C")

observables <- eqnvec(B_obs = "B + off_B")
errors <- eqnvec(B_obs = "sigma_rel*B_obs + sigma_abs")

model <- odemodel(f)
p <- P(
  trafo = eqnvec(
    A = "1", B = "0", C = "0", k1 = "exp(log_k1)", k2 = "exp(log_k2)", # dynamic parameters
    off_B = "exp(log_off_B)", # observation parameters
    sigma_rel = "exp(log_sigma_rel)", # error parameters
    sigma_abs = "exp(log_sigma_abs)" # error parameters
  ),
  condition = "C1"
)
x <- Xs(model)
g <- Y(observables, f)
err <- Y(errors, c(observables, as.eqnvec(f)), attach.input = FALSE, compile = TRUE, modelname = "err")
err <- Y(errors, states = "B_obs", parameters = c("A", "B", "C", "sigma_rel", "sigma_abs", "off_B", "k1", "k2"), attach.input = FALSE, compile = TRUE, modelname = "err")

## Simulate data
ptrue <- c(log_k1 = -1, log_k2 = -2, log_off_B = -1, log_sigma_rel = log(0.1), log_sigma_abs = log(.1))
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

