## Observation function
fn <- eqnvec(
  sine = "1 + sin(6.28*omega*time)",
  cosine = "cos(6.28*omega*time)"
)
g <- Y(fn, parameters = "omega")

## Prediction function for time
x <- Xt()

## Parameter transformations to split conditions
p <- NULL
for (i in 1:3) {
  p <- p + P(trafo = c(omega = paste0("omega_", i)), condition = paste0("frequency_", i))
}

## Evaluate prediction
times <- seq(0, 1, .01)
pars <- structure(seq(1, 2, length.out = 3), names = attr(p, "parameters"))

prediction <- (g*x*p)(times, pars)

## Plotting prediction
plot(prediction)
plotPrediction(prediction)
plotPrediction(prediction, scales = "fixed")
plotPrediction(prediction, facet = "grid")
plotPrediction(prediction, 
               scales = "fixed",
               transform = list(sine = "x^2", cosine = "x - 1"))

## Simulate data
dataset <- wide2long(prediction)
dataset <- dataset[seq(1, nrow(dataset), 5),]
dataset$value <- dataset$value + rnorm(nrow(dataset), 0, .1)
dataset$sigma <- 0.1
data <- as.datalist(dataset, split.by = "condition")

## Plotting data
plot(data)
plotData(data)

## Plotting data and prediction with subsetting
plot(prediction, data)
plotCombined(prediction, data)
plotCombined(prediction, data, 
             time <= 0.5 & condition == "frequency_1")
plotCombined(prediction, data, 
             time <= 0.5 & condition != "frequency_1", 
             facet = "grid")
