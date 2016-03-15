## Generate a prediction function

times <- 0:5
grid <- data.frame(name = "A", time = times, row.names = paste0("p", times))
x <- Xd(grid, condition = "C1")

pars <- structure(rep(0, nrow(grid)), names = row.names(grid))

## Simulate data
data.list <- lapply(1:3, function(i) {
  prediction <- x(times, pars + rnorm(length(pars), 0, 1))
  cbind(wide2long(prediction), sigma = 1)
})

data <- as.datalist(do.call(rbind, data.list))

## Generate objective function based on data and model
## Then fit the data and plot the result
obj <- normL2(data, x)
myfit <- trust(obj, pars, rinit = 1, rmax = 10)
plot(x(times, myfit$argument), data)
