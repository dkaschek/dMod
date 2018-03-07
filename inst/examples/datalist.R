## Generate datalist from scratch
mydata1 <- data.frame(name = "A",
                      time = 0:5,
                      value = 0:5,
                      sigma = .1)

mydata2 <- data.frame(name = "A",
                      time = 0:5,
                      value = sin(0:5),
                      sigma = .1)

data <- datalist(C1 = mydata1, C2 = mydata2)
print(data)
plot(data)

## Generate datalist from singla data.frame
times <- seq(0, 2*pi, length.out = 20)
mydata <- data.frame(name = "A", 
                     time = times, 
                     value = c(sin(times), 1.5 * sin(times)), 
                     sigma = .1, 
                     stage = rep(c("upper", "lower"), each = 10),
                     phase = rep(c("first", "second"), each = 20),
                     amplitude = rep(c(1,1.5), each = 20))

data <- as.datalist(mydata, split.by = c("stage", "phase"), keep.covariates = "amplitude")
print(data)
plot(data)

condition.grid <- attr(data, "condition.grid")
print(condition.grid)
