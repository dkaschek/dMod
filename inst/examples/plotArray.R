## Generate an objective function with a parameter transformation
prior <- structure(rep(0, 5), names = letters[1:5])
obj1 <- constraintL2(mu = prior)
obj2 <- constraintL2(mu = prior, sigma = 3)
p <- P(trafo = structure(paste0("exp(", names(prior), ")"), names = names(prior)))


pouter <- prior + rnorm(length(prior))

## Compute the profile around the prior
myfit <- trust(obj1*p + obj2, prior, rinit = 1, rmax = 10)
myprof <- profile(obj1*p + obj2, myfit$argument, "a", limits = c(-4, 2))
plotProfile(myprof)

grid <- data.frame(name = "parameter", 
                   time = 1:5,  
                   value = as.numeric(pouter), 
                   sigma = .1,
                   row.names = names(pouter))

x <- Xd(grid, condition = "C1")

plotArray(myprof, x, 1:5, data = datalist(C1 = grid))


