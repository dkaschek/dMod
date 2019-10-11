## Generate a prediction function
regfn <- c(y = "sin(a+b+c+d+e+f+g+h+i+j+k+l+m+n+o+p+q+r+s+t+u+v+w+x+y+z*time)")

g <- Y(regfn, parameters = letters)
x <- Xt(condition = "C1")

## Generate data
data <- datalist(
  C1 = data.frame(
    name = "y",
    time = 1:5,
    value = sin(1:5) + rnorm(5, 0, .1),
    sigma = .1
  )
)

## Initialize parameters and time 
pars <- setNames(seq(0,.5, length.out = length(letters)), letters)
times <- seq(0, 5, .1)

plot((g*x)(times, pars), data)

## Do many fits from random positions and store them into parlist
out <- as.parlist(lapply(1:100, function(i) {
  trust(normL2(data, g*x), pars + rnorm(length(pars), 0, 1), rinit = 1, rmax = 10)
}))

summary(out)

## Reduce parlist to parframe
parframe <- as.parframe(out)
plotValues(parframe)

## Reduce parframe to best fit
bestfit <- as.parvec(parframe)
plot((g*x)(times, bestfit), data)

identical(as.parframe(out), as.parframe2.parlist(out))

as.parframe(out) %>% head(6)
as.parframe2.parlist(out)%>% head(6)

as.parframe(out) %>% str1
as.parframe2.parlist(out)%>% str1


rbenchmark::benchmark(as.parframe(out))
rbenchmark::benchmark(as.parframe2.parlist(out))



m_parframe <- data.frame(index = m_idx, 
                         value = vapply(x[m_idx], function(.x) .x$value, 1.0),
                         converged = vapply(x[m_idx], function(.x) .x$converged, TRUE),
                         iterations = vapply(x[m_idx], function(.x) as.integer(.x$iterations), 1L))
m_parframe <- cbind(m_parframe,
                    as.data.frame(t(vapply(x[m_idx], function(.x) .x$argument, x[[1]]$argument))))
