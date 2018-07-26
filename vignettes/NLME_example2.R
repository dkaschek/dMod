
plot(0, 0, xlim = c(1, 50), ylim = c(0, 250))

startvec <- c(Asym = 200, xmid = 725, scal = 350)
nfun2 <- function(input, Asym, xmid, scal) {
  
  #output <<- list(input, Asym, xmid, scal)
  
  #print(unique(cbind(Asym, xmid, scal)))
  
  value <- Asym/(1+exp((xmid-input)/scal))
  grad <- cbind(Asym=1/(1+exp((xmid-input)/scal)),
                xmid=-Asym/(1+exp((xmid-input)/scal))^2*1/scal*
                  exp((xmid-input)/scal),
                scal=-Asym/(1+exp((xmid-input)/scal))^2*
                  -(xmid-input)/scal^2*exp((xmid-input)/scal))
  attr(value,"gradient") <- grad
  
  counter <<- counter + 1
  
  #if (counter == 1) { 
    #print(cbind(input, Asym, xmid, scal))
  #  matplot(as.vector(grad), pch = 1, col = 1)
  #  mycall <<- 1
  #}
  
  return(value)
  
  
  
}

mycall <- 0
counter <- 0
output <- NULL
nm1c <- nlmer(circumference ~ nfun2(age, Asym, xmid, scal)  ~ Asym | Tree, data = as.data.frame(Orange), start = startvec)



x <- eqnlist() %>% 
  addReaction("0", "circumference", "k1*circumference", "Gain term") %>%
  addReaction("circumference", "0", "k2*circumference^2", "Loss term") %>%
  odemodel() %>% Xs()
controls(x, NULL, "optionsSens") <- list(method = "lsodes", rtol = 1e-8, atol = 1e-8)

parameters <- getParameters(x)
p <- eqnvec() %>% 
  define("x~x", x = parameters) %>%
  insert("circumference~(Asym/(exp(xmid/scal)+1))") %>%
  insert("k1~(1/scal)") %>%
  insert("k2~(1/(Asym*scal))") %>%
  P()

estimate <- getParameters(p)
ini <- c(Asym = 192, xmid = 728, scal = 348)
times <- seq(0, 2000, 10)

(x*p)(times, ini) %>% plot()

# Set up function for NLME
parnames <- names(ini)
nlmeModel <- function(time, ...) {
  
  #output <<- c(list(time), list(...))
  
  
  pars <- as.data.frame(c(list(time, "circumference"), list(...)))
  names(pars) <- c("time", "name", parnames)
  
  # if (counter == 2) {
  #   #print(pars)
  # }
  # 
  
  id <- cumsum(Reduce("|", lapply(pars[-(1:2)], function(x) !duplicated(x))))
  pars <- split(pars, id)
  
  output <- lapply(pars, function(sub) {
    timesD <- unique(sub$time)
    parsD <- unlist(sub[1, parnames])
    #print(parsD)
    prediction <- (x*p)(timesD, parsD)[[1]]
    template <- data.frame(name = sub$name, time = sub$time, value = 0, sigma = 1)
    
    myres <- res(template, prediction)
    
    output <- as.numeric(myres$prediction)
    deriv <- as.matrix(attr(myres, "deriv")[, -(1:2)])
    
    list(output, deriv)
  })
  
  gradient <- do.call(rbind, lapply(output, function(x) x[[2]]))
  value <- as.numeric(unlist(lapply(output, function(x) x[[1]])))
  
  attr(value, "gradient") <- gradient
  
  counter <<- counter + 1
 
 #if (counter == 1) {
#   matplot(x = as.vector(gradient), col = 2, pch = 1, add = TRUE)
   #print(value)
 #  mycall <<- 1
 #}
 
  return(value)
  
  
  
}

mycall <- 0
counter <- 0
nm2c <- nlmer(circumference ~ nlmeModel(age, Asym, xmid, scal)  ~ Asym | Tree, 
              data = as.data.frame(Orange), start = startvec) 

#              control = glmerControl(optimizer = "optimx", calc.derivs = TRUE, 
#                                     optCtrl = list(method = "nlminb", starttests = FALSE, kkt = FALSE)))


fit <- nlme(circumference ~ nlmeModel(age, Asym, xmid, scal), data = Orange,
                  fixed = Asym+xmid+scal ~ 1, random = Asym ~ 1 | Tree,
                  start = c(100, 200, 200), verbose = TRUE)

out <- stats::profile(nm1c)
out <- stats::profile(nm2c)
