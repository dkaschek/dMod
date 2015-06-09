similarity <- function(fitlist, x, times, fixed, deriv = FALSE) {
  
  combinations <- combn(1:nrow(fitlist), 2)
  out <- do.call(rbind, lapply(1:ncol(combinations), function(j) {
    
    i1 <- combinations[1, j]
    i2 <- combinations[2, j]
    
    par1 <- unlist(fitlist[i1, -1])
    par2 <- unlist(fitlist[i2, -1])
    
    pred1 <- wide2long(x(times, par1, fixed = fixed, deriv = deriv))
    pred2 <- wide2long(x(times, par2, fixed = fixed, deriv = deriv))
    
    rss <- sum(((pred1$value - pred2$value)/(pred1$value + pred2$value + 1))^2)

    data.frame(i1 = c(i1, i2), i2 = c(i2, i1), rss = rss)    
    
  }))
  
  return(out)
  
}



coupleReactions <- function(what, to, f, couplingParameter = "alpha") {
  
  
  
}

mergeReactions <- function(what, f) {
  
  
  description <- attr(f, "description")
  rates <- attr(f, "rates")
  S <- attr(f, "SMatrix")
  volumes <- attr(f, "volumes")
  
  S[is.na(S)] <- 0
  s <- S[what,]
  s.merged <- apply(s, 2, sum)
  
  rest <- (1:length(rates))[-unique(what)[-1]]
  description <- description[rest]
  rates <- rates[rest]
  volumes <- volumes[rest]
  S <- rbind(S[rest[rest < unique(what)[1]], ],
             s.merged,
             S[rest[rest > unique(what)[1]], ])
  
  S <- S[, apply(S, 2, function(v) any(v != 0))]
  S[S == 0] <- NA
  rownames(S) <- NULL
  
  data <- data.frame(Description = description, Rate = rates, S)

  generateEquations(data, volumes = volumes)
  
  
  
  
}

myexpr <- expression({
  for(i in 1:1e3) {
    M <- matrix(rnorm(1e4, 0, 1), 100, 100)
    m <- solve(M)
  }
})

runbg <- function(expr, filename = "tmp", machine = "localhost") {
  
  # Save current workspace
  save.image(file = paste0(filename, ".RData"))
  
  # Get loaded packages
  pack <- sapply(strsplit(search(), "package:", fixed = TRUE), function(v) v[2])
  pack <- pack[!is.na(pack)]
  pack <- paste(paste0("library(", pack, ")"), collapse = "\n")
  
  # Write program into character
  program <- paste(
    #pack,
    paste0("setwd('~/", filename, "_folder')"),
    #paste0("load('", filename, ".RData')"),
    as.character(expr),
    paste0("save.image(file = '", filename, ".RData')"),
    sep = "\n"
  )
  
  # Write program code into file
  cat(program, file = paste0(filename, ".R"))
  
  # Copy files to temporal folder
  system(paste0("scp -r ", getwd(), " ", machine, ":", filename, "_folder"))
  
  # Run in background
  system(paste0("ssh ", machine, " R CMD BATCH ", filename, "_folder/", filename, ".R --vanilla"), intern = FALSE, wait = FALSE)
  
  out <- structure(vector("list", 2), names = c("get", "purge"))
  
  out[[1]] <- parse(text = paste(sep = "\n",
                            paste0("system('scp ", machine, ":", filename, "_folder/", filename, ".RData ", getwd(), "/')"),
                            paste0("load('", filename, ".RData')")))
  
  out[[2]] <- parse(text = paste(sep = "\n",
                                 paste0("system('ssh ", machine, " rm -r ", filename, "_folder')")
                                 ))

  attr(out, "code") <- program
  
  return(out)
  
}

myexpr <- expression({
  
  sigSq <- signature(n="integer", x="numeric")
  codeSq <- "
  for (int i=0; i < *n; i++) {
    x[i] = x[i]*x[i];
  }"
  sigQd <- signature(n="integer", x="numeric")
  codeQd <- "
  squarefn(n, x);
  squarefn(n, x);
"

fns <- cfunction( list(squarefn=sigSq, quadfn=sigQd), 
                  list(codeSq, codeQd), 
                  convention=".C")

squarefn <- fns[["squarefn"]]
quadfn <- fns[["quadfn"]]

x <- squarefn(10, 3)
})
