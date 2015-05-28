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

runbg <- function(expr, filename = "tmp") {
  
  # Save current workspace
  save.image(file = paste0(filename, ".RData"))
  
  # Get loaded packages
  pack <- sapply(strsplit(search(), "package:", fixed = TRUE), function(v) v[2])
  pack <- pack[!is.na(pack)]
  pack <- paste(paste0("library(", pack, ")"), collapse = "\n")
  
  # Write program into character
  program <- paste(
    pack,
    paste0("load('", filename, ".RData')"),
    as.character(expr),
    paste0("save.image(file = '", filename, ".RData')"),
    sep = "\n"
  )
  
  # Write program code into file
  cat(program, file = paste0(filename, ".R"))
  
  # Run in background
  system(paste0("R CMD BATCH ", filename, ".R"), intern = FALSE, wait = FALSE)
  
  out <- parse(text = paste0("load('", filename, ".RData')"))
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
