##
##  h o o k e j e e v e s . R  Hooke-Jeeves Minimization Algorithm
##


dfoptim_hjk <- function(par, fn, control = list(), ...) {
    if (!is.numeric(par))
        stop("Argument 'par' must be a numeric vector.", call. = FALSE)
    n <- length(par)
    if (n == 1)
        stop("For univariate functions use some different method.", call. = FALSE)
   
    #-- Control list handling ----------
    cntrl <- list(tol      = 1.e-06,
                  maxfeval = Inf,       # set to Inf if no limit wanted
                  maximize = FALSE,     # set to TRUE  for maximization
                  target   = Inf,       # set to Inf for no restriction
                  info     = FALSE)     # for printing interim information
    nmsCo <- match.arg(names(control), choices = names(cntrl), several.ok = TRUE)
    if (!is.null(names(control))) cntrl[nmsCo] <- control

    tol      <- cntrl$tol;      
    maxfeval <- cntrl$maxfeval
    maximize <- cntrl$maximize
    target   <- cntrl$target
    info     <- cntrl$info

	scale <- if (maximize) -1 else 1
    fun <- match.fun(fn)
    f <- function(x) scale * fun(x, ...)

    #-- Setting steps and stepsize -----
    nsteps <- floor(log2(1/tol))        # number of steps
    steps  <- 2^c(-(0:(nsteps-1)))      # decreasing step size
    dir <- diag(1, n, n)                # orthogonal directions

    x <- par                            # start point
    fx <- fbest <- f(x)                 # smallest value so far
    fcount <- 1                         # counts number of function calls

    if (info) cat("step\tnofc\tfmin\txpar\n")   # info header

    #-- Start the main loop ------------
    ns <- 0
    while (ns < nsteps && fcount < maxfeval && abs(fx) < target) {
        ns <- ns + 1
        hjs    <- .hjsearch(x, f, steps[ns], dir, fcount, maxfeval, target)
        x      <- hjs$x
        fx     <- hjs$fx
        sf     <- hjs$sf
        fcount <- fcount + hjs$finc

        if (info)
            cat(ns, "\t",  fcount, "\t", fx/scale, "\t", x[1], "...\n")
    }

    if (fcount > maxfeval) {
        warning("Function evaluation limit exceeded -- may not converge.")
        conv <- 1
    } else if (abs(fx) > target) {
        warning("Function exceeds min/max value -- may not converge.")
        conv <- 1
    } else {
        conv <- 0
    }

    fx <- fx / scale                    # undo scaling
    return(list(par = x, value = fx,
                convergence = conv, feval = fcount, niter = ns))
}

##  Search with a single scale -----------------------------
.hjsearch <- function(xb, f, h, dir, fcount, maxfeval, target) {
    x  <- xb
    xc <- x
    sf <- 0
    finc <- 0
    hje  <- .hjexplore(xb, xc, f, h, dir)
    x    <- hje$x
    fx   <- hje$fx
    sf   <- hje$sf
    finc <- finc + hje$numf

    # Pattern move
    while (sf == 1) {
        d  <- x-xb
        xb <- x
        xc <- x+d
        fb <- fx
        hje  <- .hjexplore(xb, xc, f, h, dir, fb)
        x    <- hje$x
        fx   <- hje$fx
        sf   <- hje$sf
        finc <- finc + hje$numf

        if (sf == 0) {  # pattern move failed
           hje  <- .hjexplore(xb, xb, f, h, dir, fb)
           x    <- hje$x
           fx   <- hje$fx
           sf   <- hje$sf
           finc <- finc + hje$numf
        }
        if (fcount + finc > maxfeval || abs(fx) > target) break
    }

    return(list(x = x, fx = fx, sf = sf, finc = finc))
}

##  Exploratory move ---------------------------------------
.hjexplore <- function(xb, xc, f, h, dir, fbold) {
    n <- length(xb)
    x <- xb

    if (missing(fbold)) {
        fb <- f(x)
        numf <- 1
    } else {
        fb <- fbold
        numf <- 0
    }

    fx <- fb
    xt <- xc
    sf <- 0                             # do we find a better point ?
    dirh <- h * dir
    fbold <- fx
    for (k in sample.int(n, n)) {       # resample orthogonal directions
        p1 <- xt + dirh[, k]
        ft1 <- f(p1)
        numf <- numf + 1

        p2 <- xt - dirh[, k]
        ft2 <- f(p2)
        numf <- numf + 1

        if (min(ft1, ft2) < fb) {
            sf <- 1
            if (ft1 < ft2) {
                xt <- p1
                fb <- ft1
            } else {
                xt <- p2
                fb <- ft2
            }
        }
    }

    if (sf == 1) {
        x <- xt
        fx <- fb
    }
    return(list(x = x, fx = fx, sf = sf, numf = numf))
}



##
##  h o o k e j e e v e s . R  Hooke-Jeeves Minimization Algorithm
##


dfoptim_hjkb <- function(par, fn, lower = -Inf, upper = Inf, control = list(), ...) {
    if (!is.numeric(par))
        stop("Argument 'par' must be a numeric vector.", call. = FALSE)
    n <- length(par)
    if (n == 1)
        stop("For univariate functions use some different method.", call. = FALSE)

    if(!is.numeric(lower) || !is.numeric(upper))
        stop("Lower and upper limits must be numeric.", call. = FALSE)
    if (length(lower) == 1) lower <- rep(lower, n)
    if (length(upper) == 1) upper <- rep(upper, n)
    if (!all(lower <= upper))
        stop("All lower limits must be smaller than upper limits.", call. = FALSE)
    if (!all(lower <= par) || !all(par <= upper))
        stop("Infeasible starting values -- check limits.", call. = FALSE)
    
   
    #-- Control list handling ----------
    cntrl <- list(tol      = 1.e-06,
                  maxfeval = Inf,       # set to Inf if no limit wanted
                  maximize = FALSE,     # set to TRUE for maximization
                  target   = Inf,       # set to Inf for no restriction
                  info     = FALSE)     # for printing interim information
    nmsCo <- match.arg(names(control), choices = names(cntrl), several.ok = TRUE)
    if (!is.null(names(control))) cntrl[nmsCo] <- control

    tol      <- cntrl$tol;      
    maxfeval <- cntrl$maxfeval
    maximize <- cntrl$maximize
    target   <- cntrl$target
    info     <- cntrl$info

	scale <- if (maximize) -1 else 1
    fun <- match.fun(fn)
    f <- function(x) scale * fun(x, ...)

    #-- Setting steps and stepsize -----
    nsteps <- floor(log2(1/tol))        # number of steps
    steps  <- 2^c(-(0:(nsteps-1)))      # decreasing step size
    dir <- diag(1, n, n)                # orthogonal directions

    x <- par                            # start point
    fx <- fbest <- f(x)                 # smallest value so far
    fcount <- 1                         # counts number of function calls

    if (info) cat("step\tnofc\tfmin\txpar\n")   # info header

    #-- Start the main loop ------------
    ns <- 0
    while (ns < nsteps && fcount < maxfeval && abs(fx) < target) {
        ns <- ns + 1
        hjs    <- .hjbsearch(x, f, lower, upper,
                            steps[ns], dir, fcount, maxfeval, target)
        x      <- hjs$x
        fx     <- hjs$fx
        sf     <- hjs$sf
        fcount <- fcount + hjs$finc

        if (info)
            cat(ns, "\t",  fcount, "\t", fx/scale, "\t", x[1], "...\n")
    }

    if (fcount > maxfeval) {
        warning("Function evaluation limit exceeded -- may not converge.")
        conv <- 1
    } else if (abs(fx) > target) {
        warning("Function exceeds min/max value -- may not converge.")
        conv <- 1
    } else {
        conv <- 0
    }

    fx <- fx / scale                    # undo scaling
    return(list(par = x, value = fx,
                convergence = conv, feval = fcount, niter = ns))
}

##  Search with a single scale -----------------------------
.hjbsearch <- function(xb, f, lo, up, h, dir, fcount, maxfeval, target) {
    x  <- xb
    xc <- x
    sf <- 0
    finc <- 0
    hje  <- .hjbexplore(xb, xc, f, lo, up, h, dir)
    x    <- hje$x
    fx   <- hje$fx
    sf   <- hje$sf
    finc <- finc + hje$numf

    # Pattern move
    while (sf == 1) {
        d  <- x-xb
        xb <- x
        xc <- x+d
        xc <- pmax(pmin(xc, up), lo)
        fb <- fx
        hje  <- .hjbexplore(xb, xc, f, lo, up, h, dir, fb)
        x    <- hje$x
        fx   <- hje$fx
        sf   <- hje$sf
        finc <- finc + hje$numf

        if (sf == 0) {  # pattern move failed
           hje  <- .hjbexplore(xb, xb, f, lo, up, h, dir, fb)
           x    <- hje$x
           fx   <- hje$fx
           sf   <- hje$sf
           finc <- finc + hje$numf
        }
        if (fcount + finc > maxfeval || abs(fx) > target) break
    }

    return(list(x = x, fx = fx, sf = sf, finc = finc))
}

##  Exploratory move ---------------------------------------
.hjbexplore <- function(xb, xc, f, lo, up, h, dir, fbold) {
    n <- length(xb)
    x <- xb

    if (missing(fbold)) {
        fb <- f(x)
        numf <- 1
    } else {
        fb <- fbold
        numf <- 0
    }

    fx <- fb
    xt <- xc
    sf <- 0                             # do we find a better point ?
    dirh <- h * dir
    fbold <- fx
    for (k in sample.int(n, n)) {       # resample orthogonal directions
        p1 <- xt + dirh[, k]
        if ( p1[k] <= up[k] ) {
            ft1 <- f(p1)
            numf <- numf + 1
        } else {
            ft1 <- fb
        }

        p2 <- xt - dirh[, k]
        if ( lo[k] <= p2[k] ) {
            ft2 <- f(p2)
            numf <- numf + 1
        } else {
            ft2 <- fb
        }

        if (min(ft1, ft2) < fb) {
            sf <- 1
            if (ft1 < ft2) {
                xt <- p1
                fb <- ft1
            } else {
                xt <- p2
                fb <- ft2
            }
        }
    }

    if (sf == 1) {
        x  <- xt
        fx <- fb
    }

    return(list(x = x, fx = fx, sf = sf, numf = numf))
}


dfoptim_nmk <-
function(par, fn, control=list(), ...) {
    ctrl <- list(tol=1.e-06, maxfeval = min(5000, max(1500, 20*length(par)^2)), regsimp=TRUE, maximize=FALSE, restarts.max=3, trace=FALSE)
	namc <- match.arg(names(control), choices=names(ctrl), several.ok=TRUE)
    if (!all(namc %in% names(ctrl))) 
        stop("unknown names in control: ", namc[!(namc %in% names(ctrl))])
    if (!is.null(names(control))) ctrl[namc] <- control
    ftol <- ctrl$tol
    maxfeval <- ctrl$maxfeval
    regsimp <- ctrl$regsimp
    restarts.max <- ctrl$restarts.max
    maximize <- ctrl$maximize
    trace <- ctrl$trace

	if (maximize) fnm <- function(par, ...) -fn(par, ...) else fnm <- function(par, ...) fn(par, ...) 

	x0 <- par
	n <- length(par)
	if (n == 1) stop(call. = FALSE, "Use `optimize' for univariate optimization")
 	if (n > 30) warning("Nelder-Mead should not be used for high-dimensional optimization")

	V <- cbind(rep(0, n), diag(n))
	f <- rep(0, n+1)
	f[1] <- fnm(x0, ...)
	V[, 1] <- x0
	scale <- 1 # max(1, sqrt(sum(x0^2)))

	if (regsimp) {
		alpha <- scale / (n * sqrt(2)) * c(sqrt(n+1) + n - 1, sqrt(n+1) -1)
		V[, -1] <- (x0 + alpha[2])
		diag(V[, -1]) <- x0[1:n] + alpha[1]
		for (j in 2:ncol(V)) f[j] <- fnm(V[,j], ...) 
	} else {
		V[, -1] <- x0 + scale * V[, -1] 
		for (j in 2:ncol(V)) f[j] <- fnm(V[,j], ...) 
	}

	f[is.nan(f)] <- Inf

	nf <- n + 1
	ord <- order(f)
	f <- f[ord]
	V <- V[, ord]
	
	rho <- 1
	gamma <- 0.5
	chi <- 2
	sigma <- 0.5
	conv <- 1
	oshrink <- 1
	restarts <- 0
	orth <- 0
	dist <- f[n+1] - f[1]
	
	v <- V[, -1] - V[, 1]
	delf <- f[-1] - f[1]
	diam <- sqrt(colSums(v^2))
	sgrad <- c(solve(t(v), delf))
	alpha <- 1.e-04 * max(diam) / sqrt(sum(sgrad^2))
	simplex.size <- sum(abs(V[, -1] - V[, 1])) / max(1, sum(abs(V[, 1])))

	itc <- 0
	conv <- 0
	message <- "Succesful convergence"

	while (nf < maxfeval & restarts < restarts.max & dist > ftol & simplex.size > 1.e-06) {

		fbc <- mean(f)
		happy <- 0
		itc <- itc + 1
		xbar <- rowMeans(V[, 1:n])
		xr <- (1 + rho) * xbar - rho * V[, n+1]
		fr <- fnm(xr, ...)
		nf <- nf + 1
		if(is.nan(fr)) fr <- Inf

		if (fr >= f[1] & fr < f[n]) {
			happy <- 1
			xnew <- xr
			fnew <- fr
		} else if (fr < f[1]) {
			xe <- (1 + rho * chi) * xbar - rho * chi * V[, n+1]
			fe <- fnm(xe, ...)
		if(is.nan(fe)) fe <- Inf
			nf <- nf + 1
			if (fe < fr) {
				xnew <- xe
				fnew <- fe
				happy <- 1
			} else {
				xnew <- xr
				fnew <- fr
				happy <- 1
			}
		} else if (fr >= f[n] & fr < f[n+1]) {
			xc <- (1 + rho * gamma) * xbar - rho * gamma * V[, n+1]
			fc <- fnm(xc, ...)
		if(is.nan(fc)) fc <- Inf
			nf <- nf + 1
			if (fc <= fr) {
				xnew <- xc
				fnew <- fc
				happy <- 1
			}
		} else if (fr >= f[n+1]) {
			xc <- (1 - gamma) * xbar + gamma * V[, n+1]
			fc <- fnm(xc, ...)
		if(is.nan(fc)) fc <- Inf
			nf <- nf + 1
			if (fc < f[n+1]) {
				xnew <- xc
				fnew <- fc
				happy <- 1
			}
		}

		if (happy == 1 & oshrink == 1) {
			fbt <- mean(c(f[1:n], fnew))
			delfb <- fbt - fbc
			armtst <- alpha * sum(sgrad^2)
			if (delfb > - armtst/n) {
			if (trace) cat("Trouble - restarting: \n")
				restarts <- restarts + 1
				orth <- 1
				diams <- min(diam)
				sx <- sign(0.5 * sign(sgrad))
				happy <- 0
				V[, -1] <- V[, 1]
				diag(V[, -1]) <- diag(V[, -1]) - diams * sx[1:n]
			}
		}

		if (happy == 1) {
			V[, n+1] <- xnew
			f[n+1] <- fnew
			ord <- order(f)
			V <- V[, ord]
			f <- f[ord]
		} else if (happy == 0 & restarts < restarts.max) {
			if (orth == 0) orth <- 1
			V[, -1] <- V[, 1] - sigma * (V [, -1] - V[, 1])
			for (j in 2:ncol(V)) f[j] <- fnm(V[,j], ...)  ## kmm change
			nf <- nf + n
			ord <- order(f)
			V <- V[, ord]
			f <- f[ord]
		}

		v <- V[, -1] - V[, 1]
		delf <- f[-1] - f[1]
		diam <- sqrt(colSums(v^2))
		simplex.size <- sum(abs(v)) / max(1, sum(abs(V[, 1])))

		f[is.nan(f)] <- Inf

		dist <- f[n+1] - f[1]
		sgrad <- c(solve(t(v), delf))
		if (trace & !(itc %% 2)) cat("iter: ", itc, "\n", "value: ", f[1], "\n")
	}

	if (dist <= ftol | simplex.size <= 1.e-06) {
		conv <- 0
		message <- "Successful convergence"
		} else if (nf >= maxfeval) {
		conv <- 1
		message <- "Maximum number of fevals exceeded"
		} else if (restarts >= restarts.max) {
		conv <- 2
		message <- "Stagnation in Nelder-Mead"
		}	

	return(list(par = V[, 1], value=f[1]*(-1)^maximize, feval=nf, restarts=restarts, convergence=conv, message=message))
}



# Nelder-Mead with Box constraints
dfoptim_nmkb <- function (par, fn, lower = -Inf, upper = Inf, control = list(), ...) 
{
    ctrl <- list(tol = 1e-06, maxfeval = min(5000, max(1500, 
        20 * length(par)^2)), regsimp = TRUE, maximize = FALSE, 
        restarts.max = 3, trace = FALSE)
    namc <- match.arg(names(control), choices = names(ctrl), 
        several.ok = TRUE)
    if (!all(namc %in% names(ctrl))) 
        stop("unknown names in control: ", namc[!(namc %in% names(ctrl))])
    if (!is.null(names(control))) 
        ctrl[namc] <- control
    ftol <- ctrl$tol
    maxfeval <- ctrl$maxfeval
    regsimp <- ctrl$regsimp
    restarts.max <- ctrl$restarts.max
    maximize <- ctrl$maximize
    trace <- ctrl$trace
    n <- length(par)

g <- function(x) {
gx <- x
gx[c1] <- atanh(2 * (x[c1] - lower[c1]) / (upper[c1] - lower[c1]) - 1)
gx[c3] <- log(x[c3] - lower[c3])
gx[c4] <- log(upper[c4] - x[c4])
gx
}

ginv <- function(x) {
gix <- x
gix[c1] <- lower[c1] + (upper[c1] - lower[c1])/2 * (1 + tanh(x[c1]))
gix[c3] <- lower[c3] + exp(x[c3])
gix[c4] <- upper[c4] - exp(x[c4])
gix
}

if (length(lower) == 1) lower <- rep(lower, n)
if (length(upper) == 1) upper <- rep(upper, n)

if (any(c(par < lower, upper < par))) stop("Infeasible starting values!", call.=FALSE)

low.finite <- is.finite(lower)
upp.finite <- is.finite(upper)
c1 <- low.finite & upp.finite  # both lower and upper bounds are finite 
c2 <- !(low.finite | upp.finite) # both lower and upper bounds are infinite
c3 <- !(c1 | c2) & low.finite # finite lower bound, but infinite upper bound
c4 <- !(c1 | c2) & upp.finite  # finite upper bound, but infinite lower bound

if (all(c2)) stop("Use `nmk()' for unconstrained optimization!", call.=FALSE)

    if (maximize) 
        fnmb <- function(par) -fn(ginv(par), ...)
    else fnmb <- function(par) fn(ginv(par), ...)
    
    x0 <- g(par)
    if (n == 1) 
        stop(call. = FALSE, "Use `optimize' for univariate optimization")
    if (n > 30) 
        warning("Nelder-Mead should not be used for high-dimensional optimization")
    V <- cbind(rep(0, n), diag(n))
    f <- rep(0, n + 1)
    f[1] <- fnmb(x0)
    V[, 1] <- x0
    scale <- 1 # max(1, sqrt(sum(x0^2)))
    if (regsimp) {
        alpha <- scale/(n * sqrt(2)) * c(sqrt(n + 1) + n - 1, 
            sqrt(n + 1) - 1)
        V[, -1] <- (x0 + alpha[2])
        diag(V[, -1]) <- x0[1:n] + alpha[1]
        for (j in 2:ncol(V)) f[j] <- fnmb(V[, j])
    }
    else {
        V[, -1] <- x0 + scale * V[, -1]
        for (j in 2:ncol(V)) f[j] <- fnmb(V[, j])
    }
    f[is.nan(f)] <- Inf
    nf <- n + 1
    ord <- order(f)
    f <- f[ord]
    V <- V[, ord]
    rho <- 1
    gamma <- 0.5
    chi <- 2
    sigma <- 0.5
    conv <- 1
    oshrink <- 1
    restarts <- 0
    orth <- 0
    dist <- f[n + 1] - f[1]
    v <- V[, -1] - V[, 1]
    delf <- f[-1] - f[1]
    diam <- sqrt(colSums(v^2))
    sgrad <- c(solve(t(v), delf))
    alpha <- 1e-04 * max(diam)/sqrt(sum(sgrad^2))
    simplex.size <- sum(abs(V[, -1] - V[, 1]))/max(1, sum(abs(V[, 
        1])))
    itc <- 0
    conv <- 0
    message <- "Succesful convergence"
    while (nf < maxfeval & restarts < restarts.max & dist > ftol & 
        simplex.size > 1e-06) {
        fbc <- mean(f)
        happy <- 0
        itc <- itc + 1
        xbar <- rowMeans(V[, 1:n])
        xr <- (1 + rho) * xbar - rho * V[, n + 1]
        fr <- fnmb(xr)
        nf <- nf + 1
        if (is.nan(fr)) 
            fr <- Inf
        if (fr >= f[1] & fr < f[n]) {
            happy <- 1
            xnew <- xr
            fnew <- fr
        }
        else if (fr < f[1]) {
            xe <- (1 + rho * chi) * xbar - rho * chi * V[, n + 
                1]
            fe <- fnmb(xe)
            if (is.nan(fe)) 
                fe <- Inf
            nf <- nf + 1
            if (fe < fr) {
                xnew <- xe
                fnew <- fe
                happy <- 1
            }
            else {
                xnew <- xr
                fnew <- fr
                happy <- 1
            }
        }
        else if (fr >= f[n] & fr < f[n + 1]) {
            xc <- (1 + rho * gamma) * xbar - rho * gamma * V[, 
                n + 1]
            fc <- fnmb(xc)
            if (is.nan(fc)) 
                fc <- Inf
            nf <- nf + 1
            if (fc <= fr) {
                xnew <- xc
                fnew <- fc
                happy <- 1
            }
        }
        else if (fr >= f[n + 1]) {
            xc <- (1 - gamma) * xbar + gamma * V[, n + 1]
            fc <- fnmb(xc)
            if (is.nan(fc)) 
                fc <- Inf
            nf <- nf + 1
            if (fc < f[n + 1]) {
                xnew <- xc
                fnew <- fc
                happy <- 1
            }
        }
        if (happy == 1 & oshrink == 1) {
            fbt <- mean(c(f[1:n], fnew))
            delfb <- fbt - fbc
            armtst <- alpha * sum(sgrad^2)
            if (delfb > -armtst/n) {
                if (trace) 
                  cat("Trouble - restarting: \n")
                restarts <- restarts + 1
                orth <- 1
                diams <- min(diam)
                sx <- sign(0.5 * sign(sgrad))
                happy <- 0
                V[, -1] <- V[, 1]
                diag(V[, -1]) <- diag(V[, -1]) - diams * sx[1:n]
            }
        }
        if (happy == 1) {
            V[, n + 1] <- xnew
            f[n + 1] <- fnew
            ord <- order(f)
            V <- V[, ord]
            f <- f[ord]
        }
        else if (happy == 0 & restarts < restarts.max) {
            if (orth == 0) 
                orth <- 1
            V[, -1] <- V[, 1] - sigma * (V[, -1] - V[, 1])
            for (j in 2:ncol(V)) f[j] <- fnmb(V[, j])
            nf <- nf + n
            ord <- order(f)
            V <- V[, ord]
            f <- f[ord]
        }
        v <- V[, -1] - V[, 1]
        delf <- f[-1] - f[1]
        diam <- sqrt(colSums(v^2))
        simplex.size <- sum(abs(v))/max(1, sum(abs(V[, 1])))
        f[is.nan(f)] <- Inf
        dist <- f[n + 1] - f[1]
        sgrad <- c(solve(t(v), delf))
        if (trace & !(itc%%2)) 
            cat("iter: ", itc, "\n", "value: ", f[1], "\n")
    }
    if (dist <= ftol | simplex.size <= 1e-06) {
        conv <- 0
        message <- "Successful convergence"
    }
    else if (nf >= maxfeval) {
        conv <- 1
        message <- "Maximum number of fevals exceeded"
    }
    else if (restarts >= restarts.max) {
        conv <- 2
        message <- "Stagnation in Nelder-Mead"
    }
    return(list(par = ginv(V[, 1]), value = f[1] * (-1)^maximize, feval = nf, 
        restarts = restarts, convergence = conv, message = message))
}

