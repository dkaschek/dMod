# Generate a parameter vector
v <- parvec(a = 2, b = 3)
print(v)
print(derivatives(v))

# Parameter vector from a named numeric
M <- matrix(c(1, 1, 0, 1), 
	    nrow = 2, ncol = 2, 
	    dimnames = list(c("a", "b"), c("A", "B"))
    )
v <- as.parvec(x = c(a = 2, b = 3), deriv = M)
print(v)
print(derivatives(v))

# Subsetting of parameter vectors
# Case 1: Dependencies in the Jacobian are maintained
w <- v[1]
print(w)
print(derivatives(w))

# Case 2: Dependencies are dropped
w <- v[1, drop = TRUE]
print(w)
print(derivatives(w))

# Concatenating parameter vectors
w <- parvec(c = 4, d = 5)
print(c(v, w))
print(derivatives(c(v, w)))
