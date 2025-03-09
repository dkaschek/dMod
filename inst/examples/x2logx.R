library(reticulate)
library(dMod)

# Define an example equation vector
eqnvec <- as.eqnvec(
  c("-k1*A/(Km+A) + k2*B", "k1*A/(Km+A) - k2*B"),
  names = c("A", "B")
)

# Apply the corrected transformation
log_transformed_eqns <- x2logx(eqnvec, c("B"))

# Inspect the transformed equations
print(log_transformed_eqns)
