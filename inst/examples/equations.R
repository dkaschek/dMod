## Generate another equation list
eq <- eqnlist()
eq <- addReaction(eq, "A", "pA", "act_A * A * stimulus", "Phosphorylation of A")
eq <- addReaction(eq, "pA", "A", "deact_A * pA", "Deposphorylation of pA")
eq <- addReaction(eq, "B", "pB", "act_B * B * pA", "Phosphorylation of B")
eq <- addReaction(eq, "pB", "B", "deact_B * pB", "Deposphorylation of pB")

## Extract data.frame of reactions
getReactions(eq)

## Get conserved quantities
conservedQuantities(eq$smatrix)

## Get fluxes
getFluxes(eq)

## Subsetting of equation list
subset(eq, "pB" %in% Product)
subset(eq, grepl("Phosphorylation", Description))

## Time derivatives of observables
observables <- eqnvec(pA_obs = "s1*pA", tA_obs = "s2*(A + pA)")
dobs <- dot(observables, eq)

## Combined equation vector for ODE and observables
f <- c(as.eqnvec(eq), dobs)
print(f)
