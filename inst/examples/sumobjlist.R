# library(dMod)
nm1 <- letters[1:3]
nm2 <- letters[2:4]

ol1 <- objlist(1, `names<-`(rep(1,3), nm1), matrix(1, 3,3, F, list(nm1, nm1)))
ol2 <- objlist(2, `names<-`(rep(2,3), nm2), matrix(2, 3,3, F, list(nm2, nm2)))
# same elements, same grad and hessian names
ol1 + ol1
# same elements, different gradient and hessian names
ol1 + ol2
