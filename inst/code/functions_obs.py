
from sympy import *
from sympy.parsing.sympy_parser import parse_expr
try:
	import readline
	readlineAvailable = True
except:
	readlineAvailable = False
	

var('epsilon')
var('t')

#returns a matrix of base vectors of the null space given a matrix in rref
#the base vectors are in the columns of the matrix
def nullSpace(matrix, pivots):
	m = matrix.cols

	notPivots = []
	solutions = zeros(m, m-len(pivots))

	i, k, l = m-1, 0, matrix.rows-1
	while i >= 0:
		if i in pivots:
			for h in range(len(notPivots)):
				solutions[i,h] = - matrix[l,notPivots[h]]
			l -= 1
		else:	
			notPivots.append(i)
			solutions[i,k] = 1
			k += 1
		i -= 1
	
	return solutions

#from stoichiometry matrix, calculate conserved quantities
def conservedQuantities(stoichiometry):
	stoiSpace, pivots = (stoichiometry.transpose()).rref()
	stoiSpace = stoiSpace[0:len(pivots),:]#base vectors in rows
	conservedBase = nullSpace(stoiSpace,pivots)#base vectors in columns
	
	return conservedBase

#enables use of simplify in R
def simplifyWrapper(expr):
	if type(expr) == type([]):
		for e in range(len(expr)):
			expr[e] = expr[e].replace('.','_6174823504_').replace('^','**')
			#expr[e] = str(simplify(parse_expr(expr[e])))
			expr[e] = str(together(parse_expr(expr[e])))
			expr[e] = expr[e].replace('_6174823504_','.').replace('**','^')
		return expr
	else:
		expr = expr.replace('.','_6174823504_').replace('^','**')
		#expr = str(simplify(parse_expr(expr)))
		expr = str(together(parse_expr(expr)))
		expr = expr.replace('_6174823504_','.').replace('**','^')
		return expr


	
	

