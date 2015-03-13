import sys
import readline

from sympy import *
from sympy.parsing.sympy_parser import parse_expr

### solve the linear system defined by a matrix (adapted from sympys solve_linear_system)
def rrefFF(system):
	matrix = system[:, :]
	swaps = []

	i, m = 0, matrix.cols

	while i < matrix.rows:
		if i == m:
			# remove trailing rows
			matrix = matrix[:i, :]
			break

		if not matrix[i, i]:
			# there is no pivot in current column
			# so try to find one in other columns
			for k in xrange(i + 1, m):
				if matrix[i, k]:
					break
			else:
				# zero row so we can safely skip it
				matrix.row_del(i)
				if not matrix:
					return None
				continue

			# we want to change the order of colums so
			# the order of variables must also change
			swaps.append((i,k))
			matrix.col_swap(i, k)

		pivot_inv = S.One/matrix[i, i]

		# divide all elements in the current row by the pivot
		matrix.row(i, lambda x, _: simplify(x * pivot_inv))

		for k in xrange(matrix.rows):
			if k != i:
				if matrix[k, i]:
					coeff = matrix[k, i]

					# subtract from the current row the row containing
					# pivot and multiplied by extracted coefficient
					matrix.row(k, lambda x, j: simplify(x - matrix[i, j]*coeff))
		i += 1
	
	return matrix, range(matrix.rows), swaps

### returns a matrix of base vectors of the null space given a matrix in rref
### the base vectors are the columns of the matrix
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


### check if the given infinitesimal generator is indeed admitted
def checkSolution(numerators, denominators, derivativesNum,
			obsNumerators, obsDenominatros, obsDerivativesNum,
			allVariables, infis, ansatz):
	n = len(allVariables)
	m = len(numerators)
	h = len(obsNumerators)
	
	admitted = True
	#for every differential equation
	for k in range(m):
		#calculate polynomial
		if ansatz == 'uni' or ansatz == 'par':
			polynomial = diff(infis[k],allVariables[k])*numerators[k]/denominators[k]
			for l in range(n):
				polynomial -= infis[l]*derivativesNum[k][l]/denominators[k]**2
		elif ansatz == 'multi':
			polynomial = 0
			for i in range(m):
				polynomial += diff(infis[k],allVariables[i])*numerators[i]/denominators[i]
			for l in range(n):
				polynomial -= infis[l]*derivativesNum[k][l]/denominators[k]**2
		polynomial = simplify(polynomial)
		if polynomial != 0:
			return False, k, polynomial
			
	#for every observable
	for k in range(h):
		#calculate polynomial
		polynomial = 0
		for l in range(n):
			polynomial += obsDerivativesNum[k][l]/denominators[k]**2 * infis[l]
		polynomial = simplify(polynomial)
		if polynomial != 0:
			return False, -k, polynomial

	return True, 0, 0

### determine known transformations from infinitesimals
def buildTransformation(infis, allVariables):
	n = len(allVariables)
	var('epsilon')
	
	transformations = [0]*n
	tType = [False]*6 #0: unknown, 1: scaling, 2: translation, 3: MM-like, 4: p>2, 5: generalized translation
	for i in range(n):
		if infis[i] == 0:
			transformations[i] = allVariables[i]
		else:
			poly = Poly(infis[i], allVariables).as_dict()
			monomials = poly.keys()
			coefs = poly.values()
			if len(monomials) == 1:
				p = None
				for j in range(n):
					if monomials[0][j] != 0:
						if j == i and p == None: # p Symmetry
							p = monomials[0][i]
						elif p == None and monomials[0][j] == 1: # 
							p = -1-j							
						else:
							transformations[i] = '-?-'
							tType[0] = True
							break
				else:
					if p == None: # translation
						transformations[i] = allVariables[i] + epsilon*coefs[0]
						tType[2] = True
					elif p <= 0: #
						transformations[i] = allVariables[i] + epsilon*coefs[0] * allVariables[-p-1]
						tType[5] = True
					elif p == 1: # scaling
						transformations[i] = exp(epsilon*coefs[0])*allVariables[i]
						tType[1] = True
					else: # p Symmetry
						transformations[i] = simplify(allVariables[i]/(1-(p-1)*epsilon*allVariables[i]**(p-1))**(sympify(1)/(p-1)))
						if p == 2: tType[3] = True
						else: tType[4] = True
			else:	
				transformations[i] = '-?-'
				tType[0] = True

	string = 'Type: '
	if tType[0]: string += 'unknown, '
	if tType[1]: string += 'scaling, '
	if tType[2]: string += 'translation, '
	if tType[3]: string += 'MM-like, '
	if tType[4]: string += 'p>2, '
	if tType[5]: string += 'gen. tanslation, '

	string = string[0:(len(string)-2)]

	return transformations, string

### print found transformations
def printTransformations(infisAll, allVariables):
	n = len(infisAll[0])

	length1 = 8
	length2 = 13
	length3 = 14
	transformations = [0]*len(infisAll)
	types = [0]*len(infisAll)
	for l in range(len(infisAll)):
		for i in range(n):
			infisAll[l][i] = nsimplify(infisAll[l][i])
		transformations[l], types[l] = buildTransformation(infisAll[l], allVariables)
		
		# count length of strings
		for i in range(n):
			if len(str(allVariables[i])) > length1:
				length1 = len(str(allVariables[i]))
			infisAll[l][i] = nsimplify(infisAll[l][i])
			if len(str(infisAll[l][i])) > length2:
				length2 = len(str(infisAll[l][i]))
			if len(str(transformations[l][i])) > length3:
				length3 = len(str(transformations[l][i]))

	print ('{0:'+str(length1)+'s} : ').format('variable') \
		+ ('{0:'+str(length2)+'s} : ').format('infinitesimal')\
		+ str('transformation')

	for l in range(len(infisAll)):
		print '-'*(length1+length2+length3+6)
		print '#' + str(l+1) + ': ' + types[l]

		for i in range(n):
			if infisAll[l][i] != 0: print ('{0:'+str(length1)+'s} : ').format(str(allVariables[i])) \
							+ ('{0:'+str(length2)+'s} : ').format(str(infisAll[l][i]))\
							+ str(transformations[l][i])

# recursive function to construct a multidimensional polynomial
# vars: variables, i: position in vars, p: degree left for other variables
# summand: current monom under construction, poly: full polynomial
# num: umber of coefficients, k: ansatz for which variable, rs: list of coefficiets
def giveDegree(vars, i, p, summand, poly, num, k, rs):

	if i == len(vars)-1:
		var('r_'+str(vars[k])+'_'+str(num))
		rs.append(parse_expr('r_'+str(vars[k])+'_'+str(num)))
		poly += parse_expr('r_'+str(vars[k])+'_'+str(num))*summand*vars[i]**p
		return poly, num+1
	else:
		for j in range(p+1):
			poly, num = giveDegree(vars, i+1, p-j, summand*vars[i]**j, poly, num, k, rs)

	return poly, num

# make infinitesimal ansatz
def makeAnsatz(ansatz, allVariables, m, q, pMax, fixed):
	n = len(allVariables)

	if ansatz == 'uni':
		#construct polynomial
		rs = []
		infis = []
		for k in range(n):
			infis.append(sympify(0))
			if allVariables[k] in fixed: continue #if in fixed, ansatz is 0
			for p in range(pMax+1):
				var('r_'+str(allVariables[k])+'_'+str(p))
				rs.append(parse_expr('r_'+str(allVariables[k])+'_'+str(p)))
				infis[-1] += rs[-1] * allVariables[k]**p

		#calculate derivatives
		diffInfis = [[0]*n]
		for i in range(n):
			diffInfis[0][i] = diff(infis[i],allVariables[i])

	elif ansatz == 'par':
		rs = []
		infis = []
		for k in range(n):
			infis.append(sympify(0))
			if allVariables[k] in fixed: continue #if in fixed, ansatz is 0
			num = 0
			for p in range(pMax+1): #for every degree for 0 to pMax
				vari = allVariables[m+q:] #all parameters
				if k < (m+q): #if ansatz is not for a 
					vari.append(allVariables[k])
					kp = len(vari)-1
				else:	
					kp = k-(m+q)
				degree, num = giveDegree(vari, 0, p, 1, 0, num, kp, rs)
				infis[-1] += degree

		#calculate derivatives
		diffInfis = [[0]*n]
		for i in range(n):
			diffInfis[0][i] = diff(infis[i],allVariables[i])

	elif ansatz == 'multi':
		rs = []
		infis = []
		for k in range(n):
			infis.append(sympify(0))
			if allVariables[k] in fixed: continue #if in fixed, ansatz is 0
			num = 0
			for p in range(pMax+1): #for every degree for 0 to pMax
				if k < m:  #if ansatz is for a dynamic variable
					vari = allVariables[:m] + allVariables[m+q:]
					kp = k
				elif k < m+q:  #if ansatz is for an input
					vari = allVariables[:]
					kp = k
				else:  #if ansatz is for a parameter	
					vari = allVariables[m+q:]  #all parameters
					kp = k-(m+q)
				degree, num = giveDegree(vari, 0, p, 1, 0, num, kp, rs)
				infis[-1] += degree

		#calculate derivatives
		diffInfis = [0]*n
		for i in range(n):
			diffInfis[i] = [0]*n
		for i in range(n):
			for j in range(n):
				diffInfis[i][j] = diff(infis[i],allVariables[j])

	return infis, diffInfis, rs

### efficient class for polynomial calculations
class Apoly:
	def __init__(self, expr, vars):
		if expr == None:
			self.coefs = []
			self.exps = []
		else:
			poly = Poly(expr, vars).as_dict()
			self.coefs = poly.values()
			self.exps = poly.keys()
			for i in xrange(len(self.exps)):
				self.exps[i] = Matrix(1,len(self.exps[i]),self.exps[i])

	def __repr__(self):
		return str(self.coefs) + '\n' + str(self.exps)
	def __str__(self):
		return str(self.coefs) + '\n' + str(self.exps)

	### add a second polynomial and return the result
	def add(self, otherPoly):
		newPoly = Apoly(None,None)
		newPoly.exps = self.exps[:]
		newPoly.coefs = self.coefs[:]
		for i in xrange(len(otherPoly.exps)):
			for j in xrange(len(newPoly.exps)):
				if otherPoly.exps[i] == newPoly.exps[j]:
					newPoly.coefs[j] = newPoly.coefs[j] + otherPoly.coefs[i]
					if newPoly.coefs[j] == 0:
						newPoly.coefs.pop(j)
						newPoly.exps.pop(j)
					break
			else:
				newPoly.coefs.append(otherPoly.coefs[i])
				newPoly.exps.append(otherPoly.exps[i])

		return newPoly

	### substract a second polynomial and return the result
	def sub(self, otherPoly):
		newPoly = Apoly(None,None)
		newPoly.exps = self.exps[:]
		newPoly.coefs = self.coefs[:]
		for i in xrange(len(otherPoly.exps)):
			for j in xrange(len(newPoly.exps)):
				if otherPoly.exps[i] == newPoly.exps[j]:
					newPoly.coefs[j] = newPoly.coefs[j] - otherPoly.coefs[i]
					if newPoly.coefs[j] == 0:
						newPoly.coefs.pop(j)
						newPoly.exps.pop(j)
					break
			else:
				newPoly.coefs.append(-otherPoly.coefs[i])
				newPoly.exps.append(otherPoly.exps[i])


		return newPoly

	### multiply a second polynomial and return the result
	def mul(self, otherPoly):
		newPoly = Apoly(None,None)
		newPoly.coefs = [0]*(len(self.coefs)*len(otherPoly.coefs))
		newPoly.exps = [0]*(len(self.coefs)*len(otherPoly.coefs))
		k = 0
		for i in xrange(len(otherPoly.exps)):
			for j in xrange(len(self.exps)):
				newPoly.coefs[k] = otherPoly.coefs[i] * self.coefs[j]
				newPoly.exps[k] = otherPoly.exps[i] + self.exps[j]
				k += 1

		i = 0
		while i < len(newPoly.coefs):
			j = i+1
			while j <len(newPoly.coefs):
				if newPoly.exps[i] == newPoly.exps[j]:
					newPoly.exps.pop(j)
					newPoly.coefs[i] = newPoly.coefs[i] + newPoly.coefs.pop(j)
				j += 1
			i += 1
		
		return newPoly

