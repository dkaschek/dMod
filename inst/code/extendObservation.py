from sympy import *
from sympy.parsing.sympy_parser import *
from sympy.parsing.sympy_tokenize import *

def getObservation(observation, variables, stoichiometry, flows, conserved):
	m = len(variables)
	inversion = [0]*m

	stoichiometry = Matrix(len(stoichiometry)/m,m,stoichiometry)
	stoichiometry = stoichiometry.transpose()

	for i in range(len(flows)):
		flows[i] = flows[i].replace('^','**') 
		flows[i] = parse_expr(flows[i])
	
	diffEquations = list(stoichiometry*Matrix(len(flows),1,flows))

	for v in range(m):
		variables[v] = parse_expr(variables[v])

	#extract observation function from read string
	obsFunctions = []
	for o in range(len(observation)):
		observation[o] = str(observation[o])
		obsFunctions.append(parse_expr(observation[o][(observation[o].find('=')+1):len(observation[o])]))

	#extract observables and parameters from read observation
	global newLine, l , observables, parameters
	newLine = True
	observables = []
	parameters = []
	global l
	l = -1

	def read():
		global l
		l += 1
		if l < len(observation):
			return observation[l] + '\n'
		else:
			raise StopIteration

	def useToken(key, value, Coord1, Coord2, fullLine):
		global newLine, observables, obsFunctions, parameters
		
		if key == 1: #1: NAME  2: NUMBER  51: OP   4: NEWLINE  0: ENDMARKER  
			var(value)
			if newLine == True:
				observables.append(parse_expr(value))
			else:
				parameters.append(parse_expr(value))
			newLine = False
		elif key == 4:
			newLine = True

	tokenize(read, useToken)
	variablesMatrix = Matrix(m,1,variables)
	h = len(observables)

	#calculate conserved quantities
	conservedBase = conservedQuantities(stoichiometry) #base vectors in columns of matrix
	conservedBase = conservedBase.transpose()

	#calculate jacobian of observation
	jacobian = zeros(h, m)
	for i in range(h):
		for j in range(m):
			jacobian[i,j] = diff(obsFunctions[i], variables[j])

	#check if observation functions are linearly dependent
	_, pivots = jacobian.rref(simplify = True)
	if len(pivots) < h:
		print 'Error: Observation functions are linearly dependent.'

	#first extend with conserved quantities
	if conserved:
		for i in range(conservedBase.rows):
			jacobianTemp = jacobian.col_join(conservedBase[i,:])

			_, pivots = jacobianTemp.rref(simplify = True)
			if len(pivots) < jacobianTemp.rows:
				continue
			else:
				jacobian = jacobianTemp
				observables.append(parse_expr('CONST'+str(i+1)))

	#finally extend with ones
	s = len(observables)
	l = jacobian.rows
	jacobian = jacobian.col_join(zeros(m-jacobian.rows, m))
	for i in range(m):
		if not i in pivots:
			jacobian[l, i] = 1
			l += 1
			observables.append(parse_expr(str(variables[i])+'OBS'))
			inversion[i] = observables[-1]

	obsFunctions = obsFunctions + list(jacobian[h:,:]*variablesMatrix)
	
	#substitute trivial inversions in observation functions
	obsFunctionsTemp = obsFunctions[:]
	for i in range(s):
		for j in range(s,m):
			obsFunctionsTemp[i] = obsFunctionsTemp[i].subs(obsFunctionsTemp[j], observables[j])

	#invert nontrivial part
	vars = []
	for i in range(len(inversion)):
		if inversion[i] == 0:
			vars.append(variables[i])
	eq = []
	for i in range(s):
		eq.append(obsFunctionsTemp[i]-observables[i])

	result = solve(eq, vars, dict=True)
	if len(result)>1: 
		print 'Warning: Inversion of observation not unique'
	result = result[0]

	for v in vars:
		inversion[variables.index(v)] = result[v]

	#build new differential equations
	newDiffEquations = list(jacobian*Matrix(m,1,diffEquations))
	for i in range(m):
		for j in range(m):
			newDiffEquations[i] = newDiffEquations[i].subs(variables[j], inversion[j])
		newDiffEquations[i] = simplify(newDiffEquations[i])

	#convert output to strings
	for i in range(m):
		observables[i] = str(observables[i])
		obsFunctions[i] = str(obsFunctions[i])
		newDiffEquations[i] = str(newDiffEquations[i])
		newDiffEquations[i] = newDiffEquations[i].replace('**','^')
		inversion[i] = str(inversion[i])

	return observables, observables[:h], obsFunctions, newDiffEquations, inversion

