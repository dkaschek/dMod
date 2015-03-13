import sys
import argparse
import time

from sympy import *
from sympy.parsing.sympy_parser import parse_expr


try:
	import scipy.linalg
	scipyAvailable = True
except:
	scipyAvailable = False

t0 = time.time()

#sys.setrecursionlimit(10000)
var('epsilon')

def symmetryDetection(model, observation, prediction, initial, ansatz, pMax, inputs, fixed, Parallel, allTrafos):
	
	if fixed == None:
		fixed = []
	elif isinstance(fixed, basestring):
		fixed = [str(fixed)]
	if len(fixed) != 0:
		s = 'Fixed variables: '
		for v in range(len(fixed)):
			s = s + str(fixed[v]) + ', '
			fixed[v] = parse_expr(fixed[v])

		sys.stdout.write(s[0:len(s)-2] + '\n')
		sys.stdout.flush()

	if inputs == None:
		inputs = []
	elif isinstance(inputs, basestring):
		inputs = [str(inputs)]
	if len(inputs) != 0:
		s = 'Input variables: '
		for v in range(len(inputs)):
			s = s + str(inputs[v]) + ', '
			inputs[v] = parse_expr(inputs[v])


		sys.stdout.write(s[0:len(s)-2] + '\n')
		sys.stdout.flush()

	###########################################################################################
	##########################     read data from files    ####################################
	###########################################################################################
	sys.stdout.write('Reading files...')
	sys.stdout.flush()

	# read model
	try:
		variables, diffEquations, parameters = readEqVector(model)
	except IOError:
		print '\nError: No such file or directory: ' + str(model)
		return

	# read observation
	try:
		observables, obsFunctions, obsParameters = readEqVector(observation)
	except IOError:
		print '\nError: No such file or directory: ' + str(observation)
		return

	# read initial values
	if initial != False:
		if isinstance(initial, basestring):
			initial = [str(initial)]

		try:
			initVars, initFunctions, initParameters = readEqVector(initial)
		except IOError:
			print '\nError: No such file or directory: ' + str(initial)
			return

		i = 0
		while i < len(initParameters):
			if initParameters[i] in variables:
				initParameters.pop(i)
			elif initParameters[i] in parameters:
				initParameters.pop(i)
			elif initParameters[i] in obsParameters:
				initParameters.pop(i)
			else:
				i += 1

		o = len(initVars)
		for i in range(o):
			if initVars[i] == initFunctions[i]:
				initFunctions[i] = parse_expr(str(initVars[i])+'_0')
				initParameters.append(parse_expr(str(initVars[i])+'_0'))
		substituted = True
		counter = 0
		while substituted:
			substituted = False
			for k in range(o):
				for j in range(o):
					if initFunctions[k].has(initVars[j]):
						initFunctions[k] = initFunctions[k].subs(initVars[j],initFunctions[j])
						substituted = True
			counter += 1
			if counter > 100:
				print '\nError: There seems to be an infinite recursion in the initial value functions'
				return

		parameters = parameters + initParameters
		
	else:
		initVars, initFunctions, initParameters = [], [], []
		o = 0

	# read predictions
	if not prediction == False:
		if isinstance(prediction, basestring):
			prediction = [str(prediction)]
		try:
			predictions, predFunctions, predParameters = readEqVector(prediction)
		except IOError:
			print '\nError: No such file or directory: ' + str(prediction)
			return

		i = 0
		while i < len(predParameters):
			if predParameters[i] in variables:
				predParameters.pop(i)
			elif predParameters[i] in parameters:
				predParameters.pop(i)
			elif predParameters[i] in obsParameters:
				predParameters.pop(i)
			else:
				i += 1
		if len(predParameters) != 0:
			print '\nError: New parameters occured in the predictions: ' + str(predParameters)
			return

	# remove inputs from parameters
	for par in inputs:
		if par in parameters:
			parameters.remove(par)
		if par in obsParameters:
			obsParameters.remove(par)

	#remove dynamic parameters and variables from observation Parameters
	for var in variables:
		if var in obsParameters:
			obsParameters.remove(var)
	for par in parameters:
		if par in obsParameters:
			obsParameters.remove(par)

	#define some stuff
	allVariables = variables + inputs + parameters + obsParameters
	n = len(allVariables)
	m = len(diffEquations)
	h = len(observables)

	sys.stdout.write('done\n')
	sys.stdout.flush()

	###########################################################################################
	#############################     prepare equations    ####################################
	###########################################################################################
	sys.stdout.write('Preparing equations...')
	sys.stdout.flush()

	# make infinitesimal ansatz
	infis, diffInfis, rs = makeAnsatz(ansatz, allVariables, m, len(inputs), pMax, fixed)

	### extract numerator and denominator of equations
	#differential equations
	numerators = []
	denominators = []
	for k in range(m):
		rational = together(diffEquations[k])
		numerators.append(numer(rational))
		denominators.append(denom(rational))

	#observation functions
	obsNumerators = []
	obsDenominatros = []
	for k in range(h):
		rational = together(obsFunctions[k])
		obsNumerators.append(numer(rational))
		obsDenominatros.append(denom(rational))

	#initial functions
	initNumerators = []
	initDenominatros = []
	for k in range(o):
		rational = together(initFunctions[k])
		initNumerators.append(numer(rational))
		initDenominatros.append(denom(rational))


	### calculate numerator of derivatives of equations
	#differential equatioins
	derivativesNum = [0]*m
	for i in range(m):
		derivativesNum[i] = [0]*n
	for k in range(m):
		for l in range(n):
			if numerators[k].has(allVariables[l]) or denominators[k].has(allVariables[l]):
				derivativesNum[k][l] = diff(numerators[k],allVariables[l]) * denominators[k] - numerators[k] * diff(denominators[k],allVariables[l])
			else:
				derivativesNum[k][l] = 0

	#observation functions
	obsDerivativesNum = [0]*h
	for i in range(h):
		obsDerivativesNum[i] = [0]*n
	for k in range(h):
		for l in range(n):
			if obsNumerators[k].has(allVariables[l]) or obsDenominatros[k].has(allVariables[l]):
				obsDerivativesNum[k][l] = diff(obsNumerators[k],allVariables[l]) * obsDenominatros[k] - obsNumerators[k] * diff(obsDenominatros[k],allVariables[l])
			else:
				obsDerivativesNum[k][l] = 0

	#initial functions
	initDerivativesNum = [0]*len(initFunctions)
	for i in range(o):
		initDerivativesNum[i] = [0]*n
	for k in range(o):
		for l in range(n):
			if initNumerators[k].has(allVariables[l]) or initDenominatros[k].has(allVariables[l]):
				initDerivativesNum[k][l] = diff(initNumerators[k],allVariables[l]) * initDenominatros[k] - initNumerators[k] * diff(initDenominatros[k],allVariables[l])
			else:
				initDerivativesNum[k][l] = 0

	sys.stdout.write('done\n')
	sys.stdout.flush()


	###########################################################################################
	############################     build linear system    ###################################
	###########################################################################################
	sys.stdout.write('\nBuilding system...')
	sys.stdout.flush()

	rSystem = buildSystem(numerators, denominators, derivativesNum, obsDerivativesNum,
				initVars, initDenominatros, initDerivativesNum, initFunctions, 
				infis, diffInfis, allVariables, rs, int(Parallel), ansatz)

	sys.stdout.write('done\n')
	sys.stdout.flush()


	###########################################################################################
	##############################     solve system    ########################################
	###########################################################################################
	sys.stdout.write('\nSolving system of size ' + str(rSystem.rows) + 'x' + str(rSystem.cols) + '...')
	sys.stdout.flush()

	if scipyAvailable:
		#get LU decomposition from scipy
		rSystem = Matrix(scipy.linalg.lu(rSystem, permute_l=True)[1])

		#calculate reduced row echelon form
		pivots = []
		pivotLines = []
		i = -1
		for j in xrange(rSystem.cols):
			if rSystem[j,j] == 0:
				k = 1
				while j-k > i:
					if rSystem[j-k,j] != 0:
						i = j-k
						break
					k += 1
				else:
					k = i-1
					while k >= 0:
						if rSystem[k,j] != 0 and (not k in pivotLines):
							rSystem.row_swap(j,k)
							i = j
							break
						k -= 1
					else: continue
			else:
				i = j

			pivots.append(j)
			pivotLines.append(i)

			coeff = rSystem[i,j]
			rSystem.row(i, lambda x, _: x/coeff)

			for k in xrange(i):
				coeff = rSystem[k,j]
				if coeff != 0:
					rSystem.row(k, lambda x, l: simplify(x - coeff*rSystem[i,l]))

		j = rSystem.rows-1
		while j >= 0:
			if not j in pivotLines:
				rSystem.row_del(j)
			j -= 1
	else:
		rSystem, pivots, swaps = rrefFF(rSystem)
		rSystem = rSystem[0:len(pivots),:]

	sys.stdout.write('done\n')
	sys.stdout.flush()

	###########################################################################################
	#############################     process results    ######################################
	###########################################################################################
	sys.stdout.write('\nProcessing results...')
	sys.stdout.flush()

	# calculate solution space
	sys.stdout.write('\n  calculating solution space')
	sys.stdout.flush()
	baseMatrix = nullSpace(rSystem, pivots)

	# rrefFF switched the order of the rs, so reverse this
	if not scipyAvailable:
		for i in range(len(swaps)):
			swap = swaps.pop()
			baseMatrix.row_swap(swap[0],swap[1])

	#convert infinitesimals to explicit sympy expressions
	sys.stdout.write('\n  convert infinitesimals')
	sys.stdout.flush()
	for i in range(len(infis)):
		infis[i] = infis[i].as_expr()

	#substitute solutions into infinitesimals
	#(and remove the ones with common parameter factors)
	sys.stdout.write('\n  substituting solutions')
	sys.stdout.flush()
	infisAll = []
	nA = 0
	for l in range(baseMatrix.cols):
		infisTmp = infis[:]
		for i in range(len(allVariables)):
			for r in range(len(rs)):
				infisTmp[i] = infisTmp[i].subs(rs[r],baseMatrix[r,l])

		if not allTrafos:
			#extract all factors from first infinitesimal
			for i in range(len(allVariables)):
				if infisTmp[i] != 0:
					fac = factor(infisTmp[i])
					if type(fac) == type(epsilon+1):
						factors = [infisTmp[i]]
					elif type(fac) == type(epsilon):
						factors = [fac]
					else:
						factors = list(fac.args)
				
					break

			i = 0
			while i < len(factors):
				if factors[i].is_number:
					factors.pop(i)
				elif factors[i] in variables:
					factors.pop(i)
				elif type(factors[i]) == type(epsilon+1):
					factors.pop(i)				
				elif type(factors[i]) == type(epsilon**2):
					if type(factors[i].args[0]) != type(epsilon+1):
						factors[i] = factors[i].args[0]
						i += 1
					else:
						factors.pop(i)					
				else:
					i += 1
		
			#check which of the factors is in all other infinitesimals
			for i in range(1,len(infisTmp)):
				if infisTmp[i] == 0: continue
				fac = factor(infisTmp[i])
				if type(fac) == type(epsilon+1):
					factorsTmp = [fac]
				elif type(fac) == type(epsilon):
					factorsTmp = [fac]
				else:
					factorsTmp = list(fac.args)

				j = 0
				while j < len(factors):
					k = 0
					while k < len(factorsTmp):
						if factorsTmp[k].is_number:
							factorsTmp.pop(k)
						elif factorsTmp[k] in variables:
							factorsTmp.pop(k)
						elif type(factorsTmp[k]) == type(epsilon+1):
							factorsTmp.pop(k)				
						elif type(factorsTmp[k]) == type(epsilon**2):
							if type(factorsTmp[k].args[0]) != type(epsilon+1):
								factorsTmp[k] = factorsTmp[k].args[0]
								k += 1
							else:
								factorsTmp.pop(k)					
						else:
							k += 1
					if factors[j] in factorsTmp:
						j += 1
						continue
					else:
						factors.pop(j)

				if len(factors) != 0:
					continue	#if potential common factors are left, try next ifinitesimal
				else:
					break		#otherwise treat next solution
		
			if len(factors) == 0:
				nA += 1
				infisAll.append(infisTmp)
		else:
			infisAll.append(infisTmp)
	print ''
	sys.stdout.write('done\n')
	sys.stdout.flush()



	# print transformations
	print '\n\n' + str(len(infisAll)) + ' transformation(s) found:'
	if len(infisAll) != 0: printTransformations(infisAll, allVariables)

	# check solutions
	if True and len(infisAll) != 0:#args.pobab:
		sys.stdout.write('\nChecking solutions...')
		sys.stdout.flush()

		l = []
		for i in range(len(infisAll)):
			admitted, k, expr = checkSolution(numerators, denominators, derivativesNum,
							obsNumerators, obsDenominatros, obsDerivativesNum,
							allVariables, infisAll[i], ansatz)

			if not admitted:
				l.append(i+1)

		if not len(l) == 0:
			string = ''
			for i in range(len(l)):
				string += '#' + str(l[i]) + ' '
			print "\nError: Symmetries failed: " + string + '\n'
			return

		sys.stdout.write('done\n')
		sys.stdout.flush()

	t4 = time.time()
	print time.strftime('\nTotal time: %Hh:%Mm:%Ss', time.gmtime(t4-t0))


	# check predictions
	if not prediction == False:
		print '\nChecking predictions:'
		for i in range(len(predictions)):
			print ' Prediction: ' + str(predictions[i])

			admits = True
			for j in range(len(infisAll)):
				infiPred = 0
				for k in range(n):
					if infisAll[j][k] != 0:
						infiPred += infisAll[j][k] * diff(predFunctions[i], allVariables[k])
		
				infiPred = simplify(infiPred)
				if infiPred != 0:
					admits = False
					p = Wild('p',exclude=[0])
					c = Wild('c')
					if infiPred.match(c*predFunctions[i]**p) != None:
						matches = infiPred.match(c*predFunctions[i]**p)
						print '  #' + str(j+1) + ': ' + str((c*predictions[i]**p).subs(c,matches[c]).subs(p,matches[p]))
					elif infiPred.match(c*(-1*predFunctions[i])**p) != None:
						matches = infiPred.match(c*(-1*predFunctions[i])**p)
						print '  #' + str(j+1) + ': ' + str((c*(-1)**p*predictions[i]**p).subs(c,matches[c]).subs(p,matches[p]))
					else:
						print '  #' + str(j+1) + ': ' + str(infiPred)

			if admits: print '  Admits all transformations'
			print ''
	


