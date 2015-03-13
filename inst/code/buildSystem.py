import sys

from sympy import *
from multiprocessing import *

#from functions import *

def buildSystem(numerators, denominators, derivativesNum, obsDerivativesNum,
			initVars, initDenominators, initDerivativesNum, initFunctions, 
			infis, diffInfis, allVariables, rs, nProc, ansatz):

	n = len(allVariables)
	m = len(numerators)
	h = len(obsDerivativesNum)
	o = len(initDerivativesNum)

	### transform expressions to polynomial class
	# infinitesimals
	infisPoly = [0]*n
	diffInfisPoly = [0]*len(diffInfis)
	for i in range(len(diffInfis)):
		diffInfisPoly[i] = [0]*n
	for i in range(n):
		infisPoly[i] = Apoly(infis[i],allVariables)
		for j in range(len(diffInfisPoly)):
			diffInfisPoly[j][i] = Apoly(diffInfis[j][i],allVariables)

	# differential equations
	numeratorsPoly = [0]*m
	denominatorsPoly = [0]*m
	derivativesNumPoly = [0]*m
	for i in range(m):
		derivativesNumPoly[i] = [0]*n
	for k in range(m):
		numeratorsPoly[k] = Apoly(numerators[k],allVariables)
		denominatorsPoly[k] = Apoly(denominators[k],allVariables)
		for l in range(n):
			derivativesNumPoly[k][l] = Apoly(derivativesNum[k][l],allVariables)

	# observation equations
	obsDerivativesNumPoly = [0]*h
	for i in range(h):
		obsDerivativesNumPoly[i] = [0]*n
	for k in range(h):
		for l in range(n):
			obsDerivativesNumPoly[k][l] = Apoly(obsDerivativesNum[k][l],allVariables)

	### calculate conditions for a differential equation
	def doEquation(k, numerators, denominators, derivativesNum, infis,
			diffInfis, allVariables, rs, ansatz, queue):

		n = len(allVariables)
		m = len(numerators)
		
		if ansatz == 'uni' or ansatz == 'par':
			#calculate polynomial
			polynomial = numerators[k].mul(denominators[k]).mul(diffInfis[0][k])
			for l in range(n):
				polynomial = polynomial.sub(derivativesNum[k][l].mul(infis[l]))

		elif ansatz == 'multi':
			polynomial = Apoly(0, allVariables)
			for i in range(m):
				summand = numerators[i].mul(denominators[i])
				for j in range(m):
					if j != i:
						summand = summand.mul(denominators[j])
				summand = summand.mul(diffInfis[k][i])
				polynomial = polynomial.add(summand)
			for i in range(n):
				summand = derivativesNum[k][i]
				for j in range(m):
					if j != k:
						summand = summand.mul(denominators[j])
				summand = summand.mul(infis[i])
				polynomial = polynomial.sub(summand)

		#determine rSystem such that the coefficients vanish
		lgs = zeros(len(polynomial.coefs), len(rs))

		for c in range(len(polynomial.coefs)):
			for l in range(len(rs)):
				lgs[c,l] = diff(polynomial.coefs[c], rs[l])

		queue.put(lgs)

	### calculate conditions for an observation equation
	def doObsEquation(k, obsDerivativesNum, infis, allVariables, rs, queue):
		n = len(allVariables)

		#calculate polynomial
		polynomial = Apoly(0, allVariables)
		for l in range(n):
			polynomial = polynomial.add(obsDerivativesNum[k][l].mul(infis[l]))

		#determine rSystem such that the coefficients vanish
		lgs = zeros(len(polynomial.coefs), len(rs))

		for c in range(len(polynomial.coefs)):
			for l in range(len(rs)):
				lgs[c,l] = diff(polynomial.coefs[c], rs[l])

		queue.put(lgs)

	### calculate conditions for an initial equation
	def doInitEquation(k, initVars, initDenominators, initDerivativesNum,
				initFunctions, infis, allVariables, rs, queue):
		n = len(allVariables)
		o = len(initDerivativesNum)

		#calculate polynomial
		polynomial = 0
		for l in range(n):
			if initVars[k].has(allVariables[l]):
				polynomial = polynomial - infis[l]*initDenominators[k]**2
			polynomial = polynomial + initDerivativesNum[k][l]*infis[l]

		#substitute initial Functions into conditions
		for l in range(o):
			if polynomial.has(initVars[l]):
				polynomial = polynomial.subs(initVars[l], initFunctions[l])
	
		#determine rSystem such that the coefficients vanish
		polynomial = Apoly(polynomial, allVariables)
		lgs = zeros(len(polynomial.coefs), len(rs))

		for c in range(len(polynomial.coefs)):
			for l in range(len(rs)):
				lgs[c,l] = diff(polynomial.coefs[c], rs[l])

		queue.put(lgs)

	### start the calculations for the first equations
	ns = 0
	queue = Queue()
	while ns < min([m+h+o, nProc]):
		if ns < m:
			p = Process(target=doEquation, args=(ns, numeratorsPoly, denominatorsPoly, derivativesNumPoly, infisPoly,
								diffInfisPoly, allVariables, rs, ansatz, queue))
		elif ns < m+h:
			p = Process(target=doObsEquation, args=(ns-m, obsDerivativesNumPoly, infisPoly, allVariables, rs, queue))
		else:
			p = Process(target=doInitEquation, args=(ns-m-h, initVars, initDenominators, initDerivativesNum,
								initFunctions, infis, allVariables, rs, queue))
		p.start()
		ns += 1


	sys.stdout.write("\rBuilding system...0%")
	sys.stdout.flush()

	### wait till a process has finished and start the calculation for a new equation
	lgsList = []
	lgsSize = 0
	finished = 0
	while ns < m+h+o:
		lgs = queue.get()
		if ns < m:
			p = Process(target=doEquation, args=(ns,numeratorsPoly, denominatorsPoly, derivativesNumPoly, infisPoly,
								diffInfisPoly, allVariables, rs, ansatz, queue))
		elif ns < m+h:
			p = Process(target=doObsEquation, args=(ns-m, obsDerivativesNumPoly, infisPoly, allVariables, rs, queue))
		else:
			p = Process(target=doInitEquation, args=(ns-m-h, initVars, initDenominators, initDerivativesNum,
								initFunctions, infis, allVariables, rs, queue))
		p.start()
		ns += 1

		lgsList.append(lgs)
		lgsSize += lgs.rows
	
		finished += 1

		prog = int(float(finished)/(m+h+o)*100)
		sys.stdout.write("\rBuilding system...%d%%" %prog)
		sys.stdout.flush()

	### wait for all processes to finish
	while finished < m+h+o:
		lgs = queue.get()

		lgsList.append(lgs)
		lgsSize += lgs.rows
	
		finished += 1

		prog = int(float(finished)/(m+h+o)*100)
		sys.stdout.write("\rBuilding system...%d%%" %prog)
		sys.stdout.flush()
	
	sys.stdout.write("\nCombining system...")
	sys.stdout.flush()	

	### combine all conditions into one matrix
	rSystem = zeros(lgsSize, len(rs))
	pos = 0
	for lgs in lgsList:
		rSystem[pos:pos+lgs.rows, :] = lgs
		pos += lgs.rows

	return rSystem
		















