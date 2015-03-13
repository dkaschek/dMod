import csv
import token

from sympy import *
from sympy.parsing.sympy_parser import *
from sympy.parsing.sympy_tokenize import *

def readEqVector(vector):
	global l
	l = 0
	def read():
		global l
		if l < len(vector):
			line = str(vector[l])
			l += 1
			line = line.replace('"','')
			line = line.replace(',','')
			line.strip()	
			return line + '\n'
		else:
			raise StopIteration

	global newLine
	newLine = True
	global observables
	observables = []
	global obsFunctions
	obsFunctions = []
	global parameters
	parameters = []
	
	def useToken(key, value, Coord1, Coord2, fullLine):
		global newLine, observables, obsFunctions, parameters
		if key == 1: #1: NAME  2: NUMBER  51: OP   4: NEWLINE  0: ENDMARKER
			var(value)
			if newLine == True:
				observables.append(parse_expr(value))
				obsFunctions.append(parse_expr(fullLine[(fullLine.find('=')+1):len(fullLine)]))
			else:
				parameters.append(parse_expr(value))
			newLine = False
		elif key == 4:
			newLine = True
	
	tokenize(read,useToken)

	parameters = sorted(list(set(parameters)), key=default_sort_key)

	for entry in observables:
		if entry in parameters:
			parameters.remove(entry)	
	
	return observables, obsFunctions, parameters

