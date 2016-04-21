from sympy import *
from numpy import concatenate
from sympy.parsing.sympy_parser import *
import csv

def SolveSymbLES(A,b):
    dim=shape(A)[0]
    Asave=A[:]
    Asave=Matrix(dim, dim, Asave)
    determinant=Asave.det()
    if(determinant==0):
        return([])
    result=[]
    for i in range(dim):
        A=Matrix(dim,dim,Asave)
        A.col_del(i)
        A=A.col_insert(i,b)
        result.append(simplify(A.det()/determinant))
    
    return(result)
    
def difftotal(expr, diffby, diffmap):
    #Take the total derivative with respect to a variable.

    #Example:

    #    theta, t, theta_dot = symbols("theta t theta_dot")
    #    difftotal(cos(theta), t, {theta: theta_dot})

    #returns

    #    -theta_dot*sin(theta)
    
    # Replace all symbols in the diffmap by a functional form
    fnexpr = expr.subs({s:s(diffby) for s in diffmap})
    # Do the differentiation
    diffexpr = diff(fnexpr, diffby)
    # Replace the Derivatives with the variables in diffmap
    derivmap = {Derivative(v(diffby), diffby):dv 
                for v,dv in diffmap.iteritems()}
    finaldiff = diffexpr.subs(derivmap)
    # Replace the functional forms with their original form
    return finaldiff.subs({s(diffby):s for s in diffmap})

def getIndOfParticipatingSpecies(SM, F, X, fastreact):
    liste=[]
    for f in F:
        for fr in fastreact: 
            if(fr in str(f)):
                testcol=SM.col(list(F).index(f))
                for i in range(len(testcol)):
                    if(testcol[i]!=0):
                        if(i not in liste):
                            liste.append(i)
    return(liste)
    
def getIndOfFastReactions(F, fastreact):    
    foundA=False
    for fr in fastreact:
        for i in range(len(F)):         
            if(fr in str(F[i])):
                if(foundA):
                    b=i
                else:
                    a=i
                    foundA=True
    if(foundA):
      return(a,b)
    else:
      print('Wrong rate constants specified!')
      return(0)

def QSS(filename,
        fastreact=[],
        state2Remove='A',
          SM=False,
          X=[],
          F=[],
          outputFormat='R'):
    if(filename==None):
      print('Use specified stoichiometry matrix ...')
    if(filename!=None):
        filename=str(filename)
        file=csv.reader(open(filename), delimiter=',')
        print('Reading csv-file ...')
        L=[]
        nrrow=0
        nrcol=0
        for row in file:
            nrrow=nrrow+1
            nrcol=len(row)
            L.append(row)
            
    #print("Test")
##### Define flux vector F
    if(filename!=None):
        F=[]    
        for i in range(1,len(L)):
            F.append(L[i][1])
            #print(F)
            F[i-1]=F[i-1].replace('^','**')
            F[i-1]=parse_expr(F[i-1])
        F=Matrix(F)
    else:
        #print(F)
        if(F!=[]):
            flist=[]
            for f in F:
                flist.append(parse_expr(f))
            F=Matrix(flist)
        else:
            print("You have to specify a flux vector or a model file!")
    #print(F)
##### Define state vector X
    if(filename!=None):
        X=[]
        X=L[0][2:]
        for i in range(len(X)):
            X[i]=parse_expr(X[i])               
        X=Matrix(X)
    else:
        if(X!=[]):
            xlist=[]
            for x in X:
                xlist.append(parse_expr(x))
            X=Matrix(xlist)
        else:
            print("You have to specify a state vector or a model file!")
    #print(X)    
##### Define stoichiometry matrix SM
    #print(SM)
    if(filename!=None):
        SM=[]
        for i in range(len(L)-1):
          SM.append(L[i+1][2:])        
        for i in range(len(SM)):
        	for j in range(len(SM[0])):
        		if (SM[i][j]==''):
        			SM[i][j]='0'
        		SM[i][j]=parse_expr(SM[i][j])    
        SM=Matrix(SM)
        SM=SM.T
    else:
        if(SM):
            SMfile=csv.reader(open("smatrix.csv"), delimiter=',')
            nrrow=0
            nrcol=0
            L=[]
            for row in SMfile:
                nrrow=nrrow+1
                nrcol=len(row)
                L.append(row)
            
            SM=[]
            for i in range(len(L)-1):
                SM.append(L[i+1][:])
            for i in range(len(SM)):
              for j in range(len(SM[0])):
        		    if (SM[i][j]=='NA'):
        			    SM[i][j]='0'
        		    SM[i][j]=parse_expr(SM[i][j])
            SM=Matrix(SM)
            SM=SM.T
        else:
            print("You have to specify a stoichiometry matrix or a model file.")
    #print(SM)     
    print('Simplifying System ...')
    PS=getIndOfParticipatingSpecies(SM, F, X, fastreact)
    i1,i2 = getIndOfFastReactions(F, fastreact)
    fastreactsymb=[parse_expr(fastreact[0]),parse_expr(fastreact[1])]
    mapping={}
    variables=[parse_expr('fastflux')]
    for el in X:
        if(list(X).index(el) in PS):
            mapping[el]=parse_expr(str(el)+'_dot')
            variables.append(parse_expr(str(el)+'_dot'))
    fastEq=F[i1]-F[i2]
    t = symbols("t")
    fastEqDiff=difftotal(F[i1], t, mapping)-difftotal(F[i2], t, mapping)
    #F.row_del(i1)
    #F=F.row_insert(i1,Matrix([parse_expr('fastflux')]))                   
    #F.row_del(i2)
    #F=F.row_insert(i2,Matrix([0]))
    F[i1]=parse_expr('fastflux')
    F[i2]=0
    eqs=[]
    for ps in PS:
        eqs.append((SM.row(ps)*F)[0]-parse_expr(str(X[ps])+'_dot'))
        #eqs.append(SM[ps]*F-parse_expr(str(X[ps])+'_dot'))
    eqs.append(fastEqDiff)
    #CM=getCM(eqs, variables)
    #print("Test")
    sol=solve(eqs, variables)
    if(sol!=[]):
      fastEqVar=parse_expr(state2Remove)
      fastEqExpr=solve(fastEq,fastEqVar)
      if(fastEqExpr!=[]):
        fastEqExpr=fastEqExpr[0]
      else:
        print("Did not find a solution for the fast reaction.")
        return([])
    else:
        print("Did not find a solution for the equation system.")
        return([])
    ausgabe=[]
    for var in variables:
        term=sol[var].subs(fastEqVar, fastEqExpr)
        term=term.subs(fastreactsymb[0],fastreactsymb[1]*parse_expr('ratio_'+str(fastreact[0])+'_'+str(fastreact[1])))
        term=simplify(term)
        ausgabe.append(str(var)+' = '+str(term))
        #print(str(var)+' = '+str(term))
    F[i1]=parse_expr(ausgabe[0].split(' = ')[1])
    fastEqExpr=fastEqExpr.subs(fastreactsymb[0],fastreactsymb[1]*parse_expr('ratio_'+str(fastreact[0])+'_'+str(fastreact[1])))
    fastEqExpr=simplify(fastEqExpr)
    print('Use the following replacement for state '+state2Remove+'!')
    print('    '+state2Remove+'='+str(fastEqExpr))
    return([state2Remove+'='+str(fastEqExpr),'flux'+str(i1)+'='+str(F[i1]),'flux'+str(i2)+'='+str(F[i2])])
