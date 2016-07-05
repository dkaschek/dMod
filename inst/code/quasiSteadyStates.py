from sympy import *
from numpy import concatenate
from numpy.linalg import matrix_rank
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
    liste=[]
    for fr in fastreact:
        for i in range(len(F)):         
            if(fr in str(F[i]) and i not in liste):
                liste.append(i)
    return(liste)
    
def findMin(vec):
    m=max(max(vec),max(-vec))    
    for el in vec:
        if(abs(el)<m and abs(el)>1e-8):
            m=abs(el)
    return(m)
    
def nullZ(A):
  ret = A.rref() # compute reduced row echelon form of A
  R = ret[0]  # matrix A in rref 
  pivcol = ret[1] #columns in which a pivot was found

  n = len(A.row(0)) # number of columns of A
  r = len(pivcol) # rank of reduced row echelon form
  nopiv = range(n)
  nopiv2 = [nopiv[i] for i in range(n) if i not in pivcol]  # columns in which no pivot was found
#  print(ret)
#  print(nopiv2)
#  print(n)
#  print(r)
  if(n > r):
    Z=eye(n-r)
    if(r>0):
        Z=concatenate((-R[pivcol, nopiv2], Z), axis=0)
  return(Z) 

def QSS(filename,
        fastreact=[],
        state2Remove=[],
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
    index_list = getIndOfFastReactions(F, fastreact)
    frsymb_list=[parse_expr(fastreact[i]) for i in range(len(fastreact))]
    mapping={}
    #F_list=[parse_expr('fastflux'+str(i)) for i in range(len(PS))]
    #variables=F_list[:]
    variables=[]
    for el in X:
        if(list(X).index(el) in PS):
            mapping[el]=parse_expr(str(el)+'_dot')
            #variables.append(parse_expr(str(el)+'_dot'))
    #print(index_list)
    FF=Matrix([F[i] for i in index_list])
    SMF=SM[PS,index_list]
    SMStimesFS=SM[PS,:]*F-SMF*FF
    #print(SMF)
    #print(SMStimesFS)
    for i in index_list:
        F[i]=parse_expr('F_'+str(i))
        #variables.append(parse_expr('F_'+str(i)))
    #F_red2=Matrix([F[i] for i in index_list])
    eqs=[]
    #print(PS)
    for ps in PS:
        #print(SMStimesFS[list(PS).index(ps)])
        eqs.append(parse_expr(str(X[ps])+'_dot')-parse_expr(str(X[ps])+'_tilde')-SMStimesFS[list(PS).index(ps)])
        variables.append(parse_expr(str(X[ps])+'_dot'))
        #eqs.append((SM_red*F_red2)[PS.index(ps)]-parse_expr('G_'+str(X[ps])))
        variables.append(parse_expr(str(X[ps])+'_tilde'))
    if(max(SMF.shape) > matrix_rank(SMF)):
        ns=nullZ(SMF.T)
        #print(ns)
        for i in range(ns.shape[1]):
            eq=0
            factor=1/findMin(ns[:,i])
            for ps in PS:
                eq=eq+ns[:,i][list(PS).index(ps)]*factor*parse_expr(str(X[ps])+'_tilde')
            eqs.append(eq)
    #print(ns)        
    t = symbols("t")
    fastEqDiff_list=[difftotal((SMF*FF)[i], t, mapping) for i in range(len(PS))]
    eqs=eqs+fastEqDiff_list
    #print(eqs)
    #print(variables)
    sol=solve(eqs, variables)
    if(sol==[]):
        print("Did not find a solution for the equation system.")
        return([])    
    #print(frsymb_list)
    #print(state2Remove)
    if(state2Remove==[]):
        varfast=[X[ps] for ps in PS]
        pivcol = SMF.rref()[1] #columns in which a pivot was found
        varfast=[varfast[i] for i in pivcol]
        state2Remove=[str(v) for v in varfast]
    else:
        if(matrix_rank(SMF)==len(state2Remove)):
            varfast=[parse_expr(state) for state in state2Remove]
        else:
            print("Rank of the fast stoichiometry matrix equals {}. Please specify {} states to remove from the system!".format(matrix_rank(SMF), matrix_rank(SMF)) )
            return([])  
    #print(varfast)
    #print(SMF*FF)
    solfast=solve(SMF*FF,varfast)
    #print(isinstance(solfast, list))
    #print(isinstance(solfast[0], list))
    if(not isinstance(solfast, list)):
        varfast=solfast.keys()
        solfast=[solfast[el] for el in solfast.keys()]
        #print(isinstance(solfast, list))
    else:
        if(isinstance(solfast[0], tuple)):
            liste=[]
            for i in range(len(solfast[0])):
                liste.append(solfast[0][i])
            solfast=liste
    #print(varfast)
    #print(solfast)
    #print(solfast[0])
    ausgabe=[]
    for var in variables:
        if('tilde' not in str(var) and str(var).split('_dot')[0] not in state2Remove):
            term=sol[var]
            for el in range(len(varfast)):
                term=term.subs(varfast[el], solfast[el])
            for i in range(1,len(frsymb_list)):
                term=term.subs(frsymb_list[i],frsymb_list[0]*parse_expr('r_'+str(frsymb_list[i])+'_'+str(frsymb_list[0])))
            term=simplify(term)        
            ausgabe.append(str(var)+' = '+str(term))
    
    #for ps in PS:
    #    if(str(X[ps])!=statenot2Remove):
    print('Use the following observation functions!')
    #print(solfast[0])
    for el in range(len(varfast)):
        term=solfast[el]
        for i in range(1,len(frsymb_list)):
            term=term.subs(frsymb_list[i],frsymb_list[0]*parse_expr('r_'+str(frsymb_list[i])+'_'+str(frsymb_list[0])))
        term=simplify(term)
        print('   '+str(varfast[el])+' = '+str(term))
    #print(solfast)
    #print(len(solfast))
    #print(ausgabe)
    #print(varfast)
    for i in range(len(varfast)):
        ausgabe.append(str(varfast[i]))
    ausgabe.append(str(len(varfast)))
    #print(ausgabe.append(varfast))
    print("Done")
    return(ausgabe)
