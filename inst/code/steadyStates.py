import sympy
from sympy import *
import numpy
from numpy import *
from numpy.linalg import matrix_rank
from sympy.parsing.sympy_parser import *
from sympy.matrices import *
import csv
import os
import random
from random import shuffle

def Reduce(eq):
    for el in ['kPDK1', 'kAkt']:
        el=parse_expr(el)
        if(eq.subs(el,0)==0):
            eq=expand(simplify(eq/el))
    return expand(eq)


def LCS(s1, s2):
    m = [[0] * (1 + len(s2)) for i in xrange(1 + len(s1))]
    longest, x_longest = 0, 0
    for x in xrange(1, 1 + len(s1)):
        for y in xrange(1, 1 + len(s2)):
            if s1[x - 1] == s2[y - 1]:
                m[x][y] = m[x - 1][y - 1] + 1
                if m[x][y] > longest:
                    longest = m[x][y]
                    x_longest = x
            else:
                m[x][y] = 0
    return s1[x_longest - longest: x_longest]

def SolveSymbLES(A,b):
    dim=shape(A)[0]
    Asave=A[:]
    Asave=Matrix(dim, dim, Asave)
    #printmatrix(Asave)
    #print(b)
    determinant=Asave.det()
    if(determinant==0):
        #print('Determinant of LCL-calculation is zero! Try to specify LCLs yourself!')
        return([])
    result=[]
    for i in range(dim):
        A=Matrix(dim,dim,Asave)
        A.col_del(i)
        A=A.col_insert(i,b)
        result.append(simplify(A.det()/determinant))
    
    return(result)

def CutStringListatSymbol(liste, symbol):
    out=[]    
    for el in liste:
        if(symbol in el):
            add=el.split(symbol)
        else:
            add=[el]
        out=out+add
    return(out)

def FillwithRanNum(M):
    dimx=len(M.row(0))
    dimy=len(M.col(0))
    ranM=zeros(dimy, dimx)
    parlist=[]
    ranlist=[]
    for i in M[:]:
        if(i!=0):
            if(str(i)[0]=='-'):
                parlist.append(str(i)[1:])
            else:
                parlist.append(str(i))
    parlist=list(set(parlist))
    for symbol in [' - ', ' + ', '*', '/', '(',')']:
        parlist=CutStringListatSymbol(parlist,symbol)
    parlist=list(set(parlist))
    temp=[]    
    for i in parlist:
        if(i!=''):
            if(not is_number(i)):
                temp.append(i)
                ranlist.append(random.random())
    parlist=temp
    for i in range(dimy):
        for j in range(dimx):
            ranM[i,j]=M[i,j]
            if(ranM[i,j]!=0):
                for p in range(len(parlist)):
                   ranM[i,j]=ranM[i,j].subs(parse_expr(parlist[p]),ranlist[p])
    return(ranM)

def FindLinDep(M, tol=1e-12):
    ranM=FillwithRanNum(M)
    Q,R=numpy.linalg.qr(ranM)
    for i in range(shape(R)[0]):
        for j in range(shape(R)[1]):
            if(abs(R[i,j]) < tol):
                R[i,j]=0.0
                
    LinDepList=[]
    for i in range(shape(R)[0]):
        if(R[i][i]==0):
            LinDepList.append(i)
    
    return(LinDepList)

def FindLCL(M, X):
    LCL=[]    
    LinDepList=FindLinDep(M)
    i=0
    counter=0
    deleted_rows=[]
    states=Matrix(X[:])
    while(LinDepList!=[]):
        i=LinDepList[0]
        testM=FillwithRanNum(M)
        rowliste=list(numpy.nonzero(testM[:,i])[0])
        colliste=[i]        
        for z in range(i):
            for k in rowliste:        
                for j in range(i):
                    jliste=list(numpy.nonzero(testM[:,j])[0])
                    if(k in jliste):
                        rowliste=rowliste+jliste
                        colliste=colliste+[j]
            rowliste=list(set(rowliste))
            colliste=list(set(colliste))
        rowliste.sort()
        colliste.sort()
        colliste.pop()        
        rowlisteTry=rowliste[0:(len(colliste))]
        vec=SolveSymbLES(M[rowlisteTry,colliste],M[rowlisteTry,i])
        shufflecounter=0
        while(vec==[] and shufflecounter < 100):
            shuffle(rowliste)
            shufflecounter=shufflecounter+1
            rowlisteTry=rowliste[0:(len(colliste))]
            vec=SolveSymbLES(M[rowlisteTry,colliste],M[rowlisteTry,i])
        if(shufflecounter==100):
            print('Problems while finding conserved quantities!')
            return(0,0)
        counter=counter+1
        try:
            mat=[states[l] for l in colliste]
            test=parse_expr('0')
            for v in range(0,len(vec)):
                test=test-parse_expr(str(vec[v]))*parse_expr(str(mat[v]))
        except:
            return([],0)
        partStr=str(test)+' + '+str(states[i])
        partStr=partStr.split(' + ')
        partStr2=[]
        for index in range(len(partStr)):
            partStr2=partStr2+partStr[index].split('-')
        partStr=partStr2
        if(len(partStr) > 1):        
            CLString=LCS(str(partStr[0]),str(partStr[1]))
            for ps in range(2,len(partStr)):
                CLString=LCS(CLString,str(partStr[ps]))
        else:
            CLString=str(partStr[0])
        if(CLString==''):
            CLString=str(counter)
        LCL.append(str(test)+' + '+str(states[i])+' = '+'total'+CLString)
        M.col_del(i)
        states.row_del(i)
        deleted_rows.append(i+counter-1)
        LinDepList=FindLinDep(M)
    return(LCL, deleted_rows)

def printmatrix(M):    
    lengths=[]
    for i in range(len(M.row(0))):
        lengths.append(0)
        for j in range(len(M.col(0))):
            lengths[i]=max(lengths[i],len(str(M.col(i)[j])))          
    string=''.ljust(5)
    string2=''.ljust(5)
    for j in range(len(M.row(0))):
        string=string+(str(j)).ljust(lengths[j]+2)
        for k in range(lengths[j]+2):        
            string2=string2+('-')        
    print(string)
    print(string2)
    for i in range(len(M.col(0))):
        string=str(i).ljust(4) + '['
        for j in range(len(M.row(0))):
            if(j==len(M.row(0))-1):
                string=string+str(M.row(i)[j]).ljust(lengths[j])
            else:
                string=string+(str(M.row(i)[j])+', ').ljust(lengths[j]+2)        
        print(string+']')    
    return()
    
def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False
    
def checkNegRows(M):
    NegRows=[]
    if((M==Matrix(0,0,[])) | (M==Matrix(0,1,[])) | (M==Matrix(1,0,[]))):
        return(NegRows)
    else:        
        for i in range(len(M.col(0))):
            foundPos=False
            for j in range(len(M.row(i))):
                if(M[i,j]>0):
                    foundPos=True
            if(foundPos==False):
                NegRows.append(i)    
        return(NegRows)
    
def checkPosRows(M):
    PosRows=[]
    if((M==Matrix(0,0,[])) | (M==Matrix(0,1,[])) | (M==Matrix(1,0,[]))):
        return(PosRows)
    else: 
        for i in range(len(M.col(0))):
            foundNeg=False
            for j in range(len(M.row(i))):
                if(M[i,j]<0):
                    foundNeg=True
            if(foundNeg==False):
                PosRows.append(i)    
        return(PosRows)             
    
def find_cycle(graph, start, end, path=[]):
    path = path + [start]
    if not graph.has_key(start):
        return None
    if ((start == end) & (path!=[start])):
        return path    
    for node in graph[start]:
        if node==end: return (path+[end])
        if node not in path:
            #print(node)
            newpath = find_cycle(graph, node, end, path)
            if newpath: return newpath
    return None

def ChooseOptimalFlux(eq, SM, F, compromise1, compromise2):    
    anzneg=str(eq.args).count('-')
    anzpos=len(eq.args)-anzneg
    #print(anzneg)
    #print(anzpos)
    #print(eq)
    eq=expand(eq)
    if(anzneg==1):
        for arg in eq.args:
            counter=0
            if((str(arg)[0:2]=='2*') or (str(arg)[0:3]=='-2*')):
                    arg=arg/2
            for i in range(len(SM.col(0))):
                if(SM[i,list(F).index(parse_expr(str(arg).replace('-','')))]!=0):
                    counter=counter+1
            #print(counter)
            if(('-' in str(arg)) & (counter==1)):
                return(arg)
    if(anzpos==1):
        for arg in eq.args:
            #print(arg)
            counter=0
            if((str(arg)[0:2]=='2*') or (str(arg)[0:3]=='-2*')):
                    arg=arg/2
            #printmatrix(SM)
            #print(F)
            for i in range(len(SM.col(0))):
                if(SM[i,list(F).index(parse_expr(str(arg).replace('-','')))]!=0):
                    counter=counter+1
            #print(counter)
            if(('-' not in str(arg)) & (counter==1)):
                return(arg)
    if(compromise1):
        for arg in eq.args:
            counter=0
            #print(arg)
            for i in range(len(SM.col(0))):
                if(SM[i,list(F).index(parse_expr(str(arg).replace('-','')))]!=0):
                    counter=counter+1
            if(counter==1):
                return(arg)
    if(compromise1 & compromise2):
        return(eq.args[0])
    return None

def FindNodeToSolve(graph):
    for el in graph:
        if(graph[el]==[]):
            return(el)
    return(None)

def CountNZE(V):
    counter=0
    for v in V:
        if(v!=0):
            counter=counter+1
    return(counter)
    
def Sparsify(M, level):
    if(level==3):
        ncol=len(M.row(0))
        for i in range(ncol):
            tobeat=CountNZE(M.col(i))
            for j in range(ncol):
                if(i!=j):
                    for factor_j in [1,2,-1,-2,0]:
                        for k in range(ncol):
                            if(k!=i and k!=j):
                                for factor_k in [1,2,-1,-2,0]:
                                    for l in range(ncol):
                                        if(l!=i and l!=j and l!=k):
                                            for factor_l in [1,2,-1,-2,0]:
                                                test=M.col(i)+factor_j*M.col(j)+factor_k*M.col(k)+factor_l*M.col(l)
                                                if(tobeat > CountNZE(test)):
                                                    M.col_del(i)
                                                    M=M.col_insert(i,test)
                                                    tobeat=CountNZE(test)
            print(str(i+1)+' columns of '+str(ncol) +' done')
    if(level==2):
        ncol=len(M.row(0))
        for i in range(ncol):
            tobeat=CountNZE(M.col(i))
            for j in range(ncol):
                if(i!=j):
                    for factor_j in [1,2,-1,-2,0]:
                        for k in range(ncol):
                            if(k!=i and k!=j):
                                for factor_k in [1,2,-1,-2,0]:
                                    test=M.col(i)+factor_j*M.col(j)+factor_k*M.col(k)
                                    if(tobeat > CountNZE(test)):
                                        M.col_del(i)
                                        M=M.col_insert(i,test)
                                        tobeat=CountNZE(test)
            print(str(i+1)+' columns of '+str(ncol) +' done')
    if(level==1):
        ncol=len(M.row(0))
        for i in range(ncol):
            tobeat=CountNZE(M.col(i))
            for j in range(ncol):
                if(i!=j):
                    for factor_j in [1,2,-1,-2,0]:
                        test=M.col(i)+factor_j*M.col(j)
                        if(tobeat > CountNZE(test)):
                            M.col_del(i)
                            M=M.col_insert(i,test)
                            tobeat=CountNZE(test)
    return(M)
    
def SuggestTrafos(eqOut, UsedVars):
    out=[]
    for i in range(len(eqOut)):
        eq=eqOut[i]
        #print(eq)
        if('-' in eq):
            print('Finding trafo...')
            foundTrafo=False
            rs=parse_expr(eq.split(' = ')[1])
            if(('1/' in str(rs.args[0])) & (len(rs.args)>=2)):
                args=(rs.args[-1]).args
            else:
                args=rs.args
            #print(args)
            anzneg=str(args).count('-')
            anzpos=len(args)-anzneg
            #print(anzpos)
            if(anzneg==1):
                for arg in args:
                    if('-' in str(arg)):
                        for var in arg.atoms():
                            if((not is_number(str(var))) & (var not in UsedVars) & (not foundTrafo)):
                                sol=solve(rs, var, simplify=False)[0]
                                foundTrafo=True
                                trafo=parse_expr('('+str(sol)+')/(1+r_'+str(var)+')')
                                trafoVar=var
                                out.append(str(trafoVar)+' = '+str(trafo))
            if((anzpos==1) & (not foundTrafo)):
                for arg in args:
                    if('-' not in str(arg)):
                        for var in arg.atoms():
                            if((not is_number(str(var))) & (var not in UsedVars) & (not foundTrafo)):
                                sol=solve(rs, var, simplify=False)[0]
                                foundTrafo=True
                                trafo=parse_expr('('+str(sol)+')*(1+r_'+str(var)+')')
                                trafoVar=var
                                out.append(str(trafoVar)+' = '+str(trafo))
            
            if((anzpos > 1) & (anzneg > 1)):
                negliste=[]
                posliste=[]
                posges=''
                anzmal=0
                for arg in args:
                    if('-' in str(arg)):
                        negliste.append(arg)
                    else:
                        posges=posges+str(arg)
                        posliste.append(arg)
                        anzmal=anzmal+str(arg).count('*')
                if(str(simplify(parse_expr(posges))).count('*')!=anzmal):
                    var=parse_expr(posges.split('*')[0])
                    if(var not in UsedVars):
                        sol=solve(rs, var)[0]
                        foundTrafo=True
                        trafo=parse_expr('('+str(sol)+')*(1+r_'+str(var)+')')
                        trafoVar=var
                        out.append(str(trafoVar)+' = '+str(trafo))
                        
            if(not foundTrafo):
                print('Alles Mist')
            else:
                #print(trafoVar)
                #print(out[-1])
                for j in range(len(eqOut)):
                    if(j >= i):
                        ls,rs=eqOut[j].split(' = ')
                        eqOut[j]=ls+' = '+str(simplify(parse_expr(rs).subs(trafoVar, trafo)))
            
        out.append(eqOut[i])
        #print(out)
    return(out)

def ODESS(filename,
          injections=[],
          forbidden=[],
          ensurePositivity=True,
          sparsifyLevel = 2,
          outputFormat='R'):
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
        
    nrspecies=nrcol-2
    
##### Remove injections  
    counter=0
    for i in range(1,len(L)):
        if(L[i-counter][1] in injections):
            L.remove(L[i-counter])
            counter=counter+1       
    
##### Define flux vector F	
    F=[]
    
    for i in range(1,len(L)):
        F.append(L[i][1])
        #print(F)
        F[i-1]=F[i-1].replace('^','**')
        F[i-1]=parse_expr(F[i-1])
        for inj in injections:
            F[i-1]=F[i-1].subs(parse_expr(inj),0)
    F=Matrix(F)
    #print(F)
##### Define state vector X
    X=[]
    X=L[0][2:]
    for i in range(len(X)):
        X[i]=parse_expr(X[i])               
    X=Matrix(X)
    #print(X)
    Xo=X.copy()
        
##### Define stoichiometry matrix SM
    SM=[]
    for i in range(len(L)-1):
    	SM.append(L[i+1][2:])        
    for i in range(len(SM)):
    	for j in range(len(SM[0])):
    		if (SM[i][j]==''):
    			SM[i][j]='0'
    		SM[i][j]=parse_expr(SM[i][j])    
    SM=Matrix(SM)
    
    
##### Increase Sparsity of stoichiometry matrix SM
    print('Sparsify Stoichiometry Matrix with sparsify-level '+str(sparsifyLevel)+'!')
    SMorig=(SM.copy()).T
    SM=Sparsify(SM, level=sparsifyLevel)
    SM=SM.T
    #return(SM)
    
##### Read forbidden rates
    UsedRC=[]
    for el in forbidden:
        UsedRC.append(parse_expr(el))
    
##### Check for zero fluxes
    icounter=0
    jcounter=0
    for i in range(len(F)):
        if(F[i-icounter]==0):
            F.row_del(i-icounter)
            for j in range(len(SM.col(i-icounter))):
                if(SM[j-jcounter,i-icounter]!=0):
                    #UsedRC.append(X[j-jcounter])
                    X.row_del(j-jcounter)
                    SM.row_del(j-jcounter)
                    SMorig.row_del(j-jcounter)
                    jcounter=jcounter+1
            SM.col_del(i-icounter)
            SMorig.col_del(i-icounter)
            icounter=icounter+1
    
    print('Removed '+str(icounter)+' fluxes that are a priori zero!')
    nrspecies=nrspecies-icounter
    #printmatrix(SM)
    #print(F)
    #print(X)
    #print(UsedRC)
#####Check if some species are zero and remove them from the system
    zeroStates=[]
    NegRows=checkNegRows(SM)
    PosRows=checkPosRows(SM)
    #print(PosRows)
    #print(NegRows)
    while((NegRows!=[]) | (PosRows!=[])):
        #print(PosRows)
        #print(NegRows)
        if(NegRows!=[]):        
            row=NegRows[0]
            zeroStates.append(X[row])
            counter=0    
            for i in range(len(F)):
                if(F[i-counter].subs(X[row],1)!=F[i-counter] and F[i-counter].subs(X[row],0)==0):
                    F.row_del(i-counter)
                    SM.col_del(i-counter)                    
                    SMorig.col_del(i-counter)
                    counter=counter+1
                else:
                    if(F[i-counter].subs(X[row],1)!=F[i-counter] and F[i-counter].subs(X[row],0)!=0):
                        F[i-counter]=F[i-counter].subs(X[row],0)
            X.row_del(row)
            SM.row_del(row)
            SMorig.row_del(row)
        else:
            row=PosRows[0]
            zeroFluxes=[]
            for j in range(len(SM.row(row))):
                if(SM.row(row)[j]!=0):
                    zeroFluxes.append(F[j])
            for k in zeroFluxes:
                StateinFlux=[]
                for state in X:
                    if(k.subs(state,1)!=k):
                        StateinFlux.append(state)
                if(len(StateinFlux)==1):
                    zeroStates.append(StateinFlux[0])
                    row=list(X).index(StateinFlux[0])
                    counter=0            
                    for i in range(len(F)):
                        if(F[i-counter].subs(X[row],1)!=F[i-counter]):
                            if(F[i-counter].subs(X[row],0)==0):
                                F.row_del(i-counter)
                                SM.col_del(i-counter)
                                SMorig.col_del(i-counter)
                            else:
                                F[i-counter]=F[i-counter].subs(X[row],0)                            
                            counter=counter+1
        #printmatrix(SM)
        NegRows=checkNegRows(SM)      
        PosRows=checkPosRows(SM)
    #printmatrix(SM)
    #print(F)
    #print(X)
    nrspecies=nrspecies-len(zeroStates)
    if(nrspecies==0):
        print('All states are zero!')
        return(0)
    else:
        if(zeroStates==[]):
            print('No states found that are a priori zero!')
        else:
            print('These states are zero:')
            for state in zeroStates:
                print('\t'+str(state))
    
    nrspecies=nrspecies+len(zeroStates)

##### Identify linearities, bilinearities and multilinearities        
    Xsquared=[]
    for i in range(len(X)):
        Xsquared.append(X[i]*X[i])        
    Xsquared=Matrix(Xsquared)
      
    BLList=[]
    MLList=[]
    for i in range(len(SM*F)):
        LHS=str(expand((SM*F)[i]))
        LHS=LHS.replace(' ','')
        LHS=LHS.replace('-','+')
        LHS=LHS.replace('**2','tothepowerof2')
        LHS=LHS.replace('**3','tothepowerof3')
        exprList=LHS.split('+')
        for expr in exprList:
            VarList=expr.split('*')
            counter=0
            factors=[]
            for j in range(len(X)):
                anz=0
                if(str(X[j]) in VarList):
                    anz=1
                    factors.append(X[j])
                if((str(X[j])+'tothepowerof2') in VarList):
                    anz=2 
                    factors.append(X[j])
                    factors.append(X[j])
                if((str(X[j])+'tothepowerof3') in VarList):
                    anz=3
                    factors.append(X[j])
                    factors.append(X[j])
                    factors.append(X[j])
                counter=counter+anz
            if(counter==2):
                string=''            
                for l in range(len(factors)):
                    if(l==len(factors)-1):
                        string=string+str(factors[l])
                    else:
                        string=string+str(factors[l])+'*'
                if(not(string in BLList)):
                    BLList.append(string)
            if(counter>2):
                string=''            
                for l in range(len(factors)):
                    if(l==len(factors)-1):
                        string=string+str(factors[l])
                    else:
                        string=string+str(factors[l])+'*'
                if(not(string in MLList)):
                    MLList.append(string)
        
    COPlusLIPlusBL=[]
    for i in range(len(SM*F)):
        COPlusLIPlusBL.append((SM*F)[i])
        for j in range(len(MLList)):
            ToSubs=expand((SM*F)[i]).coeff(MLList[j])
            COPlusLIPlusBL[i]=expand(COPlusLIPlusBL[i]-ToSubs*parse_expr(MLList[j]))
            
    COPlusLI=[]
    for i in range(len(COPlusLIPlusBL)):
        COPlusLI.append(COPlusLIPlusBL[i])
        for j in range(len(BLList)):
            ToSubs=expand((COPlusLIPlusBL)[i]).coeff(BLList[j])
            COPlusLI[i]=expand(COPlusLI[i]-ToSubs*parse_expr(BLList[j]))
    
##### C*X contains linear terms
    C=zeros(len(COPlusLI),len(X))  
    for i in range(len(COPlusLI)):
    	for j in range(len(X)):
    		C[i*len(X)+j]=expand((COPlusLI)[i]).coeff(X[j])
        
##### ML contains multilinearities
    ML=expand(Matrix(SM*F)-Matrix(COPlusLIPlusBL))
##### BL contains bilinearities
    BL=expand(Matrix(COPlusLIPlusBL)-Matrix(COPlusLI))    
#### CM is coefficient matrix of linearities
    CM=C        
#####CMBL gives coefficient matrix of bilinearities
    CMBL=[]
    if(BLList!=[]):
        for i in range(len(BLList)):
            CVBL=[]
            for k in range(len(BL)):
                CVBL.append(BL[k].coeff(BLList[i]))
            CMBL.append(CVBL)            
    else:
        CVBL=[]
        for k in range(len(BL)):
            CVBL.append(0)
        CMBL.append(CVBL)
    
    CMBL=Matrix(CMBL).T 
    
#####CMML gives coefficient matrix of multilinearities
#####Summarize multilinearities and bilinearities 
    if(MLList!=[]):
        CMML=[]
        for i in range(len(MLList)):
            CVML=[]
            for k in range(len(ML)):
                CVML.append(expand(ML[k]).coeff(MLList[i]))
            CMML.append(CVML)    
        CMML=Matrix(CMML).T  
        BLList=BLList+MLList
        CMBL=Matrix(concatenate((CMBL,CMML),axis=1))
      
    for i in range(len(BLList)):
        BLList[i]=parse_expr(BLList[i])
       
    if(BLList!=[]):    
        CMbig=Matrix(concatenate((CM,CMBL),axis=1))
    else:
        CMbig=Matrix(CM)      

#### Save ODE equations for testing solutions at the end    
    print('Rank of SM is '+str(SM.rank()) + '!')
    ODE=SMorig*F
#### Find conserved quantities
    print('\nFinding conserved quantities ...')
    #printmatrix(CMbig)
    #print(X)
    LCLs, rowsToDel=FindLCL(CMbig.transpose(), X)
    if(LCLs!=[]):
        print(LCLs)
    else:
        print('System has no conserved quantities!')
#### Define graph structure
    print('\nDefine graph structure ...\n')
    graph={}    
    for i in range(len(SM*F)):
        liste=[]
        for j in range(len(X)):
            if((SM*F)[i]!=((SM*F)[i]).subs(X[j],1)):
                if(j==i):
                    In=((SM*F)[i]).subs(X[j],0)
                    Out=simplify(((SM*F)[i]-In)/X[j])
                    if(Out!=Out.subs(X[j],1)):
                        liste.append(str(X[j]))
                else:
                    liste.append(str(X[j]))
        graph[str(X[i])]=liste
        
#### Remove cycles step by step    
    gesnew=0
    noCycle=False
    #UsedRC=[]
    UsedStates=[]
    eqOut=[]
    DependentRates=[]
    DependentOn=[]
    counter=0
    while(not noCycle):    
        i=0
        changeposneg=False
        counter=counter+1
        foundCycle=False 
        #print(graph)
        while((not foundCycle) & (not noCycle)):
            cycle=find_cycle(graph, str(X[i]), str(X[i]))
            if cycle is None:
                i=i+1
                if(i==len(X)):
                    noCycle=True
                    print('There is no cycle in the system!\n')
            else:
                print('Removing cycle '+str(counter))
                foundCycle=True
        CycleCanBeDoneByCL=False
        if(not noCycle):
            print(cycle)
            for node in cycle:
                for LCL in LCLs:
                    if(not CycleCanBeDoneByCL):
                        ls=parse_expr(LCL.split(' = ')[0])
                        if(ls.subs(parse_expr(node),1)!=ls):
                            CycleCanBeDoneByCL=True
                            LCLToRemove=LCL
                            UsedStates.append(parse_expr(node))
                            nodeToRep=node
                            print('   '+str(nodeToRep)+' --> '+'Done by CL')
        
        if(CycleCanBeDoneByCL):
            LCLs.remove(LCLToRemove)
            indexToRep=list(X).index(parse_expr(nodeToRep))
            eqOut.append(str(nodeToRep)+' = '+str(nodeToRep))
            X.row_del(indexToRep)
            SM.row_del(indexToRep)
                                
        if((not noCycle) & (not CycleCanBeDoneByCL)):
            compromise=False
            compromise2=False
            foundOptNodeAndFlux=False
            for node in cycle:
                if((not foundOptNodeAndFlux)  & (parse_expr(node) in list(X))):
                    nodeToRep=node
                    eq=Reduce((SM*F)[list(X).index(parse_expr(nodeToRep))])
                    FluxToTake=ChooseOptimalFlux(eq,SM,F, compromise1=False, compromise2=False) 
                    
                    if(not FluxToTake is None):
                        foundOptNodeAndFlux=True
                    #print(foundOptNodeAndFlux)
            if(not foundOptNodeAndFlux):
                eq=Reduce(expand((SM*F)[i]))
                l0=len(str(eq))
                nodeToRep=str(X[i])
                for node in cycle:
                            l=len(str((SM*F)[list(X).index(parse_expr(node))]))
                            if(l<l0):
                                eq=Reduce((SM*F)[list(X).index(parse_expr(node))])
                                nodeToRep=node                
                print('   Do compromise for: '+nodeToRep)
                print(eq)
                compromise=True
                trafoList=[]
                negs=parse_expr('0')
                poss=parse_expr('0')
                anzneg=0
                anzpos=0
                for arg in eq.args:
                    if('-' in str(arg)):
                        anzneg=anzneg+1
                        negs=negs-arg
                    else:
                        anzpos=anzpos+1
                        poss=poss+arg
                if(anzpos == 1 and anzneg == 1):
                    trafoList.append(str(poss)+'='+str(negs))
                else:
                    if(anzneg==1):
                        nenner=1
                        gesnew=gesnew+(anzneg-1)
                        for j in range(anzneg):
                            if(j>0):
                                nenner=nenner+parse_expr('r_'+nodeToRep+'_'+str(j))
                        trafoList.append(str(negs)+'=('+str(poss)+')*1/('+str(nenner)+')')
                    else:
                        if(anzpos==1):
                            nenner=1
                            gesnew=gesnew+(anzpos-1)
                            for j in range(anzpos):
                                if(j>0):
                                    nenner=nenner+parse_expr('r_'+nodeToRep+'_'+str(j))
                            trafoList.append(str(poss)+'=('+str(negs)+')*1/('+str(nenner)+')')                    
                        else:
                            if((anzpos > anzneg) & (not changeposneg)):
                                nenner=1
                                gesnew=gesnew+(anzneg-1)
                                for j in range(anzneg):
                                    if(j>0):
                                        nenner=nenner+parse_expr('r_'+nodeToRep+'_'+str(j))
                                for j in range(anzneg):
                                    if(j==0):
                                        trafoList.append(str(negs.args[j])+'=('+str(poss)+')*1/('+str(nenner)+')')
                                    else:
                                        trafoList.append(str(negs.args[j])+'=('+str(poss)+')*'+'r_'+nodeToRep+'_'+str(j)+'/('+str(nenner)+')')
                            else:
                                nenner=1
                                #print(poss)
                                #print(poss.args)
                                gesnew=gesnew+(anzpos-1)
                                for j in range(anzpos):
                                    if(j>0):
                                        nenner=nenner+parse_expr('r_'+nodeToRep+'_'+str(j))
                                for j in range(anzpos):
                                    if(j==0):
                                        trafoList.append(str(poss.args[j])+'=('+str(negs)+')*1/('+str(nenner)+')')
                                    else:
                                        trafoList.append(str(poss.args[j])+'=('+str(negs)+')*'+'r_'+nodeToRep+'_'+str(j)+'/('+str(nenner)+')')
                         
                tempcounter=0
                #print(trafoList)
                for trafo in trafoList:
                    #print(trafo)
                    lstrafo,rstrafo=trafo.split('=')
                    lstrafo=parse_expr(lstrafo)
                    rstrafo=parse_expr(rstrafo)
                    foundRateConst=False
                    for var in list(lstrafo.atoms()):
                        if((var not in X) & (var not in UsedStates) & (var not in UsedRC) & (not is_number(str(var))) & (not foundRateConst)):
                            RC=var
                            UsedRC.append(RC)
                            foundRateConst=True
                    if(foundRateConst):
                        rest=lstrafo/RC
                        sol=rstrafo/rest
                        eqOut.append(str(RC)+' = '+str(sol))
                        DependentRates.append(RC)
                        DependentOn.append(list(sol.atoms()))
                        tempcounter=tempcounter+1
                    else:                        
                        compromise2=True
                if(compromise2):
                    print('Did not find appropriate transformation!')
                    for temp in range(tempcounter):
                        eqOut.pop()
                        DependentRates.pop()
                        DependentOn.pop()
                        UsedRC.pop() 
                    
                if(not compromise2):
                    indexToRep=list(X).index(parse_expr(nodeToRep))
                    X.row_del(indexToRep)
                    SM.row_del(indexToRep)
                    for temp in range(tempcounter):
                            eq=eqOut[-(temp+1)]
                            RC=parse_expr(eq.split('=')[0])
                            sol=parse_expr(eq.split('=')[1])
                            for f in range(len(F)):
                                F[f]=F[f].subs(RC, sol)
                            print('    '+str(RC)+' --> '+str(sol))
                leaveCompromiseLoop=False
                lcycle=len(cycle)
                cyclecounter=0
                while((leaveCompromiseLoop==False) & (cyclecounter < lcycle) & compromise2):                    
                    nodeToRep=cycle[cyclecounter]
                    eq=expand((SM*F)[list(X).index(parse_expr(nodeToRep))])
                    #print(eq)
                    cyclecounter=cyclecounter+1
                    #print(cycle)
                    #print(nodeToRep)
                    print(cyclecounter)
                    print('   Do compromise 2')
                    compromise2=False
                    trafoList=[]
                    negs=parse_expr('0')
                    poss=parse_expr('0')
                    anzneg=0
                    anzpos=0
                    for arg in eq.args:
                        if('-' in str(arg)):
                            anzneg=anzneg+1
                            negs=negs-arg
                        else:
                            anzpos=anzpos+1
                            poss=poss+arg
                            
                    if((anzpos > anzneg) & (not changeposneg)):
                        nenner=1
                        gesnew=gesnew+(anzneg-1)
                        for j in range(anzneg):
                            if(j>0):
                                nenner=nenner+parse_expr('r_'+nodeToRep+'_'+str(j))
                        for j in range(anzneg):
                            if(j==0):
                                trafoList.append(str(negs.args[j])+'=('+str(poss)+')*1/('+str(nenner)+')')
                            else:
                                trafoList.append(str(negs.args[j])+'=('+str(poss)+')*'+'r_'+nodeToRep+'_'+str(j)+'/('+str(nenner)+')')
                    else:
                        nenner=1
                        gesnew=gesnew+(anzpos-1)
                        for j in range(anzpos):
                            if(j>0):
                                nenner=nenner+parse_expr('r_'+nodeToRep+'_'+str(j))
                        for j in range(anzpos):
                            if(j==0):
                                trafoList.append(str(poss.args[j])+'=('+str(negs)+')*1/('+str(nenner)+')')
                            else:
                                trafoList.append(str(poss.args[j])+'=('+str(negs)+')*'+'r_'+nodeToRep+'_'+str(j)+'/('+str(nenner)+')')
                         
                    tempcounter=0
                    for trafo in trafoList:
                        #print(trafo)
                        lstrafo,rstrafo=trafo.split('=')
                        lstrafo=parse_expr(lstrafo)
                        rstrafo=parse_expr(rstrafo)
                        foundRateConst=False
                        for var in list(lstrafo.atoms()):
                            if((var not in X) & (var not in UsedStates) & (var not in UsedRC) & (not is_number(str(var))) & (not foundRateConst)):
                                RC=var
                                UsedRC.append(RC)
                                foundRateConst=True
                        if(foundRateConst):
                            rest=lstrafo/RC
                            sol=rstrafo/rest
                            eqOut.append(str(RC)+' = '+str(sol))
                            DependentRates.append(RC)
                            DependentOn.append(list(sol.atoms()))
                            tempcounter=tempcounter+1
                        else:
                            compromise2=True
                    if(compromise2):
                        #print(eq)
                        #print(nodeToRep)
                        print('Did not find appropriate transformation!')
                        for temp in range(tempcounter):
                            eqOut.pop()
                            DependentRates.pop()
                            DependentOn.pop()
                            UsedRC.pop() 
                        
                    if(not compromise2):
                        leaveCompromiseLoop=True
                        indexToRep=list(X).index(parse_expr(nodeToRep))
                        X.row_del(indexToRep)
                        SM.row_del(indexToRep)
                        for temp in range(tempcounter):
                            eq=eqOut[-(temp+1)]
                            RC=parse_expr(eq.split('=')[0])
                            sol=parse_expr(eq.split('=')[1])
                            for f in range(len(F)):
                                F[f]=F[f].subs(RC, sol)
                            print('\t'+str(RC)+' --> '+str(sol))
                
                    if((lcycle==cyclecounter) & (changeposneg==False)):
                        changeposneg=True
                        cyclecounter=0
                        print('Change pos and neg!')
                if((changeposneg) & (cyclecounter==lcycle)):
                    print('Have to use equation with minus signs!')
                    print(eq)
                    sol=solve(eq,parse_expr(nodeToRep))[0]
                    eqOut.append(nodeToRep+' = '+str(sol))
                    indexToRep=list(X).index(parse_expr(nodeToRep)) 
                    X.row_del(indexToRep)
                    SM.row_del(indexToRep)
                    DependentRates.append(parse_expr(nodeToRep))
                    DependentOn.append(list(sol.atoms()))
                    for f in range(len(F)):                
                        F[f]=F[f].subs(parse_expr(nodeToRep), sol)
            
            if((not compromise) & (not compromise2)):
                indexToRep=list(X).index(parse_expr(nodeToRep))                         
                foundRateConst=False
                #print(UsedStates)
                for var in list(FluxToTake.atoms()):
                    if((var not in X) & (var not in UsedStates) & (var not in UsedRC) & (not is_number(str(var))) & (not foundRateConst)):
                        RC=var
                        UsedRC.append(RC)
                        UsedStates.append(parse_expr(nodeToRep))
                        foundRateConst=True
                        print('   '+str(nodeToRep)+' --> '+str(RC))
                if(foundRateConst):
                    #print('   '+str(len(str(eq))))
                    sol=solve(eq, RC, simplify=False)[0]
                    eqOut.append(str(RC)+' = '+str(sol))
                    X.row_del(indexToRep)
                    SM.row_del(indexToRep)
                    #print('   '+str(len(str(sol))))
                    DependentRates.append(RC)
                    DependentOn.append(list(sol.atoms()))
                    for f in range(len(F)):                
                        F[f]=F[f].subs(RC, sol)
               
        graph={}
        for i in range(len(SM*F)):
            liste=[]
            for rate in DependentRates:
                if((SM*F)[i]!=((SM*F)[i]).subs(rate,1)):
                    for var in DependentOn[DependentRates.index(rate)]:
                        if((not is_number(str(var))) & (var!=X[i])):
                            liste.append(str(var))
            for j in range(len(X)):
                    if(j==i):
                        argus=expand((SM*F)[i]).args
                        negs=[]
                        poss=[]
                        append=False
                        for arg in argus:
                            if('-' in str(arg)):
                                negs.append(arg)
                            else:
                                poss.append(arg)
                        for el in poss:
                            if(el.subs(X[j],1)!=el):
                                append=True
                        for el in negs:
                            if(el.subs(X[j],0)!=0):
                                append=True
                        
                        if(not append):
                            In=((SM*F)[i]).subs(X[j],0)
                            Out=simplify(((SM*F)[i]-In)/X[j])
                            #print(Out)
                            if(Out!=Out.subs(X[j],1)):
                                liste.append(str(X[j]))
                        else:
                            liste.append(str(X[j]))
                    else:
                        if((SM*F)[i]!=((SM*F)[i]).subs(X[j],1)):
                            liste.append(str(X[j]))
            graph[str(X[i])]=liste
                
    
#### Solve remaining equations
    eqOut.reverse()
    #print(eqOut)
    print('Solving remaining equations ...\n')
    UsedVars=UsedRC+UsedStates
    graph={}    
    for i in range(len(SM*F)):
        liste=[]
        for j in range(len(X)):
                if(j==i):
                    argus=expand((SM*F)[i]).args
                    negs=[]
                    poss=[]
                    append=False
                    for arg in argus:
                        if('-' in str(arg)):
                            negs.append(arg)
                        else:
                            poss.append(arg)
                    for el in poss:
                        if(el.subs(X[j],1)!=el):
                            append=True
                    for el in negs:
                        if(el.subs(X[j],0)!=0):
                            append=True
                    
                    if(not append):
                        In=((SM*F)[i]).subs(X[j],0)
                        Out=simplify(((SM*F)[i]-In)/X[j])
                        #print(Out)
                        if(Out!=Out.subs(X[j],1)):
                            liste.append(str(X[j]))
                    else:
                        liste.append(str(X[j]))
                else:
                    if((SM*F)[i]!=((SM*F)[i]).subs(X[j],1)):
                        liste.append(str(X[j]))
        graph[str(X[i])]=liste
    while(graph!={}):
        #print(graph)
        node=FindNodeToSolve(graph)
        #print(node)        
        index=list(X).index(parse_expr(node))
        #print((SM*F)[index])
        sol=solve((SM*F)[index],parse_expr(node), simplify=True)
        #print(sol)
        if(sol==[]):
            foundVarToSolveFor=False
            for el in ((SM*F)[index]).atoms():
                if((not is_number(str(el))) & (el not in UsedRC) & (not foundVarToSolveFor)):
                    sol=solve((SM*F)[index],el, simplify=True)[0]
                    var=el
                    foundVarToSolveFor=True
            eqOut.insert(0,str(var)+' = '+str(sol))
            for f in range(len(F)):                
                F[f]=F[f].subs(var, sol)
            #print('Inserted')
            UsedVars.append(var)
            #print(str((SM*F)[index])+' = 0')
        else:
            eqOut.insert(0,node+' = '+str(sol[0]))
            for f in range(len(F)):                
                F[f]=F[f].subs(parse_expr(node), sol[0])
            #print(node+' = '+str(sol[0]))
        for el in graph:
            if(node in graph[el]):
                graph[el].remove(node)
        graph.pop(node)
    
    
#### Find appropriate transformations to avoid negative steady state solutions
    #eqOut=SuggestTrafos(eqOut, UsedVars)

#### Test Solution  
    print('Testing Steady State...\n')
    NonSteady=False
    #print(eqOut)
    #print(ODE)
    #print(SM*F)
    for i in range(len(ODE)):
        expr=parse_expr(str(ODE[i]))
        for j in range(len(zeroStates)):
            zeroState=zeroStates[j]
            expr=expr.subs(zeroState, 0)
        #print(len(eqOut))
        for j in range(len(eqOut)):
            ls, rs = eqOut[-(j+1)].split('=')
            #print(ls)
            ls=parse_expr(ls)
            #print(rs)
            rs=parse_expr(rs)
            expr=expr.subs(ls, rs)
            #print(simplify(expr))
        expr=simplify(expr)
        #print(expr)
        if(expr!=0):
            print('\t'+str(ODE[i]))
            print('\t'+str(expr))
            NonSteady=True
    if(NonSteady):
        print('Solution is wrong!\n')
    else:
        print('Solution is correct!\n')
    
#### Print Equations
    print('I obtained the following equations:\n')
    if(outputFormat=='M'):
        for state in zeroStates:
            print('\tinit_'+str(state)+'  "0"'+'\n')
        eqOutReturn=[]
        for i in range(len(eqOut)):
            ls, rs = eqOut[i].split('=')
            ls=parse_expr(ls)
            rs=parse_expr(rs)
            for j in range(i,len(eqOut)):
                ls2, rs2 = eqOut[j].split('=')
                rs2=parse_expr(rs2)
                rs2=rs2.subs(ls,rs)
                eqOut[j]=str(ls2)+'='+str(rs2)
            for state in Xo:
                ls=ls.subs(state, parse_expr('init_'+str(state)))
                rs=rs.subs(state, parse_expr('init_'+str(state)))
            eqOut[i]=str(ls)+'  "'+str(rs)+'"'
                            
        for i in range(len(eqOut)):
            eqOut[i]=eqOut[i].replace('**','^')
                    
        for eq in eqOut:
            print('\t'+eq+'\n')
            eqOutReturn.append(eq)            
        
    else:
        for state in zeroStates:
            print('\t'+str(state)+' = 0'+'\n')
        eqOutReturn=[]
        for eq in eqOut:
            print('\t'+eq+'\n')
            ls, rs = eq.split(' = ')
            eqOutReturn.append(ls+'='+rs)
    print('Number of Species:  '+str(nrspecies))
    print('Number of Equations:  '+str(len(eqOut)+len(zeroStates)))
    print('Number of new introduced variables:  '+str(gesnew))
    return(eqOutReturn)
    
