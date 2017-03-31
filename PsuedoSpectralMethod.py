import numpy as np
import math as m
import matplotlib.pyplot as plt
import time

#

def kroneckerDelta(i,j):
    if i == j:
        return 1
    else:
        return 0
#*******************************************************************************
        
def C(j, numberofCollocationPoints):
    return 1 + kroneckerDelta(j,0) + kroneckerDelta(j,numberofCollocationPoints)
    
#*******************************************************************************

def GenerateChebyshevCollocationPoints(numberOfPoints):
    
    a = []
    a
    c=0.0
    
    while c <= numberOfPoints:
        
        a.append(m.cos(m.pi*c/numberOfPoints))
        
        c = c + 1.0
        
    return a

        
#*******************************************************************************
    
# N is numberofColPoints
def ConstructPseudoSpectralDiffMatrix(N): 
    
    #Needs to be cleaned up and tested
    #Returns np matrix
    
    columns = N
    rows = N
    

    colPoints = GenerateChebyshevCollocationPoints(N-1)
    
    D = []
    Row = []
    Row.append((2.0*(N-1)**2 + 1.0)/6.0)
    
    #First (zero) row
    k=1
    while k < columns:
    
        temp = C(0,N-1)*(-1.0)**(0.0+k)/(C(k,N-1) * (colPoints[0] - colPoints[k]))
        Row.append(temp)
        
        k = k + 1
        
    D.append(Row)
    
    #Body 
    temp=0
    c=1
    while c < rows-1:
        
        Row = []
        k=0
        while k < columns:
            
            if not c == k:
                
                temp = C(c,N-1)*(-1.0)**(c+k)/(C(k,N-1) * (colPoints[c] - colPoints[k]))
                Row.append(temp)    
                
            else:
                
                temp = -colPoints[k]/(2.0*(1.0-colPoints[k]**2.0))
                Row.append(temp)
        
            k = k + 1

        D.append(Row)
        
        c = c + 1
        
        #Last row
        
    Row = []
    k=0
    while k < columns - 1:
        
        #print(colPoints[rows - 1] - colPoints[k])
        #Should i,j in Dij be allowed to take on zero values, such as D00? Yes! - SSExample
        temp = (C(N-1,N-1)*(-1.0)**((N-1)+k))/(C(k,N-1) * (colPoints[rows - 1] - colPoints[k]))
       
        Row.append(temp)
        
        k = k + 1
    Row.append(-D[0][0])
    D.append(Row)
    
    return np.matrix(D)
   
    
    
#*******************************************************************************

def EnforceDirichletBoundaryConditionZero(D):
    
    D = D.tolist()
    
    c=0
    while c < len(D[0]):
        
        D[0][c] = 0
        D[len(D[0]) - 1][c] = 0
        D[c][0] = 0
        D[c][len(D[0]) - 1] = 0
        
        c = c + 1
        
    return np.matrix(D)
    
#*******************************************************************************
    
def EnforceDirichletBoundaryCondition(D):
    
    
    D = D.tolist()
    
    Ds = []
    
    c = 1
    while c < len(D[0])-1:
        
        k=1
        T=[]
        
        while k < len(D[0])-1:
            
            T.append(D[c][k])
            k = k + 1
            
        Ds.append(T)
            
        c = c + 1
            
    return np.matrix(Ds)
    


#*******************************************************************************
    
def chebyshevPolynomial(x,k):
    
    #do time complexity analysis
    #test
    
    TkminusOne = x
    TkminusTwo = 1
    temp = 0
    
    if k == 0:
        return 1
        
    if k == 1:
        return x
     
    c = 2
    while c <= k:
        
        temp = TkminusOne
        TkminusOne = 2*x*TkminusOne - TkminusTwo
        TkminusTwo = temp
        
        c = c + 1
    return TkminusOne
        
#*******************************************************************************
    
def chebyshevPolynomialRecursive(x,k):
    
    #do time complexity analysis
    
    if k == 0:
        return 1
        
    if k == 1:
        return x
        
    return 2*x*chebyshevPolynomialRecursive(x,k - 1) - chebyshevPolynomialRecursive(x,k - 2)
    
        
#*******************************************************************************

#Input numpy array
def PsuedoSpectralSolutiontoHeatEquation(initCoefficients, numberofIterations):
    
    L = len(initCoefficients) 
    D =  ConstructPseudoSpectralDiffMatrix(L)
    Dxx = np.dot(D,D)
    Dxx = EnforceDirichletBoundaryCondition(Dxx)
    u = [initCoefficients] 
   
    
    #Enforce Dirichlet BCs to vector
    u = u[0].tolist()  
    u.pop(0)
    u.pop(len(u)-1)
    u = [np.array(u)]
    
    
    k=1
    t=0.001
    temp=0
    c=0
    while c < numberofIterations:
        
        temp = t*k*np.dot(Dxx,u[c]) + u[c]
        # Temp will be a matrix...[[entries]]. Must convert to array...[entries]
        temp = np.array(temp.tolist()[0])
        
        #BC
        temp[0] = 0
        temp[len(temp)-1] = 0
            
        
        u.append(temp)
    
        c = c + 1
    
    return u
        


#*******************************************************************************

#Input numpy array 
def PsuedoSpectralSolutionToWaveEquation(initPosition, initMomentum, numberofIterations):
    
    
    L = len(initPosition) 
    Dx = ConstructPseudoSpectralDiffMatrix(L)
    Dx = EnforceDirichletBoundaryCondition(Dx)
    u = []
    p=[]
        
        
    #Enforce Dirichlet BCs to initPosition
    initPosition = initPosition.tolist() 
    initPosition.pop(0)
    initPosition.pop(len(initPosition)-1)
    initPosition = np.array(initPosition)
    
    
    
    #Enforce Dirichlet BCs to initMomentum
    initMomentum = initMomentum.tolist() 
    initMomentum.pop(0)
    initMomentum.pop(len(initMomentum)-1)
    initMomentum = np.array(initMomentum)
    
    
        
    u.append(initPosition)
    p.append(initMomentum)
    
    
        
        
    
    t=0.001
    tempu = 0
    tempp = 0
    k=0
    while k < numberofIterations:
        
        tempu = t*(p[k] - np.dot(Dx,u[k])) + u[k]
        #Tempu will be a matrix...[[entries]]. Must convert to array...[entries]
        tempu = np.array(tempu.tolist()[0])
        #BC U
        tempu[0] = 0
        tempu[len(tempu)-1] = 0
        
    
    
        tempp = t*np.dot(Dx,p[k]) + p[k]
        #Tempp will be a matrix...[[entries]]. Must convert to array...[entries]
        tempp = np.array(tempp.tolist()[0])
        #BC U
        tempp[0] = 0
        tempp[len(tempp)-1] = 0        
    
    
        u.append(tempu)
        p.append(tempp)
    
    
        k = k + 1
        
    
    return u
    
    
