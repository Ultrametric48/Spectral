import numpy as np
import math as m
import matplotlib.pyplot as plt

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

def EnforceBoundaryConditionZero(D):
    
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
    
def EnforceBoundaryConditionOne(D):
    
    
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
    
#Need to be numpy arrays
def waveEquationChebyshevExpansionCoefficients(initCoefficients, numberofIterations):
    
    L =len(initCoefficients)
    D =  ConstructPseudoSpectralDiffMatrix(L)
    Dxx = np.dot(D,D)
    
    temp=0
    n=2
    while n < numberofIterations:
        
        temp = (h/c)**2 * np.dot(Dxx,u[n-1]) + 2*u[n-1] - u[n-2]
        # Temp will be a matrix...[[entries]]. Must convert to array...[entries]
        temp = np.array(temp.tolist()[0])
        
        u.append(temp)
        
        n = n + 1
        
    return u


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
def PsuedoSpectralExpansionCoefficientsToHeatEquation(initCoefficients, numberofIterations):
    
    L =len(initCoefficients)
    D =  ConstructPseudoSpectralDiffMatrix(L)
    Dxx = np.dot(D,D)
    u = [initCoefficients]
    
    
    k=1
    t=0.001
    temp=0
    c=0
    while c < numberofIterations:
        
        temp = -t*k*np.dot(Dxx,u[c]) + u[c]
        # Temp will be a matrix...[[entries]]. Must convert to array...[entries]
        temp = np.array(temp.tolist()[0])
        
        u.append(temp)
    
        c = c + 1
    
    return u
        


#*******************************************************************************
# y = y' test

def PsuedoSpectralSolutionToHeatEquation(x,initCoefficients, numberofIterations):
    
    
    #Check!!!!!!!!
    
    a = []
    c = []
    u = []
    
    
    a = PsuedoSpectralExpansionCoefficientsToHeatEquation(initCoefficients, numberofIterations)
    print(a)
    print("")
    c = GenerateChebyshevCollocationPoints(len(initCoefficients)) #numberofIterations = number of points?
    print(str(len(a)) + "<---")
    
    t=0
    while t < numberofIterations:
        
       #Loop entered and run correct number of times
        
        
        n=0
        temp=0
        while n < len(initCoefficients):
            
            temp += a[t][n]*chebyshevPolynomial(x,n)
            n = n + 1
            
        u.append(temp)
        t = t + 1
    
    return u

#*******************************************************************************

        
def PsuedoSpectralSolutionToWaveEquation():
    
    #NEED TO FIX!!!
    
    #T grows very quickly as n -> inf!
    #Should x be substituted?
    
    a = waveEquationChebyshevExpansionCoefficients()
    
    colPoints = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21]
    
    U = []
    for i in a:
        u=0
        c=0
        for j in i:
            
            u += j*chebyshevPolynomial(colPoints[c],c)
            c = c + 1
            
        U.append(u)
        
    return U
        
        
        
    
    
        
        
        
    
    
    
    
        