import numpy as np
import math as m
import matplotlib.pyplot as plt
import time

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
    
    
    k=1.0
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
    
#Input numpy array 
def PsuedoSpectralSolutionToReflectingWaveEquation(Supply_init_Position, Supply_derivative_of_init_position, Supply_Pi_init_velocity, number_of_Iterations):
    
    """Begin solution with initial data from u(0,x), phi(0,x), Pi(0,x). u(0,x) should be chosen as a continuious function so that Phi is determined analytically by
    d/dx u(0,x) = Ph(0,x), the sample the two functions at the Chebyshev points for an initial data array. Pi can be chosen freely. ***numpy arrays should be input***"""
    
    L   = len(Supply_init_Position) 
    Dx  = ConstructPseudoSpectralDiffMatrix(L)
    u   = []
    pi  = []
    phi = []
    
    u.append(Supply_init_Position)
    phi.append(Supply_derivative_of_init_position)
    pi.append(Supply_Pi_init_velocity)
    
    t       = 0.001
    temp_u  = 0
    temp_pi = 0
    temp_phi= 0
    k=0
    while k < number_of_Iterations:
        
       
        temp_phi = t*np.dot(Dx,pi[k]) + phi[k]
        #Temp will be a matrix...[[entries]]. Must convert to array...[entries]
        temp_phi = np.array(temp_phi.tolist()[0])
        
        
        
        temp_pi = t*np.dot(Dx,phi[k]) + pi[k]
        #Temp will be a matrix...[[entries]]. Must convert to array...[entries]
        temp_pi = np.array(temp_pi.tolist()[0])
        
        
        
        temp_u = t*pi[k] + u[k]
        #Temp will be a matrix...[[entries]]. Must convert to array...[entries]
        temp_u = np.array(temp_u.tolist()[0])
        
        
        
        #tempu = t*(p[k] - np.dot(Dx,u[k])) + u[k]
        # Reflecting BC uxx0 = 0
        #tempu = t* np.dot(Dx,u[k])
        #Tempu will be a matrix...[[entries]]. Must convert to array...[entries]
        #tempu = np.array(tempu.tolist()[0])
        # Reflecting boundary conditions ux0 = 0 => uxx0 = 0
        # Strong boundary imposition: force u0 to be the boundary value you want
        # In this case, we're enforcing boundary conditions on the left boundary, the left most element u0.
        
        
        
        
        
        
        
    
        #tempp = t*np.dot(Dx,p[k]) + p[k]
        #Tempp will be a matrix...[[entries]]. Must convert to array...[entries]
        #tempp = np.array(tempp.tolist()[0])
        # Reflecting boundary conditions p0 -> -p0
        # Strong boundary imposition: force p0 to be the boundary value you want
        
    
    
        #u.append(tempu)
        #p.append(tempp)
    
        phi.append(temp_phi)
        pi.append(temp_pi)
        u.append(temp_u)
    
        k = k + 1
        
    
    return u
    
def FullsolnHeatEquation():
    
    
    C = GenerateChebyshevCollocationPoints(20)
    C.sort()

    
    A = []

    for i in C:
    
        A.append(m.exp(-20.0*i**2))
    
    A = np.array(A)
    SSP = PsuedoSpectralSolutiontoHeatEquation(A, 100)

    
    C.pop(0)
    C.pop()
    
    for k in SSP:
        plt.plot(C,k)
    
    
def Wavesoln():


    C = GenerateChebyshevCollocationPoints(20)
    C.sort()
    
    A = []
    for i in C:
        A.append(m.exp(-20.0*i**2))
    A = np.array(A)
    
    V = []
    for i in C:
        V.append(5.0)
        
    V = np.array(V)
    
    SSP = PsuedoSpectralSolutionToWaveEquation(A, V, 500) 
       
    
    C.pop(0)
    
    C.pop()
    
    for i in range(0,30):
        plt.plot(C,SSP[i])
    
        
    
        
        
        
    
    
        
        
        
    
    
    
    
        