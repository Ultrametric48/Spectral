
from Psuedo_spectral import *
import math as m

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
    

#*******************************************************************************    
            
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
    
        
               
#*******************************************************************************
        
#Test Area

def Reflecting_Wave_Test():
     
     
    C = GenerateChebyshevCollocationPoints(20)
    C.sort()
     
 
    u = []
    for i in C:
        u.append(m.exp(-20.0*i**2))
    u = np.array(u)
        
        
    phi = []
    for i in C:
        phi.append(-40.0*m.exp(-20.0*i**2))
    phi = np.array(phi)

        
    pi = []
    for i in range(len(C)):
        pi.append(5)
    pi = np.array(pi)
        
    
    
    return PsuedoSpectralSolutionToReflectingWaveEquation(u, phi, pi, 50)
 
        
  
        
        
        
        
        
        
        
        
        