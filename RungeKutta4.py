

class ODESolver(object):
    
    
    
    # of the form y' = f(y(t),t)
    
    f = 0
    yapprox = 0
    y1 = 0
    y2 = 0 
    y3 = 0
    y4 = 0
    
    # constructor takes in function f and I.C.
    def __init__(self, f, initialCondition):
        self.yapprox = initialCondition
        self.f = f
        
        
    def k1(self,t):
        
        return self.f(self.yapprox,t)
        
    def k2(self):
        
        y1Intermediate = yapprox + k1()*h/2.0
        
        return