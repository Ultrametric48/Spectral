

def delta(x,y):
    
    if x==y:
        return 0
    else:
        return 1

def edit_distance(s,w):
    
    if len(s) == 0 or len(w) == 0:
        return
        
    
    x = s[len(s)-1]
    y = w[len(w)-1]
    a = s[:len(s) - 1]
    b = w[:len(s) - 1]
    
    return min(edit_distance(a,b) + delta(x,y), edit_distance(s,b) + 1, edit_distance(a,w) + 1)
    
    
def edit_distance_matrix(s,w):
    
    A = []
      
    for i in range(len(s)):
        B.append(i)
        
        
        