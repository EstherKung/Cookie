import numpy as np

def secant_method(x, y1, y2, func):
    """
                Implicit Secant Method Solver
    ----------------------------------------------------------
     - For a function 'func' in the form f(x,y) = 0, returns 
       y that solves f = 0 over the provided domain x. 
     - Provide 'x' as a numpy array.  
     - Provide two guesses, 'y1' and 'y2', as initial guesses for
       secant method.
    """
    x = np.array(x, dtype = float)
    y = np.zeros_like(x)
    for i, xi in enumerate(x):
        err1 = func(xi, y1)
        err2 = func(xi, y2)
        # print('initial error', err1, err2)
        tol = 1e-4
        while (abs(err2) > tol):
            #Secant slope
            m = (err2 - err1)/(y2 - y1)
            #Remove old values
            y1 = y2
            err1 = err2
            #Get new values
            y2 = (-err2/m + y2)
            err2 = func(xi, y2)
            # print('resulting error', err2)
        y[i] = y2
    return y

            




        
