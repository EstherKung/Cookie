import numpy as np
import pandas as pd
from bisect import bisect_left
import matplotlib.pyplot as plt

#some misc/specific functions to keep the main script from being too overwhelming
def xCG_calc(item: list): #handy help calculate your xCG of a sectional component (discreised) use a list of lists [[mass, location]]
    loc = 0;  mass = 0; 
    for x in item:
        loc = loc+ (x[0] * x[1])
        mass = mass + x[0]
    xCG = loc/mass
    return xCG

def nearest(xloc, des: str, xcoords, ycoords): #returns the x, y value pair for airfoil interpolation
    xcoords = xcoords.tolist()
    ycoords = ycoords.tolist() #bruhh so this was the problem ???
    pos = bisect_left(xcoords, xloc)
    
    #print("this is the bisect location: ", pos)
    if pos == 0:
        #print("maybe")
        return xcoords[0], ycoords[0]
    if pos == len(xcoords):
        #print("um")
        return xcoords[len(xcoords)-1], ycoords[len(ycoords)-1] #you thought index of -1 would work...
    before = xcoords[pos - 1]
    if xloc == xcoords[pos]:
        after = xcoords[pos + 1]
    else: after = xcoords[pos]
    
    match des:
        case 'prev': 
            return float(before), float(ycoords[pos - 1])
        case 'next':
            return float(after), float(ycoords[pos])

def csec(afile):
    data = pd.read_table(afile,
                         sep='\s+',skiprows=[0],names=['x','y'],index_col=False)
    area = np.trapezoid(data['y'], data['x'])
    return(abs(area))

def centroid(afile):
    data = pd.read_table(afile,
                        sep='\s+', skiprows=[0],names =['x','y'], index_col=False)

    loc = np.trapezoid(data['y']*data['x'], data['x'])
    weight = np.trapezoid(data['y'], data['x'])
    centroid = loc/weight
    return(centroid)

def thickness(afile):
    data = pd.read_table(afile,
                        sep='\s+',skiprows=[0],names=['x','y'],index_col=False)
    split = int(data['x'].idxmin())
    top = data.iloc[:split + 1, :]
    top = top.sort_values(by=['x'], ascending=True)
    bottom = data.iloc[split:, :]
    bottom.sort_values(by=['x'], ascending=True)
    if len(top) > len(bottom): master = top; check = bottom
    else: master = bottom; check = top #default, bottom is the master
    
    t_max = 0
    for x, y in zip(master['x'], master['y']):
        if x == master['x'].tolist()[0]: continue
        if x == master['x'].tolist()[-1]: continue
        
        xprev, yprev = nearest(x, 'prev', check['x'], check['y'])
        xnext, ynext = nearest(x, 'next', check['x'], check['y'])
        #plt.plot([xprev, xnext], [yprev, ynext])

        yinterp = (ynext-yprev)/(xnext - xprev) * (x - xprev) + yprev

        t = abs(yinterp) + abs(y)
        if t > t_max:
            t_max = t
    return round(t_max, 4)

def camber(afile):
    print('no')

def planform(inputs: list): #to find planform parameters
    """in the order 
    [AR, b, S, taper_ratio, c_r, c_t]
    and for the ones that aren't defined, set to -1
    """
    AR = inputs[0]; b = inputs[1]; S = inputs[2]
    taper = inputs[3]; c_r = inputs[4]; c_t = inputs[5]
    #inputs = [AR, b, S, taper, c_r, c_t]

    num_inputs = inputs.count(-1)
    count_AR_param = inputs[0:3].count(-1)
    count_taper_param = inputs[3:6].count(-1)

    if num_inputs > 3:
        return "not enough information"
    if inputs.count(-1) == 3:
        if count_AR_param == 0:
            return "invalid combination"
        if count_taper_param == 0:
            return "invalid combination"
    
    #Funny checks
    no_chord = (c_r == -1 and c_t == -1)
    no_bs = (b == -1 and S == -1)
    no_ratio = (AR == -1 and taper == -1)

    missing_chord = (c_r == -1 or c_t == -1)
    missing_bs = (b == -1 or S == -1)
    missing_ratio = (AR == -1 or taper == -1)

    def AR_param2(AR, b, S): #use when given 2 out of 3 AR parameters
        if AR == -1:
            AR = b**2/S
        elif b == -1:
            b = np.sqrt(S*AR)
        elif S == -1: 
            S = b**2/AR
        
        return [AR, b, S]
    
    def taper_param2(taper, c_r, c_t): #use when given 2 out of 3 taper parameters
        if taper == -1:
            taper == c_t/c_r
        elif c_t == -1:
            c_t = taper * c_r
        elif c_r == -1:
            c_r = c_t/taper
        
        return [taper, c_r, c_t]

    def AR_param1(AR, b, S): 
        #use when given 1 out of 3 AR parameters
        #necessitates we have 2 out of 3 taper parameters
        #so let's find those first
        if b == -1:
            b = S/((c_r + c_t)/2)
        if S == -1:
            S = b*(c_r + c_t)/2
        AR = b/S
        return [AR, b, S]

    def taper_param1(taper, c_r, c_t): 
        #use when given 1 out of 3 taper parameters
        #necessitates we have 2 out of 3 AR parameters
        if c_t == -1:
            c_t = 2*S/b - c_r
        if c_r == -1:
            c_r = 2*S/b - c_t
        taper = c_t/c_r
        return[taper, c_r, c_t]
    
    if no_chord:
        [AR, b, S] = AR_param2(AR, b, S)
        c_r = 2*S/b/(taper + 1)
        c_t = taper*(2*S/b/(taper + 1))

    elif no_bs:
        [taper, c_r, c_t] =  taper_param2(taper, c_r, c_t)
        b = AR*(c_r + c_t)/2
        S = b*(c_r + c_t)/2
    
    elif no_ratio:
        if num_inputs == 3:
            if missing_chord:
                AR = b/S
                if c_t == -1:
                    c_t = 2*S/b - c_r
                if c_r == -1:
                    c_r = 2*S/b - c_t
                taper = c_t/c_r
            if missing_bs:
                taper = c_t/c_r
                if b == -1:
                    b = S/((c_r+c_t)/2)
                if S == -1:
                    S = b*(c_r + c_t)/2
                AR = b/S
        else:
            AR = b/S
            taper = c_t/c_r

    elif count_AR_param == 2:
        [taper, c_r, c_t] = taper_param2(taper, c_r, c_t)
        [AR, b, S] = AR_param1(AR, b, S)
    
    elif count_taper_param == 2:
        [AR, b, S] = AR_param2(AR, b, S)
        [taper, c_r, c_t] = taper_param1(taper, c_r, c_t)

    return float(AR), float(b), float(S), float(taper), float(c_r), float(c_t)