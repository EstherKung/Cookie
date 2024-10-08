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

def csec(afile): #unit area
    data = pd.read_table(afile,
                         sep='\s+',skiprows=[0],names=['x','y'],index_col=False)
    area = np.trapz(data['y'], data['x'])
    return(abs(area))

def centroid(afile):
    data = pd.read_table(afile,
                        sep='\s+', skiprows=[0],names =['x','y'], index_col=False)

    loc = np.trapz(data['y']*data['x'], data['x'])
    weight = np.trapz(data['y'], data['x'])
    centroid = loc/weight
    return(centroid)

def thickness(afile): #finds (thickness, % chord location)
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
    xloc = -1
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
            xloc = x
    
    if xloc == -1: return "Error in finding maximum thicness location"
    return round(t_max, 4), xloc

def camber(afile): #incomplete
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

aircraftparm = {"b"        : 1.83, #wingspan
                "c"        : 0.305, #chord
                "d"        : 0.0366, #(optional) airfoil thickness 
                "t_c"      : 0.09, #percent airfoil thickness
                "m_0"      : 12, #max takeoff weight
                "N"        : 6, #ultimate load factor
                "w"        : 0.05, #cap width
                "s_crit"   : 7e+8, #carbon fiber critical compressive stress
                "t_cap"    : 0.0001016, #cap ply thickness
                "t_skin"   : 0.0001016, #skin ply thickness
                "n_skin"   : 2, #of skin plies
                "sig_skin" : 0.093, #skin fabric weight
                "sig_tape" : 0.08, #cap fabric weight
                "rho_core" : 16.0185, #foam core density
                "G_core"   : 1.2e7, #foam core shear modulus
                "g"        : 9.81, #acceleration due to gravity
                }

#let it take in arguments of wingspan, chord, % airfoil thickness
def foamcore(b, c, t_c, m_0, N = aircraftparm['N'], w = aircraftparm['w'], 
             s_crit = aircraftparm['s_crit'], t_cap = aircraftparm['t_cap'], t_skin = aircraftparm['t_skin'], 
             n_skin = aircraftparm['n_skin'], sig_skin = aircraftparm['sig_skin'], sig_tape = aircraftparm['sig_tape'], 
             rho_core = aircraftparm['rho_core'], G_core = aircraftparm['G_core'], g = 9.81, afile = "airfoil_library/AG25.dat"):
    """code is based in metric!!!"""
    b = b * 0.0254 #inches to m
    c = c * 0.0254 #inches to m
    m_0 = m_0 * 0.453592 #lbs to kg
    
    #Load, Shear, Bending 
    y = np.linspace(0,b/2,500)
    omega = (N*m_0*g)/b
    # print(f"omega: {omega}")
    V = omega*y-(1/2)*omega*b
    M = (1/2)*omega*y**2-(1/2)*omega*b*y+(1/8)*omega*b**2

    #Cap Plies
    n = (1.5*M) / (s_crit*t_cap*w*((c*t_c)-2*(n_skin*t_skin)))
    n = np.ceil(n)

    #Shear Web Check
    d_s = (N*m_0*g*b) / (8*G_core*w*(c*t_c-2*n_skin*t_skin-n[1]*t_cap))
    # print(f"Shear Deflection: {d_s}")

    #Tabulation & Display:
    data = {
        "y": y,
        "V": V,
        "M": M,
        "n": n,
    }
    
    wing_df = pd.DataFrame(data)
    # print(wing_df)

    ply_pattern = wing_df.drop_duplicates(subset=['n'], keep='first')
    ply_pattern = ply_pattern.drop(['V', 'M'], axis=1)
    # print(ply_pattern)

    # plt.scatter(y, V)
    # plt.scatter(y, M)
    # plt.show()
    
    # plt.scatter(y, n)
    # plt.show()

    # wing_df.to_csv('wing_s.csv')

    #Wing Mass Calc. 
    ply_pattern.loc[:,'y'] = ply_pattern.y.shift(-1)
    ply_pattern = ply_pattern[:-1]
    # print(f"layup sched: {ply_pattern}")

    A_cap = ply_pattern['y'].sum() #tape length for 1 cap layup
    # print(f"1 Cap length:{A_cap}")
    A_tape = A_cap*w #surface area for 1 spar cap layup
    # print(f"1 Cap Area: {A_tape}")


    W_cap = 6*A_tape*sig_tape
    W_skin = 6*c*b*sig_skin
    #W_foam = (1/2)*(c*c*t_c*b*rho_core)*0.76 #ONLY USE IF DOING THICKNESS SWEEP
    W_foam = rho_core * csec(afile) * c** 2 * b * 0.5

    mass = {"Spar Cap Mass":  W_cap,
            "Wing Skin Mass": W_skin,
            "Foam Core Mass": W_foam}

    #mult by 2 to get both halves
    W_wing = 2*(1.1*(W_cap + W_skin + W_foam))

    W_wing = W_wing / 14.5939

    return W_wing