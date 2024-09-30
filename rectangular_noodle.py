import numpy as np
import pandas as pd
import matplotlib.pyplot as plt 

#Plane properties
aircraftparm = {"b"        : 1.83, #wingspan
                "c"        : 0.2814, #chord
                "t_c"      : 0.076, #percent airfoil thickness
                "m_0"      : 5, #max takeoff weight
                "N"        : 6, #ultimate load factor
                "w"        : 0.05, #cap width
                "s_crit"   : 7e+8, #carbon fiber critical compressive stress
                "t_cap"    : 0.0001016, #cap ply thickness
                "t_skin"   : 0.0001016, #skin ply thickness
                "n_skin"   : 2, #of skin plies
                "sig_skin" : 0.08, #skin fabric weight
                "sig_tape" : 0.08, #cap fabric weight
                "rho_core" : 35, #foam core density
                "G_core"   : 1.2e7, #foam core shear modulus
                "g"        : 9.81, #acceleration due to gravity
                "afile"    : "AG25.dat"
                }

#Airfoil Cross Section
def csec(afile):
    data = pd.read_table(afile,
                         sep='\s+',skiprows=[0],names=['x','y'],index_col=False)
    area = np.trapz(data['y'], data['x'])
    area_norm = abs(area)*aircraftparm["c"]**2
    return(area_norm)

#Wing Sizing 
def rect_noodle_sizing(b, c, t_c, m_0, N, w, s_crit, t_cap, t_skin, n_skin, sig_skin, sig_tape, rho_core, G_core, g, afile):
    
    #Load, Shear, Bending 
    y = np.linspace(0,b/2,500)
    omega = (N*m_0*g)/b
    print(f"omega: {omega}")
    V = omega*y-(1/2)*omega*b
    M = (1/2)*omega*y**2-(1/2)*omega*b*y+(1/8)*omega*b**2

    #Cap Plies
    n = (1.5*M) / (s_crit*t_cap*w*((c*t_c)-2*(n_skin*t_skin)))
    n = np.ceil(n)

    #Shear Web Check
    d_s = (N*m_0*g*b) / (8*G_core*w*(c*t_c-2*n_skin*t_skin-n[1]*t_cap))
    print(f"Shear Deflection: {d_s}")

    #Tabulation & Display:
    data = {
        "y": y,
        "V": V,
        "M": M,
        "n": n,
    }
    
    wing_df = pd.DataFrame(data)
    print(wing_df)

    ply_pattern = wing_df.drop_duplicates(subset=['n'], keep='first')
    ply_pattern = ply_pattern.drop(['V', 'M'], axis=1)
    print(ply_pattern)

    plt.scatter(y, V)
    plt.scatter(y, M)
    plt.show()
    
    plt.scatter(y, n)
    plt.show()

    wing_df.to_csv('wing_s.csv')

    #Wing Mass Calc. 
    ply_pattern.loc[:,'y'] = ply_pattern.y.shift(-1)
    ply_pattern = ply_pattern[:-1]
    print(f"layup sched: {ply_pattern}")

    A_cap = ply_pattern['y'].sum() #tape length for 1 cap layup
    print(f"1 Cap length:{A_cap}")
    A_tape = A_cap*w #surface area for 1 spar cap layup
    print(f"1 Cap Area: {A_tape}")


    W_cap = 6*A_tape*sig_tape
    W_skin = 6*c*b*sig_skin
    #W_foam = (1/2)*(c*c*t_c*b*rho_core)*0.76 #ONLY USE IF DOING THICKNESS SWEEP
    W_foam = rho_core * csec(afile) * aircraftparm["b"] * 0.5


    mass = {"Spar Cap Mass":  W_cap,
            "Wing Skin Mass": W_skin,
            "Foam Core Mass": W_foam}
    
    print(mass)

    #mult by 2 to get both halves
    W_wing = 2*(1.1*(W_cap + W_skin + W_foam))
    print(W_wing)

    return W_wing


csec("AG25.dat")
rect_noodle_sizing(**aircraftparm)
