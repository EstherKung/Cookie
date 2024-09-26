import numpy as np
import pandas as pd
import matplotlib.pyplot as plt 

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
                "g"        : 9.81, #acceleration due to gravity
                }

def rect_noodle_sizing(b, c, d, t_c, m_0, N, w, s_crit, t_cap, t_skin, n_skin, sig_skin, sig_tape, rho_core, g):
    
    #Load, Shear, Bending 
    y = np.linspace(0,b/2,500)
    omega = (N*m_0*g)/b
    print(f"omega: {omega}")
    V = omega*y-(1/2)*omega*b
    M = (1/2)*omega*y**2-(1/2)*omega*b*y+(1/8)*omega*b**2

    #Cap Plies
    n = (1.5*M) / (s_crit*t_cap*w*((c*t_c)-2*(n_skin*t_skin)))
    n = np.ceil(n)

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
    print(ply_pattern)

    A_cap = ply_pattern.y*ply_pattern.n #table summing area for each ply layer
    print(A_cap)
    A_tape = A_cap.to_numpy().sum() #surface area for 1 spar cap layup
    print(A_tape)


    W_cap = 6*A_tape*sig_tape
    W_skin = 6*c*b*sig_skin
    W_foam = (1/2)*(c*c*t_c*b*rho_core)

    mass = {"Spar Cap Mass":  W_cap,
            "Wing Skin Mass": W_skin,
            "Foam Core Mass": W_foam}
    
    print(mass)

    #mult by 2 to get both halves
    W_wing = 2*(1.1*(W_cap + W_skin + W_foam))
    print(W_wing)
    return W_wing


rect_noodle_sizing(**aircraftparm)
