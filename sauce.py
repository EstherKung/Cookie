#cooking
import numpy as np
import pandas as pd
from wing_def import *
from plane_def import *
from AVLWrapper.Code.Codebase import *
from LG_def import *
 
w = 0; 
def avl_sm(plane: Plane):
    #for avl stability analysis
    plane.to_avl()
    ad = ADaX(planename='Trainer.avl')
    ad.output_config(['s Xnp', 't Cref'])
    #ad.use_run("AVLWrapper/AVL/Planes/Funny_Run_File.run", 1)
    ad._load(); ad._oper()
    ad.vcv('d1 pm 0'); ad.vcv('d2 ym 0'); ad.vcv('a a 0')
    ad._x(); ad._save()
    results = ad.run(print_output=False)
    Xnp = float(results['Xnp']); mac = float(results['Cref']); 
    plt.plot(Xnp, 0, 'x')

    CG_chord = (CG - wing.x)/ wing.planform['cr']

    return sm(Xnp, CG, mac), CG_chord, mac

def CG_routine(wing_loc, CG): #tallies up the CG of the new configuration
    global wing, fuse_nose, fuse_main, fuse_aft, NLG, MLG, htail, vtail
    wing = Wing(wing_loc, "foam", [AR, -1, S, tr, -1, -1], [-1, 0, 0, 0], afile)
    fuse_nose = Fuselage(0, "nose", 8, [3, 3], [6, 6])

    #making the main fuselage section end where the wing ends
    main_end = wing_loc + wing.planform['cr'] - fuse_nose.length
    fuse_main = Fuselage(fuse_nose.length, "main", main_end, [6, 6])
    fuse_aft = Fuselage(fuse_main.x + fuse_main.length, "aft", aft_sec, [6, 6], [2, 2])

    total_length = fuse_nose.length + fuse_main.length + fuse_aft.length
    LG_height, MLG_aft, NLG_fwd_wheel, NLG_fwd_strut, MLG_side = LG_def(CG, MTOW, total_length)

    NLG = LandingGear(CG - NLG_fwd_wheel, "NLG", point_mass = True)
    MLG = LandingGear(CG + MLG_aft, "MLG", point_mass = True) #MLG 10% of the total length behind the CG location

    emp_xloc = fuse_aft.x + fuse_aft.length - 0.75* htail.c #assuming empennage plate starts where the fuselage ends
    l_emp = emp_xloc - (wing.x + wing.planform['cr']/4); #distance between c/4 of wing and c/4 of the tail
    S_h = V_h * S * wing.planform['cr'] / l_emp #that c is the mac, supposedly
    S_v = V_v * S * wing.planform['b'] / l_emp

    htail = Component(emp_xloc, "balsa", type = 'area', AR = ARh, area = S_h, aero = 'hstab')
    vtail = Component(emp_xloc, "balsa", type = 'area', AR = ARv, area = S_v, aero = 'vstab')

    components = [fuse_nose, fuse_main, fuse_aft,
                NLG, MLG, 
                wing, htail, vtail,
                batt, ESC, prop, motor]

    Trainer = Plane("Trainer", components, Sref = S)

    CG_chord = (CG - wing.x)/ wing.planform['cr'] 
    #print(f"CG @ {CG_chord:1.1%} chord. ")
    return Trainer, CG_chord

afile = "airfoil_library/S4233.dat"

"Point on Constraint Diagram-------------------------------"
TW = 0.2046                                          #lb/lb--
WS = 1.245            / 144 #lb/ft^2 * (ft/12in)^2 = lb / in^2
print(
    f"""Point Selected from the Constraint Diagram: 
    TW = {TW} lb/lb, WS = {WS: 1.4f} lb/sq.in"""
    )
"----------------------------------------------------------"

"Specifications--------------------------------------------"
T = 5                                                   #lbs
MTOW = T/TW                                             #lbs
print(f"initial MTOW = {MTOW: 1.2f} lbs.")
"----------------------------------------------------------"


"Tail Volume Coefficients & AR-----------------------------"
AR = 6; tr = 0.4; 
ARh = 0.5*AR; V_h = 0.40; 
ARv = 1.5; V_v = 0.04; 
"----------------------------------------------------------"

"Default Components on Plane-------------------------------"
batt = Component(4, "battery", point_mass = True)
Rx = Component(4, "Rx", point_mass = True)
Rx_batt = Component(8, "Rx batt", point_mass = True)
ESC = Component(4, "ESC", point_mass = True)
prop = Component(0, "prop", point_mass = True)
motor = Component(1, "motor", point_mass = True)
"----------------------------------------------------------"


#Routine Setup--------------------------------------------||
S =  MTOW / WS
wing_loc = 24; 
wing = Wing(wing_loc, "foam", 
            [AR, -1, S, tr, -1, -1], [-1, 0, 0, 0], afile)
#initially assume CG @ half root chord behind the wing
CG = wing.x + wing.planform['cr']/2
aft_sec = 24 #start with aft fuselage length = 24 in.
S_h = V_h * wing.planform['S'] * wing.planform['cr'] / aft_sec
htail = Component(aft_sec, "balsa", type = 'area', AR = ARh, area = S_h, aero = 'hstab') #dummy htail for loop to work

#--------------------------------------------------------||


#Fixed Point Iteration Mass Convergence Study
res = 1
print(
    "Performing Routine 1: Mass Convergence Study..."
    )

while res > 1e-3: 
    MTOW_old = MTOW
    S =  MTOW / WS

    wing_loc = 24; #this has to be iterated for CG placment in second iteration
    wing = Wing(wing_loc, "foam", [AR, -1, S, tr, -1, -1], [-1, 0, 0, 0], afile)

    fuse_nose = Fuselage(0, "nose", 8, [3, 3], [6, 6])

    #making the main fuselage section end where the wing ends
    main_end = wing_loc + wing.planform['cr'] - fuse_nose.length
    fuse_main = Fuselage(fuse_nose.length, "main", main_end, [6, 6])
    fuse_aft = Fuselage(fuse_main.x + fuse_main.length, "aft", aft_sec, [6, 6], [2, 2])

    total_length = fuse_nose.length + fuse_main.length + fuse_aft.length

    LG_height, MLG_aft, NLG_fwd_wheel, NLG_fwd_strut, MLG_side = LG_def(CG, MTOW, total_length)
    
    NLG = LandingGear(CG - NLG_fwd_wheel, "NLG", point_mass = True)
    MLG = LandingGear(CG + MLG_aft, "MLG", point_mass = True) #MLG 10% of the total length behind the CG location

    emp_xloc = fuse_aft.x + fuse_aft.length - 0.75* htail.c #assuming empennage plate starts where the fuselage ends
    l_emp = emp_xloc - (wing.x + wing.planform['cr']/4); #distance between c/4 of wing and c/4 of the tail
    S_h = V_h * S * wing.planform['cr'] / l_emp #that c is the mac
    S_v = V_v * S * wing.planform['b'] / l_emp

    htail = Component(emp_xloc, "balsa", type = 'area', AR = ARh, area = S_h, aero = 'hstab')
    vtail = Component(emp_xloc, "balsa", type = 'area', AR = ARv, area = S_v, aero = 'vstab')

    components = [fuse_nose, fuse_main, fuse_aft,
                NLG, MLG, 
                wing, htail, vtail,
                batt, ESC, prop, motor, Rx, Rx_batt]

    Trainer = Plane("Trainer", components, Sref = S)
    Trainer.plane_plot(components) #plane plot

    CG, mtot = Trainer.CG_tally()
    MTOW = Trainer.MTOW_tally()
    #print(f"MTOW = {MTOW: 1.2f} lbs.")

    res = abs(MTOW_old - MTOW)
    #print(res)
print(f"Mass Convergence Study Done. MTOW = {MTOW: 1.2f} lbs.")

#CG Placement Study by varying wing placement
print("Performing Routine 2: CG Placement Study...")
CG_chord = (CG - wing.x)/ wing.planform['cr']

#---Configure upper & Lower Limits---
CG_upper = 0.35                    #
CG_lower = 0.27                    #
#------------------------------------


while  CG_chord < CG_lower or CG_chord > CG_upper: 
    CG_chord_old = CG_chord
    if CG_chord_old > CG_upper:
        wing_loc += 1
    if CG_chord_old < CG_lower:
        wing_loc -= 1
    
    Trainer, CG_chord = CG_routine(wing_loc, CG)
    Trainer.plane_plot(Trainer.components)
    CG, mtot = Trainer.CG_tally()
    MTOW = Trainer.MTOW_tally()

    plt.savefig(f"iteration/{w}.png"); w += 1; 
print(f"CG Placement Study done. CG @ {CG_chord:1.1%} chord; MTOW = {MTOW: 1.2f} lbs.")

print("Send to AVL for stability...")
static_margin, CG_chord, mac = avl_sm(Trainer)

#Configure upper & lower limits-------
sm_upper = 0.15
sm_lower = 0.07
#-------------------------------------

while sm_lower < static_margin > sm_upper:
    if static_margin < sm_lower: #too close
        aft_sec += 0.1 #increase empennage length
        #V_h += 0.001 #increase hstab size - change volume coeff, hold off for now
        pass
    if static_margin > sm_upper: #too large
        aft_sec -= 0.1 #decrease empennage length
        #V_h -= 0.001 #decrease hstab size - change volume coeff, hold off for now
        pass
    fuse_aft = Fuselage(fuse_main.x + fuse_main.length, "aft", aft_sec, [6, 6], [2, 2])

    total_length = fuse_nose.length + fuse_main.length + fuse_aft.length
    LG_height, MLG_aft, NLG_fwd_wheel, NLG_fwd_strut, MLG_side = LG_def(CG, MTOW, total_length)
    
    NLG = LandingGear(CG - NLG_fwd_wheel, "NLG", point_mass = True)
    MLG = LandingGear(CG + MLG_aft, "MLG", point_mass = True) #MLG 10% of the total length behind the CG location

    emp_xloc = fuse_aft.x + fuse_aft.length - 0.75* htail.c #assuming empennage plate starts where the fuselage ends
    l_emp = emp_xloc - (wing.x + wing.planform['cr']/4); #distance between c/4 of wing and c/4 of the tail
    S_h = V_h * S * wing.planform['cr']/ l_emp #that c is supposed to be the mac but it looks tiny with the mac ??
    S_v = V_v * S * wing.planform['b'] / l_emp

    htail = Component(emp_xloc, "balsa", type = 'area', AR = ARh, area = S_h, aero = 'hstab')
    vtail = Component(emp_xloc, "balsa", type = 'area', AR = ARv, area = S_v, aero = 'vstab')
    components = [fuse_nose, fuse_main, fuse_aft,
                NLG, MLG, 
                wing, htail, vtail,
                batt, ESC, prop, motor, Rx, Rx_batt]
    Trainer = Plane("Trainer", components, Sref = S)
    Trainer.plane_plot(components) #plane plot
    CG, mtot = Trainer.CG_tally()
    MTOW = Trainer.MTOW_tally()
    static_margin, CG_chord, mac = avl_sm(Trainer)
    plt.savefig(f"iteration/{w}.png"); w += 1; 

    print(f"static margin: {static_margin:.1%}, CG @ {CG_chord:1.1%} chord.")


while  CG_chord < CG_lower or CG_chord > CG_upper: 
    #2nd iteration if necessary
    CG_chord_old = CG_chord
    if CG_chord_old > CG_upper:
        wing_loc += 0.1
    if CG_chord_old < CG_lower:
        wing_loc -= 0.1

    Trainer, CG_chord = CG_routine(wing_loc, CG)
    Trainer.plane_plot(Trainer.components) #plane plot
    CG, mtot = Trainer.CG_tally()
    MTOW = Trainer.MTOW_tally()
    static_margin, CG_chord, mac = avl_sm(Trainer)
    plt.savefig(f"iteration/{w}.png"); w += 1; 

    print(f"static margin: {static_margin:.1%}, CG @ {CG_chord:1.1%} chord.")

#Export Results
export = open("plane_geometry.txt", 'w')
export.write(f"""[Results]
static margin: {static_margin:.1%}, CG @ {CG_chord:1.1%} chord.
MTOW estimate = {MTOW: 1.2f} lbs.\nPlanform Area S = {wing.planform['S']: 1.4f}\n
Fuselage Definitions: Total Length {total_length: 1.2f} in.\n
\tNose Section:\t{fuse_nose.length: 1.4f} in.\n\tMain Section:\t{fuse_main.length: 1.4f} in.\n\tAft Section:\t{fuse_aft.length: 1.4f} in.\n
Wing Definitions: Location from nose: {wing.x:1.2f} (in.)\n
\tWingspan:\t{wing.planform['b']: 1.2f}\n\tRoot chord:\t{wing.planform['cr']: 1.2f}\n\tTip chord:\t{wing.planform['ct']: 1.2f}\n\tLE sweep:\t{float(np.rad2deg(wing.angles['sweep'])): 1.2f} deg.\n
Empennage Definitions: (in.)\n
\thstab:\n\tchord\t{htail.c: 1.4f}\n\tspan\t{htail.b: 1.4f}\n\tvstab:\n\tchord\t{vtail.c:1.4f}\n\tspan\t{vtail.b: 1.4f}\n
Landing Gear Placement: (dist. from nose, in.)\n
\tNLG:\t{NLG.x: 1.4f}\n\tMLG:\t{MLG.x: 1.4f}
""")
export.close()
out = open("plane_geometry.txt")
print(out.read())
print(CG)
#plt.show()