#cooking
import numpy as np
import pandas as pd
from wing_def import *
from plane_def import *
from AVLWrapper.Code.Codebase import *

afile = "C:\\Users\\eneiche\\Documents\\airfoil-sauce\\airfoil_library\\Selig\\S7055.dat"

"Point on Constraint Diagram-------------------------------"
TW = 0.225                                          #lb/lb--
WS = 1.1            / 144 #lb/ft^2 * (ft/12in)^2 = lb / in^2
print(
    f"Point Selected from the Constraint Diagram: TW = {TW} lb/lb, WS = {WS: 1.4f} lb/sq.in"
    )
"----------------------------------------------------------"

"Specifications--------------------------------------------"
T = 5                                                   #lbs
MTOW = T/TW                                             #lbs
print(f"initial MTOW = {MTOW: 1.2f} lbs.")
"----------------------------------------------------------"


"Tail Volume Coefficients & AR"
AR = 6; tr = 0.4; 
ARh = 4; V_h = 0.40; 
ARv = 1.5; V_v = 0.04; 

#Routine Setup
res = 1
S =  MTOW / WS
wing_loc = 12; 
wing = Wing(wing_loc, "foam", 
            [AR, -1, S, tr, -1, -1], [-1, 0, 0, 0], afile)
CG = wing.x + wing.planform['cr']/2 #initially assume it's a half root chord behind the wing

print(
    "Performing Routine 1: Mass Convergence Study..."
    )

#iterate one point on the constraint diagram
#to-do: cycle for optimal point on constraint diagram
while res > 1e-3: 
    MTOW_old = MTOW
    S =  MTOW / WS

    wing_loc = 12; #this has to be iterated for CG placment in second iteration
    wing = Wing(wing_loc, "foam", [AR, -1, S, tr, -1, -1], [-1, 0, 0, 0], afile)

    fuse_nose = Fuselage(0, "nose", 8, [3, 3], [6, 6])

    #making the main fuselage section end where the wing ends
    main_end = wing_loc + wing.planform['cr'] - fuse_nose.length
    fuse_main = Fuselage(fuse_nose.length, "main", main_end, [6, 6])
    fuse_aft = Fuselage(fuse_main.x + fuse_main.length, "aft", 24, [6, 6], [2, 2])

    total_length = fuse_nose.length + fuse_main.length + fuse_aft.length
    NLG = LandingGear(fuse_main.x, "NLG", point_mass = True)
    MLG = LandingGear(CG + 0.1*total_length, "MLG", point_mass = True) #MLG 10% of the total length behind the CG location

    emp_xloc = fuse_aft.x + fuse_aft.length #assuming empennage plate starts where the fuselage ends
    l_emp = emp_xloc - (wing.x + wing.planform['cr']/4); #distance between c/4 of wing and c/4 of the tail
    S_h = V_h * S * wing.planform['cr'] / l_emp #that c is the mac
    S_v = V_v * S * wing.planform['b'] / l_emp

    htail = Component(emp_xloc, "balsa", type = 'area', AR = ARh, area = S_h, aero = 'hstab')
    vtail = Component(emp_xloc, "balsa", type = 'area', AR = ARv, area = S_v, aero = 'vstab')

    batt = Component(fuse_main.x, "battery", point_mass = True)
    ESC = Component(fuse_main.x, "ESC", point_mass = True)
    prop = Component(0, "prop", point_mass = True)
    motor = Component(2, "motor", point_mass = True)

    components = [fuse_nose, fuse_main, fuse_aft,
                NLG, MLG, 
                wing, htail, vtail,
                batt, ESC, prop, motor]

    Trainer = Plane("Trainer", components, Sref = S)
    Trainer.plane_plot(components) #plane plot

    CG, mtot = Trainer.CG_tally()
    MTOW = Trainer.MTOW_tally()
    #print(f"MTOW = {MTOW: 1.2f} lbs.")

    res = abs(MTOW_old - MTOW)
    #print(res)
print(f"Mass Convergence Study Done. MTOW = {MTOW: 1.2f} lbs.")

print("Performing Routine 2: CG Placement Study...") #by varying the wing placement
CG_chord = (CG - wing.x)/ wing.planform['cr'] #target this to be 30% - 35%

#---Configure upper & Lower Limits---
CG_upper = 0.33                    #
CG_lower = 0.25                    #
#------------------------------------

while CG_lower > CG_chord < CG_upper: 
    CG_chord_old = CG_chord
    if CG_chord_old > CG_upper:
        wing_loc += 1
    if CG_chord_old < CG_lower:
        wing_loc -= 1

    wing = Wing(wing_loc, "foam", [AR, -1, S, tr, -1, -1], [-1, 0, 0, 0], afile)
    fuse_nose = Fuselage(0, "nose", 8, [3, 3], [6, 6])

    #making the main fuselage section end where the wing ends
    main_end = wing_loc + wing.planform['cr'] - fuse_nose.length
    fuse_main = Fuselage(fuse_nose.length, "main", main_end, [6, 6])
    fuse_aft = Fuselage(fuse_main.x + fuse_main.length, "aft", 24, [6, 6], [2, 2])

    total_length = fuse_nose.length + fuse_main.length + fuse_aft.length
    NLG = LandingGear(fuse_main.x, "NLG", point_mass = True)
    MLG = LandingGear(CG + 0.1*total_length, "MLG", point_mass = True) #MLG 10% of the total length behind the CG location

    emp_xloc = fuse_aft.x + fuse_aft.length #assuming empennage plate starts where the fuselage ends
    l_emp = emp_xloc - (wing.x + wing.planform['cr']/4); #distance between c/4 of wing and c/4 of the tail
    S_h = V_h * S * wing.planform['cr'] / l_emp #that c is the mac
    S_v = V_v * S * wing.planform['b'] / l_emp

    htail = Component(emp_xloc, "balsa", type = 'area', AR = ARh, area = S_h, aero = 'hstab')
    vtail = Component(emp_xloc, "balsa", type = 'area', AR = ARv, area = S_v, aero = 'vstab')

    batt = Component(fuse_main.x, "battery", point_mass = True)
    ESC = Component(fuse_main.x, "ESC", point_mass = True)
    prop = Component(0, "prop", point_mass = True)
    motor = Component(2, "motor", point_mass = True)

    components = [fuse_nose, fuse_main, fuse_aft,
                NLG, MLG, 
                wing, htail, vtail,
                batt, ESC, prop, motor]

    Trainer = Plane("Trainer", components, Sref = S)
    Trainer.plane_plot(components) #plane plot

    CG, mtot = Trainer.CG_tally()
    MTOW = Trainer.MTOW_tally()
    CG_chord = (CG - wing.x)/ wing.planform['cr'] 
    #print(f"CG @ {CG_chord:1.1%} chord. ")

print(f"CG Placement Study done. CG @ {CG_chord:1.1%} chord; MTOW = {MTOW: 1.2f} lbs.")

print("Send to AVL for stability...")

Trainer.to_avl()
#avl- for stability analysis
ad = ADaX(planename='Trainer.avl')
ad.output_config(['t CLtot', 's Xnp', 't Elevator', 't Cref'])
#ad.use_run("AVLWrapper/AVL/Planes/Funny_Run_File.run", 1)
ad._load(); ad._oper()
ad.vcv('d1 pm 0'); ad.vcv('d2 ym 0'); ad.vcv('a a 0')
ad._x(); ad._save()
results = ad.run(print_output=False)
Xnp = float(results['Xnp']); mac = float(results['Cref']); 
plt.plot(Xnp, 0, 'x')
#print(Xnp)
#elv = float(results['Elevator'])
static_margin = sm(Xnp, CG, mac) #target this to be 5% - 10%
sm_upper = 0.10
sm_lower = 0.05
CG_chord = (CG - wing.x)/ wing.planform['cr'] #target this to be 30% - 35%

while sm_lower > static_margin < sm_upper:
    pass






#Export Results
export = open("plane_geometry.txt", 'w')
export.write(f"""
static margin: {static_margin:.1%}, CG @ {CG_chord:1.1%} chord.
MTOW estimate = {MTOW: 1.2f} lbs.\n
Fuselage Definitions: Total Length {total_length: 1.2f} in.\n
\tNose Section:\t{fuse_nose.length: 1.4f} in.\n\tMain Section:\t{fuse_main.length: 1.4f} in.\n\tAft Section:\t{fuse_aft.length: 1.4f} in.\n
Wing Definitions: (in.)\n
\tWingspan:\t{wing.planform['b']: 1.2f}\n\tRoot chord:\t{wing.planform['cr']: 1.2f}\n\tTip chord:\t{wing.planform['ct']: 1.2f}\n\tLE sweep:\t{float(90 - np.rad2deg(wing.angles['sweep'])): 1.2f} deg.\n
Empennage Definitions: (in.)\n
\thstab:\n\tchord\t{htail.c: 1.4f}\n\tspan\t{htail.b: 1.4f}\n\tvstab:\n\tchord\t{vtail.c:1.4f}\n\tspan\t{vtail.b: 1.4f}\n
Landing Gear Placement: (dist. from nose, in.)\n
\tNLG:\t{NLG.x: 1.4f}\n\tMLG:\t{MLG.x: 1.4f}
""")
export.close()
out = open("plane_geometry.txt")
print(out.read())
#plt.show()