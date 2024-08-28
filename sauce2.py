#new sauce
import numpy as np
import pandas as pd
from wing_def import *
from plane_def import *
from AVLWrapper.Code.Codebase import *
from LG_def import *

def avl_np(plane: Plane, CG_chord):
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
    plt.plot(Xnp, 0, 'x', color ='tab:orange')

    return (Xnp - CG_chord)/mac #static margin percentage

def buildup(CG, MTOW, nose_sec, wing_loc, aft_sec): #tallies up the CG of the new configuration
    global wing, fuse_nose, fuse_main, fuse_aft, NLG, MLG, htail, vtail, total_length, avionics_loc
    wing = Wing(wing_loc, "foam", [AR, -1, S, tr, -1, -1], [-1, 0, 0, 0], afile)
    fuse_nose = Fuselage(0, "nose", nose_sec, [3, 3], [6, 6])

    #making the main fuselage section end where the wing ends
    main_sec = 1.6 * wing.planform['cr'] #see main_sec constraint on Notion
    fuse_main = Fuselage(fuse_nose.length, "main", main_sec, [6, 6])
    fuse_aft = Fuselage(fuse_main.x + fuse_main.length, "aft", aft_sec, [6, 6], [2, 2])

    total_length = fuse_nose.length + fuse_main.length + fuse_aft.length
    LG_height, MLG_aft, NLG_fwd_wheel, NLG_fwd_strut, MLG_side = LG_def(CG, MTOW, total_length)

    NLG = LandingGear(CG - NLG_fwd_strut, "NLG", point_mass = True)
    MLG = LandingGear(CG + MLG_aft, "MLG", point_mass = True) #MLG 10% of the total length behind the CG location

    emp_xloc = fuse_aft.x + fuse_aft.length - 0.75* htail.c #assuming empennage plate starts where the fuselage ends
    l_emp = emp_xloc - (wing.x + wing.planform['cr']/4); #distance between c/4 of wing and c/4 of the tail
    S_h = V_h * S * wing.planform['cr'] / l_emp #that c is the mac, supposedly
    S_v = V_v * S * wing.planform['b'] / l_emp

    htail = Component(emp_xloc, "balsa", type = 'area', AR = ARh, area = S_h, aero = 'hstab')
    vtail = Component(emp_xloc, "balsa", type = 'area', AR = ARv, area = S_v, aero = 'vstab')

    payload = Component(fuse_main.x + fuse_main.length/2, "payload", point_mass = True)
    avionics_loc = fuse_main.x 
    batt = Component(avionics_loc, "battery", point_mass = True)
    Rx = Component(avionics_loc, "Rx", point_mass = True)
    Rx_batt = Component(avionics_loc + 4, "Rx batt", point_mass = True)
    ESC = Component(avionics_loc, "ESC", point_mass = True)

    components = [fuse_nose, fuse_main, fuse_aft,
                NLG, MLG, 
                wing, htail, vtail,
                batt, ESC, prop, motor, 
                payload #funny payload
                ]
    

    Trainer = Plane("Trainer", components, Sref = S) #collect all component into one plane

    CG_new, mtot = Trainer.CG_tally()
    MTOW = Trainer.MTOW_tally()

    return Trainer, MTOW, CG_new

#mass convergence is still the same:
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

prop = Component(0, "prop", point_mass = True)
motor = Component(1, "motor", point_mass = True)
"----------------------------------------------------------"


#Routine Setup--------------------------------------------||
S =  MTOW / WS
    #initial fuselage definitions
nose_sec = 8; 
wing_loc = 12; 
aft_sec = 24; 
wing = Wing(wing_loc, "foam", 
            [AR, -1, S, tr, -1, -1], [-1, 0, 0, 0], afile)
#initially assume CG @ half root chord behind the wing
CG = wing.x + wing.planform['cr']/2
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
    CG_old = CG

    Trainer, MTOW, CG = buildup(CG_old, MTOW_old, nose_sec, wing_loc, aft_sec)
    #Trainer.plane_plot(Trainer.components)
    #plt.plot(CG, 0, 'x', color='tab:blue')

    res = abs(MTOW_old - MTOW)
print(f"Mass Convergence Study Done. MTOW = {MTOW: 1.2f} lbs.")
print(CG)

#CG Placement Study by varying wing placement
print("Performing Routine 2: CG Placement Study...")
CG_chord = (CG - wing.x)/ wing.planform['cr']
Trainer.plane_plot(Trainer.components)
plt.show()


#NP study

#force a CG location to iterate for NP
CG_percent_chord = 0.35
CG_chord = wing.planform['cr'] * CG_percent_chord
CG_target_loc = wing.x + CG_chord

#Configure upper & lower limits-------
sm_upper = 0.15
sm_lower = 0.10
#-------------------------------------
increment = 0.05
static_margin = avl_np(Trainer, CG_target_loc)
print(f"static margin: {static_margin: 1.2%}")

w = 0
while not (sm_lower < static_margin < sm_upper): 
    if static_margin < sm_lower: 
        aft_sec += increment
    if static_margin > sm_upper:
        aft_sec -= increment
    
    Trainer, MTOW, CG = buildup(CG_target_loc, MTOW_old, nose_sec, wing_loc, aft_sec)
    Trainer.plane_plot(Trainer.components)
    static_margin = avl_np(Trainer, CG_target_loc)
    print(f"static margin: {static_margin: 1.2%}")
    plt.plot(CG_target_loc, 0, 'x')
    plt.savefig(f"iteration/{w}.png"); w += 1; 
#plt.show()
Trainer, MTOW, CG = buildup(CG_target_loc, MTOW_old, nose_sec, wing_loc, aft_sec)

print(f"current CG discrepancy: {CG - CG_target_loc}") #if forward, negative; if aft, positive

while NLG.x < 3:
    nose_sec += 0.5
    wing_loc += 0.5
    CG_percent_chord = 0.35
    CG_chord = wing.planform['cr'] * CG_percent_chord
    CG_target_loc = wing.x + CG_chord

    Trainer, MTOW, CG = buildup(CG_target_loc, MTOW_old, nose_sec, wing_loc, aft_sec)
    #print(f"nose sec gives {nose_sec} and NLG is currently at {NLG.x}")
    Trainer.plane_plot(Trainer.components)
    #static_margin = avl_np(Trainer, CG_target_loc)
    #print(f"static margin: {static_margin: 1.2%}")
    plt.plot(CG_target_loc, 0, 'x')
    plt.savefig(f"iteration/{w}.png"); w += 1; 

    print(f"CG discrepancy {CG - CG_target_loc}")
    #print(f"confirm this stays constant over adjustment {htail.x-wing.x}")

print(MTOW)

#nose placement
#write wing section length relative to nose section
#same deal goes for the aft sec, everything should shift I think 
