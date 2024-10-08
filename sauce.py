#cooking
import numpy as np
import pandas as pd
from wing_def import *
from plane_def import *
from AVLWrapper.Code.Codebase import *
from LG_def import *
from powerplant_def import *
 
w = 0; 
def avl_np(plane: Plane, CG_chord):
    #for avl stability analysis
    plane.to_avl()
    ad = ADaX(planename='stick.avl')
    ad.output_config(['s Xnp', 't Cref'])
    #ad.use_run("AVLWrapper/AVL/Planes/Funny_Run_File.run", 1)
    ad._load(); ad._oper()
    ad.vcv('d1 pm 0'); ad.vcv('d2 ym 0'); ad.vcv('a a 0')
    ad._x(); ad._save()
    results = ad.run(print_output=False)
    Xnp = float(results['Xnp']); mac = float(results['Cref']); 
    plt.plot(Xnp, 0, 'x', color ='tab:orange')

    return (Xnp - CG_chord)/mac #static margin percentage

def buildup(CG, MTOW, wing_def: list, wing_loc, l_h, S_h = None): #tallies up the CG of the new configuration
    global wing, boom, pod, NLG, MLG, htail, vtail, total_length, payload

    ### 1. Wing
    AR, b, S, tr, cr, ct = planform(wing_def)
    
    wing = Wing(wing_loc, "other", 
                [AR, -1, S, tr, -1, -1], [[-1, -1], 0, 0, 0], 
                afile, MTOW)

    ### 2. EMPENNAGE
    stagger = 0.5 #gap distance between horizontal & vertical plate
    if S_h is None: S_h = V_h * S * wing.planform['cr'] / l_h #that c is typ. MAC; DBF 2025 uses rectangular wing
    emp_datum = wing.x + wing.MAC()/4 + l_h

    htail = Component(emp_datum, "emp", type = 'area', 
                      AR = ARh, area = S_h, aero = 'hstab')
    
    l_v = l_h + htail.c + stagger
    S_v = V_v * S * wing.planform['b'] / l_v
    vtail = Component(emp_datum + htail.c + 0.5, "emp", type = 'area', 
                      AR = ARv, area = S_v, aero = 'vstab')

    ### 3. POD + BOOM
    pod = Component(2.4, "pod", length = 25, diameter = 3, x_cg = 0.4* 25)
    total_length = vtail.x + vtail.c #total length of plane
    boom = Component(pod.c, "boom", type = 'length', 
                     length = total_length - pod.c,
                     xdim = total_length - pod.c, ydim = 1)

    ### 4. Landing Gear
    LG_height, MLG_aft, NLG_fwd_wheel, NLG_fwd_strut, MLG_side = LG_def(CG, MTOW, total_length)

    # NLG = LandingGear(CG - NLG_fwd_strut, "NLG", point_mass = True)
    # MLG = LandingGear(CG + MLG_aft, "MLG", point_mass = True) 

    # LG - initial guesses
    NLG = LandingGear(CG - NLG_fwd_strut, "other", type = 'specified', 
                      mass = MTOW * 0.1 * 0.35, x_cg = 0)
    MLG = LandingGear(CG + MLG_aft, "other", type = 'specified', 
                      mass = MTOW * 0.1 * 0.65, x_cg = 0)
    
    ### Propulsion 
    motor = Component(1.5, "motor", point_mass = True)
    prop = Component(0.25, "prop", point_mass=True)
    
    ### Avionics
    batt = Component(8.5, "battery", point_mass = True)
    Rx = Component(7, "Rx", point_mass = True)
    Rx_batt = Component(7, "Rx batt", point_mass = True)
    ESC = Component(6, "ESC", point_mass = True)
    
    components = [boom, pod, 
                NLG, MLG, 
                wing, htail, vtail,
                batt, ESC, Rx, Rx_batt,
                motor, prop]
    

    stick_empty = Plane("stick, empty weight", components, Sref = S) #collect all component into one plane
    W_empty = stick_empty.MTOW_tally()
    # print(f"look {W_empty}, look {MTOW}")

    payload = Component(wing.x + wing.planform['cr']/4, "other",
                        type = "specified", mass = MTOW - W_empty, point_mass = True)
    
    print(f"payload mass: {payload.m(): 1.2f} lbs")
    
    components.append(payload)
    print(f"Payload Mass Fraction: {payload.m()/MTOW: 1.1%}")
    
    stick = Plane("stick", components, Sref = S)

    CG_new, M_tot = stick.CG_tally()
    MTOW = stick.MTOW_tally()

    CG_chord = (CG - wing.x)/ wing.planform['cr']
    #print(f"CG @ {CG_chord:1.1%} chord. ")
    return stick, MTOW, CG_new

afile = "airfoil_library/AG25.dat"

"Point on Constraint Diagram-------------------------------"
TW = 0.20                                          #lb/lb--
T = 10 #lbs
WS = 2            / 144 #lb/ft^2 * (ft/12in)^2 = lb / in^2
print(
    f"""Point Selected from the Constraint Diagram: 
    TW = {TW} lb/lb, WS = {WS: 1.4f} lb/sq.in"""
    )
"----------------------------------------------------------"

"Tail Volume Coefficients & AR-----------------------------"
AR = 5.5; tr = 0.4; 
ARh = 0.5*AR; V_h = 0.45; 
ARv = 1.5; V_v = 0.04; 
"----------------------------------------------------------"

"Specifications--------------------------------------------"
MTOW = WS * (72 ** 2 / AR)                              #lbs
print(f"initial MTOW = {MTOW: 1.2f} lbs.")
"----------------------------------------------------------"

#Routine Setup--------------------------------------------||
wing_loc = 28
wing = Wing(wing_loc, "other", 
                [AR, 72, -1, tr, -1, -1], [[-1, -1], 0, 0, 0], 
                afile, MTOW)
#initially assume CG @ half root chord behind the wing
#CG = wing.x + wing.planform['cr']/2
#initially assume l_h is 40% span
l_h = 0.4 * wing.planform['b']
S_h = V_h * wing.planform['S'] * wing.planform['cr'] / l_h 
#--------------------------------------------------------||


#force a CG location to iterate for NP
CG_percent_chord = 0.35
CG_chord = wing.planform['cr'] * CG_percent_chord
CG_target_loc = wing.x + CG_chord

stick, MTOW, CG = buildup(wing.x + CG_chord, MTOW, [AR, 72, -1, tr, -1, -1], wing_loc, l_h)
print(MTOW)
plt.show()

#Configure upper & lower limits-------
sm_upper = 0.15
sm_lower = 0.08
#-------------------------------------
increment = 2
static_margin = avl_np(stick, CG_target_loc)
print(f"static margin: {static_margin: 1.2%}")

stick.plane_plot(stick.components)
static_margin = avl_np(stick, CG_target_loc)
print(f"static margin: {static_margin: 1.2%}")
plt.plot(CG_target_loc, 0, 'x')
plt.savefig(f"iteration/{w}.png"); w += 1; 

w = 0
while not (sm_lower < static_margin < sm_upper): 
    MTOW_old = MTOW
    if static_margin < sm_lower: 
        S_h -= increment
    if static_margin > sm_upper:
        #l_h -= increment
        S_h += increment
    print(l_h)
    stick, MTOW, CG = buildup(CG_target_loc, MTOW_old, [AR, 72, -1, tr, -1, -1], wing_loc, l_h, S_h)
    stick.plane_plot(stick.components)
    static_margin = avl_np(stick, CG_target_loc)
    print(f"static margin: {static_margin: 1.2%}")
    plt.plot(CG_target_loc, 0, 'x')
    plt.savefig(f"iteration/{w}.png"); w += 1; 
#plt.show()
stick, MTOW, CG = buildup(CG_target_loc, MTOW, [AR, 72, -1, tr, -1, -1], wing_loc, l_h)

print(f"current CG discrepancy: {CG - CG_target_loc}") #if forward, negative; if aft, positive
print(f"boom length: {boom.c}")
print(f"current Aspect Ratio: {AR}")
print(stick.mass_breakdown())

#Export Results
export = open("plane_geometry.txt", 'w')
export.write(f"""[Results]
static margin: {static_margin:.1%} forcing CG @ 35% chord. CG discrepancy: {CG - CG_target_loc} in.\n
MTOW estimate = {MTOW: 1.2f} lbs. Payload Mass Fraction: {payload.m()/MTOW: 1.1%}\n
Planform Area S = {wing.planform['S']: 1.4f} AR = {AR}\n
Available Thrust: {T: 1.2f} lbs.\n
Fuselage Definitions: Total Length {total_length: 1.2f} in. Total Mass {boom.m() + pod.m(): 1.4f} lbs.\n
\tPod Length: {pod.c: 1.2f} in.\n\tBoom Length: {boom.c: 1.2f} in.\n
Wing Definitions: Location from nose: {wing.x:1.2f} (in.)\n
\tWingspan:\t{wing.planform['b']: 1.2f}\n\tRoot chord:\t{wing.planform['cr']: 1.2f}\n\tTip chord:\t{wing.planform['ct']: 1.2f}\n\tLE sweep:\t{float(np.rad2deg(wing.angles['sweep'])): 1.2f} deg.\n
Empennage Definitions: (in.)\n
\thstab:\n\tchord\t{htail.c: 1.4f}\n\tspan\t{htail.b: 1.4f}\n\tloc\t{htail.x: 1.4f}\n
\tvstab:\n\tchord\t{vtail.c:1.4f}\n\tspan\t{vtail.b: 1.4f}\n\tloc\t{vtail.x:1.4f}\n
Landing Gear Placement: (dist. from nose, in.)\n
\tNLG:\t{NLG.x: 1.4f}\n\tMLG:\t{MLG.x: 1.4f}
""")
export.close()
# out = open("plane_geometry.txt")
# print(out.read())