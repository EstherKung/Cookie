#cooking
import numpy as np
import pandas as pd
from wing_def import *
from plane_def import *
from AVLWrapper.Code.Codebase import *

afile = "C:\\Users\\eneiche\\Documents\\airfoil-sauce\\airfoil_library\\Selig\\S7055.dat"

MTOW = 5 #lbs
TW = 0.225 #lb/lb
WS = 1.1 / 144 #lb/ft^2 
AR = 6; tr = 0.4; 
S =  MTOW / WS
ARh = 4; V_h = 0.40; 
ARv = 1.5; V_v = 0.04; 

wing_loc = 12; #this has to be iterated
wing = Wing(wing_loc, "foam", [AR, -1, S, tr, -1, -1], [-1, 0, 0, 0], afile)
CG = wing.x + wing.planform['cr']/2 #originally assume it's a half root chord behind the wing

fuse_nose = Fuselage(0, "nose", 3, [3, 3], [6, 6])
fuse_main = Fuselage(fuse_nose.length, "main", 24, [6, 6])
fuse_aft = Fuselage(fuse_main.x + fuse_main.length, "aft", 24, [6, 6], [3, 3])

total_length = fuse_nose.length + fuse_main.length + fuse_aft.length
NLG = LandingGear(fuse_main.x, "NLG", point_mass = True)
MLG = LandingGear(CG + 0.1*total_length, "MLG", point_mass = True) #MLG 10% of the total length behind the CG location

emp_xloc = fuse_aft.x + fuse_aft.length
l_emp = 36; #distance between c/4 of wing and c/4 of the tail
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


Trainer.plane_plot(components)
CG, mtot = Trainer.CG_tally()
print("guys guys our cg ", CG)
print("the total mass in slugs ", mtot)
MTOW = Trainer.MTOW_tally()
print("this is mtow in lbs ", MTOW)

Trainer.to_avl()

ad = ADaX(planename='Trainer.avl')
ad.output_config(['t CLtot', 's Xnp', 't Elevator', 't Cref'])
ad.use_run("AVLWrapper/AVL/Planes/Funny_Run_File.run", 1)
ad._load(); ad._oper()
ad.vcv('d1 pm 0'); ad.vcv('d2 ym 0'); ad.vcv('a a 0')
ad._x(); ad._save()
results = ad.run(print_output=True)
Xnp = float(results['Xnp'])
print(Xnp)
plt.plot(Xnp, 0, 'x')
elv = float(results['Elevator'])
static_margin = sm(Xnp, CG, float(results['Cref']))
print("elevator deflection: ", elv, " static margin: ", static_margin, " CLtot: ", float(results['CLtot']))

plt.show()