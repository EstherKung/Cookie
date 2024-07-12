from Codebase import *
from aircraft_sauce import *
aircraft_geometry('Juni', 0.24, 5.75, 1, inc = 1.5)

#ad = ADaX(planename='Juni.avl')
# ad.output_config(['s Xnp', 't Elevator', 't Cref'])
# ad._load(); ad._oper()
# ad.vcv('d1 pm 0'); ad.vcv('d2 ym 0'); ad.vcv('a a 0')
# ad._x(); ad._save()

# init = ad.run(print_output=False)

# Xnp = float(init['Xnp'])
# elv = float(init['Elevator'])
# mac = float(init['Cref'])

Xcg = 0.05

#tail_zerotrim(mac, 0.03, V_ht = 0.43, l = 0.77)
ad = ADaX(planename='Junitaper.avl')
ad.output_config(['t CLtot', 's Xnp', 't Elevator', 't Cref'])
ad._load(); ad._oper()
ad.vcv('d1 pm 0'); ad.vcv('d2 ym 0'); ad.vcv('a a 0')
ad._x(); ad._save()
#ad.run(print_output = True)
results = ad.run(print_output=False)
Xnp = float(results['Xnp'])
elv = float(results['Elevator'])
static_margin = sm(Xnp, Xcg, 0.2043)
print("elevator deflection: ", elv, " static margin: ", static_margin, " CLtot: ", float(results['CLtot']))

while not -0.5 <= elv <= 0.5 and not 0.10 < static_margin < 0.15:
    break
