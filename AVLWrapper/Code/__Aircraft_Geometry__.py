from Code.Codebase import *

#Planform parameters for defining your main wing
S = 8890
ar = 6.89
tr = 0.4
#Four angles for defining your wing
sweep = 30 /180*np.pi
dihedral = -3 /180*np.pi
incidence = 0 
twist = -3
#Calculating values more relevant to AVL
b = np.sqrt(S*ar)
cr = 2*S/(b*(1+tr))
ct = tr*cr

"Plane"
plane = des_plane(Name='TitanWing', Mach=0.8, Sref = 0, bref = 0, cref = 0, Xref=45.478, Yref=0.0, Zref=4.83)

"Surfaces"
af1 = 'AVL/Airfoils/sc20710.dat'
af2 = 'AVL/Airfoils/sc20012.dat'
wing = des_surf(plane, Name='Wing', YDUPLICATE=0.0, Nspanwise=31)
hstab = des_surf(plane, Name='Hstab', TRANSLATE=(141.75,0,15), YDUPLICATE=0.0)
vstab = des_surf(plane, Name='Vstab', TRANSLATE=(131.5,0.0,0.0))
plane.use_reference(wing)
hstab.assign(COMPONENT=1)
vstab.assign(COMPONENT=1)

"Sections"
wingsec000 = des_sec(wing, Xle=0.0, Yle=0.0, Zle=0.0, Chord=cr, Ainc=incidence, AFILE=af1)
wingsec200 = des_sec(wing, Xle=b/2*np.tan(sweep), Yle=b/2, Zle=b/2*np.tan(dihedral), Chord=ct, Ainc=incidence+twist, AFILE=af1) 

hstabsec000 = des_sec(hstab, Xle=0.0, Yle=0.0, Zle=0.0, Chord=28.8, AFILE=af2)
hstabsec100 = des_sec(hstab, Xle=31.50934, Yle=45, Zle=0.0, Chord=7.2, AFILE=af2)
hstabsec010 = des_sec.interp(s1=hstabsec000, s2=hstabsec100, span=2.5)
hstabsec090 = des_sec.interp(s1=hstabsec000, s2=hstabsec100, span=40.5)

vstabsec000 = des_sec(vstab, Xle=0.0, Yle=0.0, Zle=0.0, Chord=35.00161, AFILE=af2)
vstabsec100 = des_sec(vstab, Xle=27.30110, Yle=0, Zle=42.04, Chord=19.25089, AFILE=af2)

"Control Surfaces"
#Keep your control surfaces named "Aileron", "Elevator", and "Rudder" for the best results
elvstart = des_ctrl(hstabsec010, Cname='Elevator', Xhinge=0.80, SgnDup=1)
elvend = des_ctrl(hstabsec090, Cname='Elevator', Xhinge=0.80, SgnDup=1)

"""
Don't need to edit this unless necessary
"""
wing.order(orientation=HORIZONTAL)
hstab.order(orientation=HORIZONTAL)
vstab.order(orientation=VERTICAL)

plane.write_to_avl(filename='TitanWing')
