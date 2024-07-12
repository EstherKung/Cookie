from Codebase import *


def aircraft_geometry(fname: str, S, ar, tr, sw = None, dh = 0.0, inc = 0.0, tw = 0.0):
    #Planform parameters for defining your main wing
#    S = 0.2
#    ar = 6
#    tr = 1

    #Calculating values more relevant to AVL
    b = np.sqrt(S*ar)
    cr = 2*S/(b*(1+tr))
    ct = tr*cr

    #Four angles for defining your wing
    if sw is None: sweep = np.arctan(0.125*b*cr*(1-tr)) #for quarter no taper
    else: sweep = sw /180*np.pi #leading edge
    dihedral = dh /180*np.pi
    incidence = inc 
    twist = tw
    
    "Plane"
    plane = des_plane(Name= fname, Mach=0, Sref = S, bref = 0, cref = 0, Xref=0, Yref=0.0, Zref=0)

    "Surfaces"
    af1 = 'AVL/Airfoils/SD7032-099-88.dat'
    wing = des_surf(plane, Name='Wing', YDUPLICATE=0.0, Nspanwise=31)
    hstab = des_surf(plane, Name='Hstab', TRANSLATE=(0.720, 0, 0.19), YDUPLICATE=0.0)
    vstab = des_surf(plane, Name='Vstab', TRANSLATE=(0.700, 0, 0))
    plane.use_reference(wing)
    hstab.assign(COMPONENT=1)
    vstab.assign(COMPONENT=1)

    "Sections"
    wingsec000 = des_sec(wing, Xle=0.0, Yle=0.0, Zle=0.0, Chord=cr, Ainc=incidence, AFILE=af1)
    wingsec200 = des_sec(wing, Xle=b/2*np.tan(sweep), Yle=b/2, Zle=b/2*np.tan(dihedral), Chord=ct, Ainc=incidence+twist, AFILE=af1) 

    hstabsec000 = des_sec(hstab, Xle=0.0, Yle=0.0, Zle=0.0, Chord=0.07)
    hstabsec100 = des_sec(hstab, Xle=0, Yle=0.138, Zle=0.0, Chord=0.07)
    vstabsec000 = des_sec(vstab, Xle=0.0, Yle=0.0, Zle=0.0, Chord=0.10)
    vstabsec100 = des_sec(vstab, Xle=0.0, Yle=0, Zle=0.175, Chord=0.10)

    "Control Surfaces"
    #Keep your control surfaces named "Aileron", "Elevator", and "Rudder" for the best results
    elvstart = des_ctrl(hstabsec000, Cname='Elevator', Xhinge=0.6, SgnDup=1)
    elvend = des_ctrl(hstabsec100, Cname='Elevator', Xhinge=0.6, SgnDup=1)
    rudstart = des_ctrl(vstabsec000, Cname = 'Rudder', Xhinge = 0.75)
    rudend = des_ctrl(vstabsec100, Cname = 'Rudder', Xhinge = 0.75)
    """
    Don't need to edit this unless necessary
    """
    wing.order(orientation=HORIZONTAL)
    hstab.order(orientation=HORIZONTAL)
    vstab.order(orientation=VERTICAL)

    plane.write_to_avl(filename="Junitaper")

def tail_zerotrim(mac, S_h, V_ht = 0.4, V_vt = 0.07, l = None): #manually write in angles & wing def
    #Planform parameters for defining your main wing
    S = 0.2
    ar = 6
    tr = 1
    #need to know mean aerodynamic chord, so pull this from the original geometry first 't Cref's

    #Calculating values more relevant to AVL
    b = np.sqrt(S*ar)
    cr = 2*S/(b*(1+tr))
    ct = tr*cr

    #Four angles for defining your wing
    sweep = np.arctan(0.125*b*cr*(1-tr)) #for quarter no taper
    dihedral = 0 /180*np.pi
    incidence = 0 
    twist = 0


    #htail parameters
    tr_h = 1; ARh = 4
    sw_h = np.radians(0) #LE sweep

    b_h = np.sqrt(ARh * S_h)/2
    cr_h =S_h/(b_h * (1 + tr))
    ct_h = cr_h * tr_h
    xle_h = b_h * np.tan(sw_h)

    if l is None: l = (V_ht * S * mac)/S_h

    S_v = (V_ht * S * b)/l

    #vtail parameters
    tr_v = 1; ARv = 1.2
    sw_v = np.radians(0) #LE sweep

    b_v = np.sqrt(ARv *S_v)/2
    cr_v = S_v * tr_v
    ct_v = cr_v * tr_v
    xle_v = b_v * np.tan(sw_v)

    "Plane"
    plane = des_plane(Name='Juni', Mach=0, Sref = 0.2, bref = 0, cref = 0, Xref=0, Yref=0.0, Zref=0)

    "Surfaces"
    af1 = 'AVL/Airfoils/SD7032-099-88.dat'
    wing = des_surf(plane, Name='Wing', YDUPLICATE=0.0, Nspanwise=31)
    #this is for a t-tail
    hstab = des_surf(plane, Name='Hstab', TRANSLATE=(l, 0, b_v), YDUPLICATE=0.0)
    vstab = des_surf(plane, Name='Vstab', TRANSLATE=(l, 0, 0))
    plane.use_reference(wing)
    hstab.assign(COMPONENT=1)
    vstab.assign(COMPONENT=1)

    "Sections"
    wingsec000 = des_sec(wing, Xle=0.0, Yle=0.0, Zle=0.0, Chord=cr, Ainc=incidence, AFILE=af1)
    wingsec200 = des_sec(wing, Xle=b/2*np.tan(sweep), Yle=b/2, Zle=b/2*np.tan(dihedral), Chord=ct, Ainc=incidence+twist, AFILE=af1) 

    hstabsec000 = des_sec(hstab, Xle=0.0, Yle=0.0, Zle=0.0, Chord=cr_h)
    hstabsec100 = des_sec(hstab, Xle=xle_h, Yle=b_h, Zle=0, Chord=ct_h)
    vstabsec000 = des_sec(vstab, Xle=0.0, Yle=0.0, Zle=0.0, Chord=cr_v)
    vstabsec100 = des_sec(vstab, Xle=xle_v, Yle=0.0, Zle=b_v, Chord=ct_v)

    "Control Surfaces"
    #Keep your control surfaces named "Aileron", "Elevator", and "Rudder" for the best results
    elvstart = des_ctrl(hstabsec000, Cname='Elevator', Xhinge=0.75, SgnDup=1)
    elvend = des_ctrl(hstabsec100, Cname='Elevator', Xhinge=0.75, SgnDup=1)
    rudstart = des_ctrl(vstabsec000, Cname = 'Rudder', Xhinge = 0.75)
    rudend = des_ctrl(vstabsec100, Cname = 'Rudder', Xhinge = 0.75)
    """
    Don't need to edit this unless necessary
    """
    wing.order(orientation=HORIZONTAL)
    hstab.order(orientation=HORIZONTAL)
    vstab.order(orientation=VERTICAL)

    plane.write_to_avl(filename='Junitaper')

def sm(Xnp, Xcg, mac):
    sm = (Xnp - Xcg)/mac
    return sm