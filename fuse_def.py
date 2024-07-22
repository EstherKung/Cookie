#material properties used across all sections
global maFG, mlCFTow, denASA
maFG = 2.74 * (1/514.785) * ((1/36)**2); # Mass per unit area oz/yd**2 to slug/in **2
mlCFTow = 0.446 * (1/14593.9) * (1/39.37) # Mass/Length of carbon fiber tow g/m to slug/in
infill_factor = 0.6; #assume this is the correction factor (can verify with the prints or check what the)
denASA = infill_factor* (.99) * (1/14593.9) * ((2.54/1)**3); # g/cm**3 to slug/in**3
#from Dean, translated 7/21

ASA_thickness = 3/32; #in, thickness of asa aero print
glass_thickness = 0.012; #in, two layers of fiberglass

def main_fuse(mfB, mfH, mfL):
    """
    Main Fuselage
    inputs: cross section base, height, section length (in)
    outputs: mass (slugs), cg (in)
    """

    # Geometric properties of main fuselage
    mfb = mfB - (2*(ASA_thickness + glass_thickness)); # Inner base length, in 
    mfh = mfH - (2*(ASA_thickness + glass_thickness)); # Outer base length, in
    mfCSA = (mfB*mfH)-(mfb*mfh); # Cross sectional area of main fuselage

    mfIyy = ((mfB*(mfH**3))/12) - ((mfb*(mfh**3))/12) # Moment of inertia in**4
    mfIzz = ((mfH*(mfB**3))/12) - ((mfh*(mfb**3))/12) # Moment of inertia in**4
    mfryy = (mfIyy/(mfCSA))**.5 # Radius of gyration in
    mfrzz = (mfIzz/(mfCSA))**.5 # Radius of gyration in

    # Main fuselage center of gravity
    mfcgx = mfB/2; mfcgy = mfH/2; mfcgz = mfL/2

    # Code to weigh main fuselage directly
    #mfwASA = 0; # Weight of 3D printed main fuselage in g measured in scale
    #mfwASAslug = mfwASA * (1/14593.9); # Conversion of g to slugs 
    #mfmlASA = mfwASAslug/mfL # Mass/Length of ASA slug/in

    # Code for estimated main fuselage weight
        #pure ASA
    mf5inwASA = 55.8 * (1/14593.9); # Experimental weight of main fuselage 5 in test section g to slugs
    mfobASA = mfB - glass_thickness; # Outer base length of ASA aero fuselage
    mfohASA = mfH - glass_thickness; # Outer height of ASA aero fuselage
    mfibASA = mfB - glass_thickness - (2*ASA_thickness); # Inner base length of ASA aero fuselage
    mfihASA = mfH - glass_thickness - (2*ASA_thickness); # Inner height of ASA aero fuselage
    mfCSAASA = (mfobASA*mfohASA)-(mfibASA*mfihASA); # Cross sectional area of ASA print
    mfASAV = mfCSAASA * 5; # Volume of ASA main fuselage 5 in test section in**3
    denASA = mf5inwASA/mfASAV; # Experimental density of main fuselage 5 in test section slug/in**3
    mfASAV2 = mfCSAASA * mfL; # Volume of main fuselage
    mfwASA = denASA * mfASAV2 # Estimated weight of main fuselage print slugs
    mfmlASA = mfwASA/mfL # Estimated mass/length of main fuselage print slug/in

    # Material properties of carbon fiber tow
    mlCFTow = 0.446 * (1/14593.9) * (1/39.37) # Mass/Length of carbon fiber tow g/m to slug/in
    mfwCFTow = mlCFTow *mfL # Weight of carbon fiber tow slugs

    # Material properties of fiber glass main fuselage
    mfoSArl = 2 * (mfohASA*mfL); # Outer surface area 2 sides right and left 
    mfoSAtb = 2 * (mfobASA*mfL); # Outer surface area 2 sides top and bottom
    mfiSArl = 2 * (mfihASA*mfL); # Inner surface area 2 sides right and left 
    mfiSAtb = 2 * (mfibASA*mfL); # Inner surface area 2 sides top and bottom
    mfSAtot = mfoSArl + mfoSAtb + mfiSArl + mfiSAtb; # Total surface area
    maFG = 2.74 * (1/514.785) * ((1/36)**2); # Mass per unit area oz/yd**2 to slug/in **2
    mfwFG = mfSAtot * maFG # Weight of 2 layers of fiber glass
    mfmlFG = mfwFG/mfL # Mass/unit length of fiber glass

    # Material properties of Aero expoxy
    mfwAE = mfwFG # Weight of Aero Epoxy (50:50 ratio of fiber glass and aero epoxy)
    mfmlAE = mfmlFG # Mass/unit length of Aero Epoxy (50:50 ratio of fiber glass and aero epoxy)

    # Total weight and mass/length of main fuselage
    mlMainFuselageslugin = mfmlASA + (4*mlCFTow) + (mfmlFG + mfmlAE) # Main fuselage mass/length slug/in
    mlMainFuselagelbsin = mlMainFuselageslugin * 32.17 # Main fuselage mass/length lbs/in
    wMainFuselageslug = mfwASA + (4*mfwCFTow) + (mfwFG +mfwAE) # Main fuselage weight slug
    wMainFuselagelbs = wMainFuselageslug * 32.17 # Main fuselage weight lbs
    
    return wMainFuselageslug, mfcgx

def nose_fuse(nB1, nH1, nB2, nH2, nL):
    """
    Nose Section
    inputs: face 1 base & height, face 2 base & height, length (in)
    outputs: mass (slugs), cg (in) 
    """

    nA1 = nB1 * nH1; # Area of face 1 in**2
    nA2 = nB2 * nH2; # Area of face 2 in**2
    nAa = (nA2 + nA1)/2; # Average area of face 1 and 2 in**2
    noSL = nAa**.5; # Nose outer side length assuming average cross sectional shape is a square in
    niSL = noSL - (2*(ASA_thickness + glass_thickness)); # Nose inner side length assuming average cross sectional shape is a square in
    nCSA = (noSL**2) - (niSL**2); # Nose cross sectional area in**2

    nIyy = ((noSL**4)/12) - ((niSL**4)/12) # Moment of inertia in**4
    nIzz = ((noSL**4)/12) - ((niSL**4)/12) # Moment of inertia in**4
    nryy = (nIyy/(nCSA))**.5 # Radius of gyration in
    nrzz = (nIzz/(nCSA))**.5 # Radius of gyration in

    # Nose fuselage center of gravity
    ncgx = noSL/2
    ncgy = noSL/2
    ncgz = nL/2

    # Material properties of ASA Aero

    # # Code to weigh nose fuselage directly
    # nwASA = 0; # Weight of 3D printed nose fuselage in g measured in scale
    # nwASAslug = nwASA * (1/14593.9); # Conversion of g to slugs 
    # nmlASA = nwASAslug/nL # Mass/Length of ASA slug/in

    # Code for estimated nose fuselage weight
    noSLASA = noSL - glass_thickness; # Outer side length of ASA aero fuselage
    niSLASA = noSL - glass_thickness - (2*ASA_thickness); # Inner side length of ASA aero fuselage
    nCSAASA = (noSLASA**2)-(niSLASA**2); # Cross sectional area of ASA print
    nASAV = nCSAASA * nL; # Volume of nose fuselage
    nwASA = denASA * nASAV # Estimated weight of nose fuselage print slugs
    nmlASA = nwASA/nL # Estimated mass/length of nose fuselage print slug/in

    # Material properties of carbon fiber tow
    nwCFTow = mlCFTow *nL # Weight of carbon fiber tow slugs

    # Material properties of fiber glass nose fuselage
    noSA = 4 * (noSLASA*nL); # Outer surface area for nose ASA print in**2
    niSA = 4 * (niSLASA*nL); # Inner surface area for nose ASA print in**2
    nSAtot = noSA + niSA; # Total surface area in**2
    nwFG = nSAtot * maFG # Weight of 2 layers of fiber glass
    nmlFG = nwFG/nL # Mass/unit length of fiber glass

    # Material properties of Aero expoxy
    nwAE = nwFG # Weight of Aero Epoxy (50:50 ratio of fiber glass and aero epoxy)
    nmlAE = nmlFG # Mass/unit length of Aero Epoxy (50:50 ratio of fiber glass and aero epoxy)

    # Total weight and mass/length of main fuselage
    mlNoseslugin = nmlASA + (4*mlCFTow) + (nmlFG + nmlAE) # Nose fuselage mass/length slug/in
    mlNoselbsin = mlNoseslugin * 32.17 # Nose fuselage mass/length lbs/in
    wNoseslug = nwASA + (4*nwCFTow) + (nwFG +nwAE) # Nose fuselage weight slug
    wNoselbs = wNoseslug * 32.17 # Nose fuselage weight lbs

    return wNoseslug, ncgx

def aft_fuse(afB1, afH1, afB2, afH2, afL):
    """
    Aft Section
    inputs: face 1 base & height, face 2 base & height, length (in)
    outputs: mass (slugs), cg (in) 
    """

    # Geometric properties of aft fuselage
    afA1 = afB1 * afH1; # Area of face 1 in**2
    afA2 = afB2 * afH2; # Area of face 2 in**2
    afAa = (afA2 + afA1)/2; # Average area of face 1 and 2 in**2
    afoSL = afAa**.5; # Aft outer side length assuming average cross sectional shape is a square in
    afiSL = afoSL - (2*(ASA_thickness + glass_thickness)); # Aft inner side length assuming average cross sectional shape is a square in
    afCSA = (afoSL**2) - (afiSL**2); # Aft cross sectional area in**2
    afIyy = ((afoSL**4)/12) - ((afiSL**4)/12) # Moment of inertia in**4
    afIzz = ((afoSL**4)/12) - ((afiSL**4)/12) # Moment of inertia in**4
    afryy = (afIyy/(afCSA))**.5 # Radius of gyration in
    afrzz = (afIzz/(afCSA))**.5 # Radius of gyration in

    # Nose fuselage center of gravity
    afcgx = afoSL/2
    afcgy = afoSL/2
    afcgz = afL/2

    # Material properties of ASA Aero

    # Code to weigh aft fuselage directly
    afwASA = 0; # Weight of 3D printed aft fuselage in g measured in scale
    afwASAslug = afwASA * (1/14593.9); # Conversion of g to slugs 
    afmlASA = afwASAslug/afL # Mass/Length of ASA slug/in

    # Code for estimated aft fuselage weight
    afoSLASA = afoSL - glass_thickness; # Outer side length of ASA aero fuselage
    afiSLASA = afoSL - glass_thickness - (2*ASA_thickness); # Inner side length of ASA aero fuselage
    afCSAASA = (afoSLASA**2)-(afiSLASA**2); # Cross sectional area of ASA print
    afASAV = afCSAASA * afL; # Volume of aft fuselage
    afwASA = denASA * afASAV # Estimated weight of aft fuselage print slugs
    afmlASA = afwASA/afL # Estimated mass/length of aft fuselage print slug/in

    # Material properties of carbon fiber tow
    afwCFTow = mlCFTow *afL # Weight of carbon fiber tow slugs

    # Material properties of fiber glass aft fuselage
    afoSA = 4 * (afoSLASA*afL); # Outer surface area for aft ASA print in**2
    afiSA = 4 * (afiSLASA*afL); # Inner surface area for aft ASA print in**2
    afSAtot = afoSA + afiSA; # Total surface area in**2
    afwFG = afSAtot * maFG # Weight of 2 layers of fiber glass
    afmlFG = afwFG/afL # Mass/unit length of fiber glass

    # Material properties of Aero expoxy
    afwAE = afwFG # Weight of Aero Epoxy (50:50 ratio of fiber glass and aero epoxy)
    afmlAE = afmlFG # Mass/unit length of Aero Epoxy (50:50 ratio of fiber glass and aero epoxy)

    # Total weight and mass/length of main fuselage
    mlAftslugin = afmlASA + (4*mlCFTow) + (afmlFG + afmlAE) # Aft fuselage mass/length slug/in
    mlAftlbsin = mlAftslugin * 32.17 # Aft fuselage mass/length lbs/in
    wAftslug = afwASA + (4*afwCFTow) + (afwFG +afwAE) # Aft fuselage weight slug
    wAftlbs = wAftslug * 32.17 # Aft fuselage weight lbs

    return wAftslug, afcgx