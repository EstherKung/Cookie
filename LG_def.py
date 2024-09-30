import numpy as np

def LG_def(CG, MTOW, fuse_len,
           tipback_angle = 15, AOA_TO = 15, overturn_angle = 60,
           NLG_W_percent = 0.10, NLG_strut_travel = 7):
    #Fuselage Characteristics!
    Fuselage_len = fuse_len /12 #ft, converted from in
    CG = CG /12 #ft, converted from in

    deg2rad = np.pi/180

    LG_height = (np.tan(AOA_TO*deg2rad)*(Fuselage_len-CG))/(1+(np.tan(AOA_TO*deg2rad)*np.tan(tipback_angle*deg2rad)))   #for both MLG and NLG

    MLG_aft = np.tan(tipback_angle*deg2rad)*LG_height   #aft from CG

    NLG_fwd_wheel = MLG_aft*((MTOW-MTOW*NLG_W_percent)/(MTOW*NLG_W_percent))   #Fwd from CG wheel location

    NLG_fwd_strut = NLG_fwd_wheel - (LG_height*np.tan(NLG_strut_travel*deg2rad))    #Fwd from CG to strut location 

    MLG_side = (NLG_fwd_wheel+MLG_aft)*np.tan((np.arcsin(((LG_height/np.tan(overturn_angle*deg2rad))/NLG_fwd_wheel))))

    #also convert everything to inches again for output
    return LG_height *12, MLG_aft *12, NLG_fwd_wheel *12, NLG_fwd_strut *12, MLG_side *12

    #all these locations are relevant to the CG