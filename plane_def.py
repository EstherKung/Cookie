#class file for the plane

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Polygon
import math
from wing_def import * 

#constants at flight conditions
rho = 1.225 #standard sea level conditions (might want to be hotter than this??)
mu = 1.8e-5

class Component:
    global mat_property
    mat_property = {
        "balsa": 952.5, #g/m^2
        "cf_spar": 138.89, #g/m
        "cf_boom": 136.7, #g/m
        "styrofoam": 32000, #g/m^3
        "lwpla": 2400, #g/m^3
        "foamboard": 1500, #g/m^2
        "plywood": 2064, #g/m^2
        
        "fuse": 300, #fuselage mass
        "drivetrain": 218.6, #g // motor + esc + prop, we have 2
        "battery": 480, #g
        "avionics": 41.4, #g
        "tilt_servo": 66.0 #g
    }
    

    def __init__(this, x: float, mat: str, type = None, **kwargs):
        this.cfg = kwargs
        this.x = x #distance from the nose
        this.material = mat_property[mat]
        this.S_wet = None
        this.FF = 1 #default form factor

        match type: 
            case None: #single components: powertrain, fuselage
                this.mass = this.material
                this.x_cg = this.cfg['cg']
                this.c = this.cfg['length']
                if 'width' in this.cfg: 
                    this.b = this.cfg['width'] #for plotting
                if 'diameter' in this.cfg:
                    this.b = this.cfg['diameter'] #for plotting
                    this.d = this.cfg['diameter']; this.l = this.cfg['length'] #for drag buildup
                if 'xsec' in this.cfg:
                    this.d = 2 * (this.cfg['xsec'] / np.pi) ** 0.5
                    this.b = this.d; this.l = this.c # for drag buildup
            case 'length': #boom, spar
                this.mass = this.material * this.cfg['length']
                this.x_cg = this.cfg['xdim']/2
                this.b = this.cfg['ydim']; this.c = this.cfg['xdim'] #for plotting
            case 'area': #aero surfaces: htail, vtail
                this.mass = this.material * this.cfg['area']
                this.cfg['b'] = np.sqrt(this.cfg['area']*this.cfg['AR']) #span
                this.cfg['c'] = this.cfg['b']/this.cfg['AR'] #chord (along flow)
                this.b = this.cfg['b']; this.c = this.cfg['c']
                this.x_cg = this.cfg['c']/2
            case 'section': #to build a robust code for sections (e.g. fuse)
                this.mass = this.material * this.cfg['length']
                this.x_cg = this.cfg['x_cg']
    
    def m(this): #in g
        return this.mass
    
    def wetted_area(this, des: str, fuse_d = 0.0, d = 0.0, l = 0.0):
        match des:
            case 'flat':
                this.S_wet = 2 * this.cfg['area']
            case 'wing':
                t = this.thickness
                S_exp = this.planform['S'] - (fuse_d * this.planform['cr'])

                if this.planform['tr'] == 1:
                #this is safa formulation - rectangular wing?
                    this.S_wet = 2*(1 + 0.2*t) * S_exp
                else:
                #this is from HOOU @ HAW
                    this.S_wet = 2*S_exp*((1+0.25*t) * (1 + this.planform['tr'])/(1+this.planform['tr']))

                #this is from a paper for subsonic aircrafts
                #ar = this.cfg['AR']; tr = 1; #taper ratio
                #S = this.cfg['area']
                #kq = 1; #volume factor, for tapered wing = 0.95
                #Q = (kq * tc) * (1 + tr) ** -0.5 * S * (S/ar) ** 0.5
                #this.S_wet = (2+ 0.5*tc) * ((Q * (ar * (1+tr)) ** 0.5)/(kq * tc)) ** 2/3
            case 'nose':
                this.S_wet = 0.75 * 3.1415 * this.d * l
                #takes length of the nose
            case 'tail':
                this.S_wet = 0.72 * 3.1415 * this.d * l
                #takes length of the tail
            case 'fuse':
                this.S_wet = 2 * np.pi * this.d/2 * this.l + 2 * np.pi * (this.d/2)**2
    def FF_calc(this, des:str, l_n = 0.0, l_t = 0.0):
        match des:
            case 'wing':
                t = this.cfg['t']
                this.FF = 1 + 2 * t + 60 * t ** 4
            case 'fuse':
                l = this.l - l_n - l_t
                ld = l/this.d
                this.FF = 1 + 2.8 * (ld ** -1.5) + 3.8 * (ld ** -3)
    
    # function to calculate -> output drag buildup of component, consider flat plate
    def drag_visc(this, u, S_ref):
        global rho, mu
        Re = rho * u * this.c / mu
        c_f = 0.027/(Re**(1/7))
        c_d = c_f * this.FF * this.S_wet/S_ref
        this.cfg['c_d'] = c_d
        return c_d

    def geometry(this):
        return plt.Rectangle(xy=(this.x,-this.b/2), width=this.c, height=this.b, color='r', fill = False, linewidth = 1.5)

class Wing(Component):
    def __init__(this, x: float, mat: str, param: list, angles: list, afile: str, type = 'wing'):
        """
        this is assuming linear taper, single section wing only, constant airfoil profile throughout.
        angles in the order [sweep, dihedral, incidence, twist]
        if want c/4 no sweep, put sweep = -1
        param in the order in the order [AR, b, S, taper_ratio, c_r, c_t] 
        and for the ones that aren't defined, set to -1
        """
        super().__init__(x, mat, type = 'aero')
        this.planform = {}
        this.planform['AR'], this.planform['b'], this.planform['S'], this.planform['tr'], this.planform['cr'], this.planform['ct'] = planform(param)

        print(this.planform)
        if angles[0] == -1: this.planform['sweep'] = np.arctan(0.125*this.planform['b']*this.planform['cr']*(1-this.planform['tr'])) #for quarter no sweep
        else: this.planform['sweep'] = np.radians(angles[0]); #leading edge sweep
        this.planform['dihedral'] = angles[1]; this.planform['inc'] = angles[2]; this.planform['tw'] = angles[3]
        #to find mass
        x = np.linspace(0, this.planform['b']/2, 50) #your discretisation wooo
        xsec = csec(afile) * ((this.planform['ct'] - this.planform['cr'])/ (this.planform['b']/2) * x + this.planform['cr'])
        #print(xsec)
        this.mass = this.material * 2 * np.trapezoid(xsec, x)     
        xcentroid = centroid(afile) * ((this.planform['ct'] - this.planform['cr'])/ (this.planform['b']/2) * x + this.planform['cr'])
        this.x_cg = np.trapezoid(xcentroid, x)

        this.thickness = thickness(afile) #technically this is the thickness-to-chord, since dat file normalised it with c = 1
    
    def geometry(this):
        x = [this.x,
             0.5*this.planform['b']/2 *np.tan(this.planform['sweep']),
             0.5*this.planform['b']/2 * np.tan(this.planform['sweep']) + this.planform['ct'],
             this.x + this.planform['cr'],
             0.5*this.planform['b']/2 *np.tan(this.planform['sweep']) + this.planform['ct'],
             0.5*this.planform['b']/2 * np.tan(this.planform['sweep'])]
        y = [0,
             this.planform['b']/2, this.planform['b']/2,
             0,
             -this.planform['b']/2, -this.planform['b']/2]
        return plt.Polygon(xy = list(zip(x,y)), color='r', fill = False, linewidth = 1.5)

class Plane():
    def __init__(this, all: list, Sref):
        this.components = all
        this.Sref = Sref

    def plane_plot(this, components: list): #plot the plane: wing, hstab, vstab, fuse, boom
        plt.figure()
        plane_plot = plt.gcf()
        ax = plane_plot.gca()
        for geometry in components:
            ax.add_patch(geometry.geometry())
        ax.axis('equal')
        # set the figure limits
        plt.xlim([-0.5, 2]); plt.ylim([-2, 2])
    
    ### this section is for MTOW tally
    def MTOW_tally(this):
        m_tot = 0
        for component in this.components:
            m_tot = m_tot + component.mass
        MTOW = 9.81 * m_tot/1000
        return MTOW #Newtons, N

    ### this section is for CG tally
    def CG_tally(this):
        m_tot = 0 
        tally = 0 
        for component in this.components:
            loc = component.x_cg + component.x
            tally = tally + (component.mass * loc)
            m_tot = m_tot + (component.mass)
        cg = tally/m_tot
        m_tot = m_tot/1000 #convert to kg
    
        return[cg, m_tot]

    ### this section is for drag buildup - wip
    def drag_tally(this, external, u):
        global rho
        cd_tot = 0
        for component in external:
            cd_tot = cd_tot + component.drag_visc(u, this.Sref)
        D = 0.5 * rho * u **2 * this.Sref * cd_tot
        return D
    
    def to_avl():
        pass

me = Wing(0, "styrofoam", [6, -1, 0.2, 0.4, -1, -1], [0, 0, 0, 0], "NACA0012.dat")
print(me.mass)
print(me.x_cg)