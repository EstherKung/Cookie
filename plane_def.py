#class file for the plane

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Polygon
from wing_def import * 
from AVLWrapper.aircraft_sauce import *
from material_library import *
from fuse_def import *

#constants at flight conditions
rho = 0; #standard sea level conditions (might want to be hotter than this??)
mu = 0; #need these in imperial

class Component:
    def __init__(this, x: float, mat: str, type = None, point_mass = False, **kwargs):
        this.cfg = kwargs 
        #specify aero = wing, hstab, vstab for avl purposes
        this.x = x #distance from the nose
        
        this.material = materials(mat)
        this.S_wet = None
        this.FF = 1 #default form factor

        match type: 
            case None: #single components: powertrain, fuselage
                this.mass = this.material
                if point_mass: this.x_cg = 0; this.c = 0; this.b = 0; 
                else:
                    this.c = this.cfg['length']
                    if not 'cg' in this.cfg: this.x_cg = this.c/2
                    else: this.x_cg = this.cfg['cg']
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
                l = this.l - l_n - l_t #maybe we don't need to subtract?
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
        return plt.Rectangle(xy=(this.x,-this.b/2), width=this.c, height=this.b, color='tab:purple', fill = False, linewidth = 1.5)

class Wing(Component):
    def __init__(this, x: float, mat: str, param: list, angles: list, afile: str, type = 'wing'):
        """
        this is assuming linear taper, single section wing only, constant airfoil profile throughout.
        angles in the order [sweep, dihedral, incidence, twist]
        if want c/4 no sweep, put sweep = -1
        param in the order in the order [AR, b, S, taper_ratio, c_r, c_t] 
        and for the ones that aren't defined, set to -1
        """
        super().__init__(x, mat, type = 'aero', aero = 'wing')
        this.planform = {}; this.angles = {}
        this.planform['AR'], this.planform['b'], this.planform['S'], this.planform['tr'], this.planform['cr'], this.planform['ct'] = planform(param)

        #print(this.planform)
        if angles[0] == -1: this.angles['sweep'] = np.arctan(0.125*this.planform['b']*this.planform['cr']*(1-this.planform['tr'])) #for quarter no sweep
        else: this.angles['sweep'] = np.radians(angles[0]); #leading edge sweep
        this.angles['dihedral'] = np.radians(angles[1]); this.angles['inc'] = angles[2]; this.angles['tw'] = angles[3]
        #to find mass
        x = np.linspace(0, this.planform['b']/2, 50) #your discretisation wooo
        xsec = csec(afile) * ((this.planform['ct'] - this.planform['cr'])/ (this.planform['b']/2) * x + this.planform['cr'])
        #print(xsec)
        this.mass = this.material * 2 * np.trapezoid(xsec, x)   
        xcentroid = centroid(afile) * ((this.planform['ct'] - this.planform['cr'])/ (this.planform['b']/2) * x + this.planform['cr'])
        this.x_cg = np.trapezoid(xcentroid, x)/ (0.5 * this.planform['b']) #this is like average value theorem

        this.thickness = thickness(afile) #technically this is the thickness-to-chord, since dat file normalised it with c = 1
    
    def geometry(this):
        x = [this.x,
             this.x + this.planform['b']/2 *np.tan(np.pi/2 - this.angles['sweep']),
             this.x + this.planform['b']/2 * np.tan(np.pi/2 - this.angles['sweep']) + this.planform['ct'],
             this.x + this.planform['cr'],
             this.x + this.planform['b']/2 *np.tan(np.pi/2 - this.angles['sweep']) + this.planform['ct'],
             this.x + this.planform['b']/2 * np.tan(np.pi/2 - this.angles['sweep'])]
        y = [0,
             this.planform['b']/2, this.planform['b']/2,
             0,
             -this.planform['b']/2, -this.planform['b']/2]
        return plt.Polygon(xy = list(zip(x,y)), color='tab:purple', fill = False, linewidth = 1.5)

class LandingGear(Component):
    def geometry(this):
        return plt.Rectangle(xy = (this.x, -0.5), width = 1, height = 1, color = 'tab:green', fill = True, linewidth = 1.5)
    
class Fuselage(Component):
    def __init__(this, x: float, section: str, length: float, xsec1: list, xsec2 = None):
        this.cfg = {} #so the code doesn't break
        """define the length of the section (in)
        define the frontal cross section in [base, height]
        define the secondary cross section in [base, height]
        if rectangular, no need to define xsec2"""
        this.x = x

        match section:
            case "main": this.mass, this.x_cg = main_fuse(xsec1[0], xsec1[1], length)
            case "nose": this.mass, this.x_cg = nose_fuse(xsec1[0], xsec1[1], xsec2[0], xsec2[1], length)
            case "aft": this.mass, this.x_cg = aft_fuse(xsec1[0], xsec1[1], xsec2[0], xsec2[1], length)
        
        #for plotting
        this.b1 = xsec1[0]
        if xsec2 is not None: this.b2 = xsec2[0]
        else: this.b2 = this.b1
        this.length = length
        
    def geometry(this):
        x = [this.x, this.x + this.length,
             this.x + this.length, this.x]
        y = [this.b1/2, this.b2/2, 
             -this.b2/2, -this.b1/2]
        return plt.Polygon(xy = list(zip(x,y)), color='tab:purple', fill = False, linewidth = 1.5)

class Plane():
    def __init__(this, name, all: list, Sref):
        this.name = name
        this.components = all
        this.Sref = Sref
        this.avl = []
        for component in all:
            if 'aero' in component.cfg:
                this.avl.append(component)

    def plane_plot(this, components: list): #plot the plane: wing, hstab, vstab, fuse, boom
        plt.figure()
        plane_plot = plt.gcf()
        ax = plane_plot.gca()
        for geometry in components:
            ax.add_patch(geometry.geometry())
        ax.axis('equal')
        # set the figure limit
        plt.xlim([-5, 90]); plt.ylim([-100, 100])
        plt.title(this.name)
    
    def to_avl(this):
        wing_loc = []; hstab_loc = []; vstab_loc = []
        pf_vstab = {}; pf_hstab = {}
        for surf in this.avl:
            xle = 0.0; yle = 0.0; zle = 0.0; 
            xle = surf.x #distance from nose
            if 'yle' in surf.cfg:
                yle = surf.cfg['yle']
            if 'zle' in surf.cfg:
                zle = surf.cfg['yle']
            match surf.cfg['aero']:
                case 'wing':
                    wing_loc = [xle, yle, zle]
                    pf = surf.planform #dictionary
                    ang = surf.angles #dictionary
                #empennage is defined by the AR, area, and taper ratio that is typically 1 (might change this)   
                case 'hstab':
                    hstab_loc = [xle, yle, zle]
                    if not 'tr' in surf.cfg: tr = 1
                    else: tr = surf.cfg['tr']
                    param = [surf.cfg['AR'], -1, surf.cfg['area'], tr, -1, -1]
                    pf_hstab['AR'], pf_hstab['b'], pf_hstab['S'], pf_hstab['tr'], pf_hstab['cr'], pf_hstab['ct'] = planform(param)
                case 'vstab':
                    vstab_loc = [xle, yle, zle]
                    if not 'tr' in surf.cfg: tr = 1
                    else: tr = surf.cfg['tr']
                    param = [surf.cfg['AR'], -1, surf.cfg['area'], tr, -1, -1]
                    pf_vstab['AR'], pf_vstab['b'], pf_vstab['S'], pf_vstab['tr'], pf_vstab['cr'], pf_vstab['ct'] = planform(param)
        aircraft_geometry(this.name, 
                          pf = pf, #[AR, b, S, taper_ratio, c_r, c_t]
                          ang = ang, #[sweep, dihedral, incidence, twist]
                          wing_loc = wing_loc, hstab_loc = hstab_loc, vstab_loc = vstab_loc,
                          pf_h = pf_hstab, pf_v = pf_vstab,
                          afs = {'main': "C:/Users/eneiche/Documents/airfoil-sauce/airfoil_library/Selig/S7055.dat",
                                 'hstab': None, 'vstab': None}
        )

    
    ### this section is for MTOW tally
    def MTOW_tally(this):
        m_tot = 0
        for component in this.components:
            m_tot = m_tot + component.mass
        #MTOW = 9.8 * m_tot/1000 #Newtons
        MTOW = 32.17405 * m_tot #slugs to lbs
        return MTOW 

    ### this section is for CG tally
    def CG_tally(this):
        m_tot = 0 
        tally = 0 
        for component in this.components:
            loc = component.x_cg + component.x
            tally = tally + (component.mass * loc)
            m_tot = m_tot + (component.mass)
        cg = tally/m_tot
        #m_tot = m_tot/1000 #kg
        m_tot = m_tot #slugs
        plt.plot(cg, 0, 'x')
        return float(cg), float(m_tot)

    ### this section is for drag buildup - wip
    def drag_tally(this, external, u):
        global rho
        cd_tot = 0
        for component in external:
            cd_tot = cd_tot + component.drag_visc(u, this.Sref)
        D = 0.5 * rho * u **2 * this.Sref * cd_tot
        return D


