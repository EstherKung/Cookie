import numpy as np
import matplotlib.pyplot as plt
from implicit import secant_method
from std_atmosphere import *
# from labellines import labelLines

"Calculate reasonable limiTW"
"---------------------------"
"TW - Thrust Loading"
"WS - Wing Loading"

g = 32.17
Tn, Tx = 6.494, 9.363 # lbf, thrust @ v=36.1 ft/s to static thrust
bn, bx = 4, 6 # ft, min span to max span
cn, cx = 0.75, 2 # ft, for AR 5-8
Sn, Sx = bn*cn, bx*cx # ft**2
Wn, Wx = 12, 40 # lbf structural to max

WSn, WSx = Wn/Sx, Wx/Sn
TWn, TWx = Tn/Wx, Tx/Wn
# print(WSn, WSx, TWn, TWx)

WS = np.linspace(WSn, WSx, 31)

class Constraint_Diagram:
    def __init__(self):
        self.rho = 23.77e-4
        self.g = 32.17
        self.fig = plt.figure('Constraint Diagram')
        self.ax = plt.axes()
        # self.ax = self.fig.add_subplot(1, 1, 1)
        self.TWn = 0.2; self.TWx = 0.2
        self.cfg = {}
        self.cfg['rho'] = 23.77e-4
        self.cfg['g'] = 32.17

    def Constraint_Stall(self, V_stall: float = None):
        """
        \nPlots the WS constraint for the desired stall velocity V_stall
        \nINPUT: V_stall
        \nOUTPUT: None
        \n──────────────────────────────────────
        \nIndicates how slowly you can fly before stalling; lower V_stall is more constraining
        \nThis constraint assumes you are trying to fly slowly while achieving a load factor of n=1
        \nBut don't get it twisted; stall occurs at a certain angle of attack due to flow separation
        \nLower WS DOES NOT delay stall (i.e, stall at higher angle of attack). High Reynolds numbers, 
        \nslats, vortex generators, and turbulent boundary layers delay stall.
        \nLower WS DOES let you fly slower near the stall angle of attack. This is desirable for
        \nnearly all aircraft that want to fly "slow" (~150 knots = 172 mph = 277.8 km/h) at landing.
        \nFor AMAT, we would like to fly slow at landing, but at competition most planes approach
        \nthe landing strip at speeds only slightly below cruise speed.
        """
        if not V_stall: V_stall = self.cfg['Vst']
        WS_stall = 0.5*self.cfg['rho']*V_stall**2*self.cfg['CLx']
        self.ax.plot([WS_stall, WS_stall], [0, 0.5], label='Stall at {:3.1f} ft/s'.format(V_stall))
 
    def Constraint_Cruise(self, V_cruise: float = None, k1: float = 1.0):
        """
        \nPlots the cruise constraint for a desired cruise velocity.
        \nINPUT: V_cruise
        \nOUTPUT: None
        \n────────────────────────────────────────
        \nAssumes small angles. Some notes on small angle approximation:
        \n\t99% accurate up to 14 degrees (Good starting guess for stall alpha: 14 degrees)
        \n\t99.5% accurate up to 9.9 degrees
        \n\t99.8% accurate up to 6.2 degrees (Common cruise alpha: 3-5 degrees)
        \nRemember: we discovered that flying at the best L/D results in a constant
        \nT/W cruise constraint. The constant T/W is equal to the generic cruise constraint's minimum.
        """
        if not V_cruise: V_cruise = self.cfg['Vcr'] #Default case
        q = 0.5*self.rho*V_cruise**2 #Dynamic pressure
        TW = (q*self.cfg['CD0']/self.cfg['WS'] + self.cfg['k']/q*self.cfg['WS']) #Constraint equation
        TW *= k1 #TW correction
        self.ax.plot(self.cfg['WS'], TW, label='Cruise at {:3.1f} ft/s'.format(V_cruise))
    
    def Constraint_Maneuver(self, V_maneuver:float = None, phi:float = None, k1: float = 1.0):
        """
        \nPlots the maneuver constraint for given velocity and bank angle
        \nINPUT: V_maneuver, phi
        \nOUTPUT: None
        \n────────────────────────────────────────
        """
        if not V_maneuver: V_maneuver = self.cfg['Vcr'] #Default case
        if not phi: phi = self.cfg['phi']
        n = 1/np.cos(phi*np.pi/180) #Load factor from φ
        q = 0.5*self.rho*V_maneuver**2 #Dynamic pressure
        TW = (q*self.cfg['CD0']/self.cfg['WS'] + self.cfg['k']*n**2/q*self.cfg['WS']) #Constraint equation
        TW *= k1 #TW correction factor
        self.ax.plot(self.cfg['WS'], TW, label='Sustained {:2.0f}$^\circ$ cruise turn'.format(phi), linestyle='--')

    def Constraint_Climb(self, G: float = None, k1: float = 1.0):
        """
        \nPlots the climb constraint for given climb gradient
        \nINPUT: G 
        \nOUTPUT: None
        \n────────────────────────────────────────
        """
        if not G: G = self.cfg['G']
        TW = self.cfg['ks']**2*self.cfg['CD0']/ self.cfg['CLx'] + self.cfg['k']*self.cfg['CLx']/self.cfg['ks']**2 + G #Constraint equation
        TW *= k1 #Correction factor
        self.ax.plot(self.cfg['WS'], np.ones_like(self.cfg['WS'])*TW, label='Climb at {:3.1f} ft/s, G={:4.3f}'.format(self.cfg['ks']*27.1, G))

    def Constraint_Takeoff(self, takeoff_distance:float = None):
        """
        \nPlots the takeoff constraint
        \nINPUT: takeoff_distance
        \nOUTPUT: None
        \n────────────────────────────────────────
        \nUses a linear model of the thrust versus velocity
        """
        if not takeoff_distance: takeoff_distance = self.cfg['Lto']
        def takeoff_error(ws, TW):
            # print('to error call', ws, TW)
            acc = lambda V: self.g*(TW - 0.5*self.rho*V**2/ws*(self.cfg['CDto'] - 0.03*self.cfg['CLto']) - 0.03)
            Vto = np.sqrt((2*ws)/(self.rho*0.8*self.cfg['CLx'])) 
            v = np.zeros(2)
            x = np.zeros(2)
            dt = 0.005
            while x[1] < .9*takeoff_distance:
                v[0] = v[1]
                x[0] = x[1]
                v[1] = v[0] + dt*acc(v[0])
                x[1] = x[0] + dt*v[1]
            return Vto - v[0]
        #TW1 = 0.1; TW2 = 0.11 # Don't change this unless you have to
        TW1 = 0.1; TW2 = 0.07

        self.ax.plot(self.cfg['WS'], secant_method(self.cfg['WS'], TW1, TW2, takeoff_error), label = 'Takeoff')

    def Constraint_Landing(self, landing_distance: float = None, V_landing: float = None):
        """
        \nPlots the landing constraint
        \nINPUT: landing_distance, V_landing
        \nOUTPUT: None
        \n────────────────────────────────────────
        \nAs mentioned before, most RC planes will land at just below their cruise speed
        """
        if not landing_distance: landing_distance = self.cfg['Lld']
        if not V_landing: V_landing = self.cfg['Vld']
        dt = 0.05
        #WS1 = 0.9
        WS1 = 0.4
        WS2 = 0.6
        #WS2 = 1.1
        def landing_error(Vld, ws):
            x = 0
            while Vld > 0:
                Vld -= dt*self.g*( 0.03 + 0.5*self.cfg['rho']*Vld**2/ws*( self.cfg['CDto'] - 0.03*self.cfg['CLto'] ) )
                x += Vld*dt
            return 0.9*landing_distance - x
        WS = secant_method(np.array([V_landing]), WS1, WS2, landing_error)
        self.ax.plot([WS, WS], [0, 1.0], label='Landing at {:3.1f} ft/s'.format(V_landing))


cd = Constraint_Diagram()

"Altitude"
cd.cfg['rho'] = rho_std(0)
cd.cfg['g']   = 32.17 #ft/s**2

"Constraint defaults"
cd.cfg['Lto'] = 300 #Takeoff distance, ft
cd.cfg['Lld'] = 600 #Landing distance, ft
cd.cfg['G']   = 0.1 #Climb gradient
cd.cfg['ks']  = 1.25 #Climb configuration
cd.cfg['phi'] = 45  #Sustained bank angle, deg



"Speeds"
cd.cfg['Vcr'] = 80.0 #Cruise, ft/s
cd.cfg['Vld'] = 35.0 #Landing/approach, ft/s
cd.cfg['Vst'] = 23.0 #Stall, f/ts

"Aerodynamics"
cd.cfg['k']    = 1/(np.pi*0.9*6.3)
cd.cfg['CD0']  = 0.02 #Zero-lift drag
cd.cfg['CLx']  = 1.4  #Maximum 
cd.cfg['CLto'] = 0.15 #Takeoff configuration
cd.cfg['CDto'] = 0.02 #Takeoff configuration


"WS grid"
cd.cfg['WS']  = np.linspace(0.1,2.2,21)

# cd.Constraint_Stall()
cd.Constraint_Cruise()
cd.Constraint_Climb()
cd.Constraint_Landing()
cd.Constraint_Takeoff()
cd.Constraint_Maneuver(phi=45)
cd.Constraint_Cruise(V_cruise=31.0)

plt.legend(loc='upper right', fontsize=14)
plt.title('Constraint diagram', fontsize=18)
plt.xlabel(r'Wing loading $(lb/ft^2)$', fontsize=14)
plt.ylabel(r'Thrust-to-weight ratio', fontsize=14)
plt.xlim([0.1, 2.2])
plt.ylim([0.0, 0.5])
plt.grid()
plt.show()