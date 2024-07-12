import numpy as np

"""
std_atmosphere
standard imperial engineering atmosphere, data from:
https://www.engineeringtoolbox.com/standard-atmosphere-d_604.html
Call 'rho_std(altitude)'
give altitude in ft
returns rho in slug/ft**3
"""
h = np.arange(0,45001,5000)
rho_dat = np.array([23.77, 20.48, 17.56, 14.96, 12.67, 10.66, 8.91, 7.38, 5.87, 4.62])
rho_dat /= 10000
rho_std = lambda x: np.interp(x, h, rho_dat) #give x in feet, returns rho in slug/ft**3
