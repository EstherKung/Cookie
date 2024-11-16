
import numpy as np
import matplotlib.pyplot as plt
rho = 1.3756E-6 #slugs/in ^3 
mu = (3.737E-7)/12 #slug/(ft s)* 1ft/12in = slug / (ft in)

u = 12 * np.array([55, 60, 65, 70, 75, 80])
print (u)
t = 0.10
S_ref = np.array([1296, 1152, 1036.8, 942.54545, 864.0, 797.53846, 740.57143])
MAC = np.array([19.10226, 16.97979, 15.28181, 13.89255, 12.73484, 11.75524, 10.91558])
S_wet = 2*S_ref*((1+0.25*t) * (1 + 0.4)/(1+0.4))
print(S_wet)
FF = 1 + 2 * t + 60 * t ** 4
print(FF)

AR = [4, 4.5, 5, 5.5, 6, 6.5, 7]
fuse_wet = lambda d, l: 2 * np.pi * d/2 * l + 2*np.pi * (d/2)**2
pod_wet = [fuse_wet(5.5, 18.00),
 fuse_wet(5.5, 18.00),
 fuse_wet(4, 20),
 fuse_wet(4, 20),
 fuse_wet(3.5, 20.00),
 fuse_wet(3, 22.00),
 fuse_wet(3, 22.00)]
#boom diameter 0.75
boom_wet = [fuse_wet(0.75, 54.14), fuse_wet(0.75, 63.1249), fuse_wet(0.75, 50.40),
            fuse_wet(0.75, 49.42), fuse_wet(0.75, 48.83), fuse_wet(0.75, 48.33), fuse_wet(0.75, 43.62)]

c_vtail = np.array([7.7668, 7.39008, 7.1444, 6.9468, 6.6511, 6.3901, 6.4908])
c_htail = np.array([12.7114, 10.31641, 9.3542, 8.2607, 7.2509, 6.4359, 6.0703])
s_htail = np.array([25.4228, 23.21, 23.3855, 22.7388, 21.7708, 20.9167, 21.246])
s_vtail = np.array([11.6502, 11.08513, 10.7166, 10.4202, 9.9766, 9.5852, 9.97362])

S_h = c_vtail * s_vtail
S_v = c_htail * s_vtail

c_pod = np.array([18, 18, 20, 20, 20, 22, 22])
c_boom = np.array([54.14, 63.12, 50.40, 49.42, 48.83, 48.33, 43.62])
FF_fuse = np.array([18/5.5, 18/5.5, 20/4, 20/4, 20/3.5, 22/3, 22/3])
FF_boom = np.array([54.14/0.75, 63.1249/0.75, 50.40/0.75, 49.42/0.75, 48.83/0.75, 48.33/0.75, 43.62/0.75])

plt.figure()
for u in u: 
    #wing
    Re_wing = rho * u * MAC/mu  
    c_f = 0.027/(Re_wing**(1/7))
    c_d_wing = c_f * FF * S_wet/S_ref
    #print(f"parasite wing drag at {u}:", c_d_wing)

    #empennage, assumes flat plate
    Re_vtail = rho * u * c_vtail/mu
    Re_htail = rho * u * c_htail/mu

    c_fv = 0.027/(Re_vtail**(1/7))
    c_fh = 0.027/(Re_htail**(1/7))
    c_dv = c_fv * S_v/S_ref
    c_dh = c_fh * S_h/S_ref
    #print(f"parasite v plate at {u}:", c_fv)
    #print(f"parasite v plate at {u}:", c_fh)

    Re_pod = rho * u * c_pod/mu
    Re_boom = rho * u * c_boom/mu
    cf_pod = 0.027/(Re_pod**(1/7))
    cf_boom = 0.027/(Re_boom**(1/7))
    
    cd_pod = cf_pod * FF_fuse * pod_wet/S_ref
    cd_boom = cf_boom * FF_boom * boom_wet/S_ref

    cd_total = cd_pod + cd_boom + c_dv + c_dh + c_d_wing
    print(f"at speed = {u/12} ft/s, total drag coefficient: {cd_total}")
    D_total = 1/2 * 0.0023769 * ((u/12) ** 2) * S_ref/144 * cd_total
    print(f"at speed = {u/12} ft/s, total drag of plane is {D_total} lbs.")

    plt.plot(AR, D_total)





plt.legend(["55", "60", "65", "70", "75", "80"])

plt.show()























