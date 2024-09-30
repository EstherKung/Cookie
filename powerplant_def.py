import matplotlib.pyplot as plt
import numpy as np

T = [33.0472538,
25.3310838,
16.3362342,
15.8953102,
16.4685114,
16.203957,
19.2463326,
20.1061344,
7.99395212,
9.09626212,
11.41331774,
9.7444204,
6.21261916,
7.3965001,
8.20559564] #thrust in lbs

mass = [435,
590,
310,
317,
362,
370,
420,
418,
89,
88,
125,
126,
59,
83.4,
83.8] #mass in g

#plt.scatter(T, mass)
#a, b = np.polyfit(T, mass, 1)
#x = np.linspace(1, 40, 200)
#plt.plot(x, a*x + b)
#plt.xlabel("Thrust in lbs")
#plt.ylabel("mass of motor in g")
#plt.show()

### polynomial regression
c = -310
d = 56.3
f = -0.98

def propulsor_mass(thrust):
    #mass = a*thrust + b
    #mass = c + d*thrust + f* thrust**2
    mass = 3.76* thrust ** 1.54
    mass = mass * 6.85218e-5 #convert from g to slugs
    return mass * 1.25 #a 25% factor for the propeller
