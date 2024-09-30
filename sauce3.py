import numpy as np
import matplotlib.pyplot as plt

Re = 250000
AR = np.linspace(4, 8, 100)
t = 0.09
S_exp = lambda AR: (1.83 ** 2) / AR #S = b^2 / AR
S_wet = lambda t, exp: 2*(1 + 0.2*t) * exp
FF = lambda t: 1 + 2 * t + 60 * t ** 4
c_f = 0.027/(Re**(1/7))

c_d = lambda FF, wet, exp: (c_f * FF * wet)/exp

plt.figure()
plt.plot(AR, S_exp(AR))
plt.figure()
plt.plot(AR, 0.5 * 1.225 * 20 ** 2 * 0.5 * S_exp(AR) * c_d(FF(t), S_wet(t, S_exp(AR)), S_exp(AR)))
plt.xlabel("AR"); plt.ylabel("$c_d$")
plt.title("Re = 250000; b = 1.83 m, thickness = 9%, $\lambda = 1$")
plt.show()