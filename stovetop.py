from plane_def import *
afile = "airfoil_library/S4233.dat"

wing = Wing(-20, "foam", [6, -1, 130, 0.4, -1, -1], [-1, 0, 0, 0], afile)

plane = [wing]

test = Plane("haha", plane, 0.66)
test.plane_plot(plane)
plt.show()

