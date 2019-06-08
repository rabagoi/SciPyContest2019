import numpy as np
from matplotlib.cm import get_cmap

# Creates an array of points between a start point p1 and end point p2.
def Interpolation(p1, p2, N=24):
    return np.array([0.5*((p2+p1)-(p2-p1)*np.cos(np.pi/N*t)) for t in range(N)])


def TrajColorMap(inclination):
    crit_inc = 38.5      # Approximate critical angle for e=0.5
    # Colormap for the coplanar region
    polar_angle = 90-np.fabs(90-inclination)
    if 90-np.fabs(90-inclination) < crit_inc:
        return get_cmap('viridis')(0.3+polar_angle/crit_inc)
    else:
        return get_cmap('plasma_r')(0.0+(polar_angle-crit_inc)/(90-crit_inc))

