import rebound
from rebound.plotting import fading_line
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import line_collection_2d_to_3d, Line3DCollection
from matplotlib.collections import LineCollection
from matplotlib.cm import get_cmap
from AnimationTools import TrajColorMap, Interpolation
from matplotlib.animation import FuncAnimation

R = 1
N = 100
def GetLTrack(inclination):
    trajectory = []
    inc_p = []
    i_o = inclination*np.pi/180
    omega_o = np.pi/2 if i_o < np.pi/2 else -1*np.pi/2
    #Setup binary system
    sim = rebound.Simulation()
    sim.add(m=1)
    sim.add(m=1, a=2, e=0.5)    #Separation = 1
    #Add test particle at 60 degree inclination, Omega=90 degrees
    sim.add(m=0, a=6, inc=i_o, Omega=omega_o)
    sim.move_to_com()
    
    while (sim.t < 10000):
        sim.integrate(sim.t+30)
        p = sim.particles[-1]
        L_Kozai = np.cross([p.x, p.y, p.z], [p.vx, p.vy, p.vz])
        #Graph a unit angular momentum vector
        L_Kozai = R*L_Kozai/np.sqrt(np.sum(L_Kozai**2))    
        # plot the overall trajectory of L
        trajectory.append(L_Kozai)
        # Calculate and check precession angle
        # Can be done with angle range [-pi, pi]?
        prec_angle = np.arctan2(L_Kozai[1], L_Kozai[0])
        if prec_angle < 0 and len(inc_p) > 0 and inc_p[-1] > 0:
            print("Loop at t =", sim.t)
            trajectory.append(trajectory[0])
            break
        inc_p.append(prec_angle)
    return np.array(trajectory).T

def StereographicProjection(xyz, normalized=False):
    proj = [xyz[0]/(R-xyz[2]), xyz[1]/(R-xyz[2]), 0*xyz[2]]
    if normalized:
        mag = np.sqrt(proj[0]**2+proj[1]**2)
        norm = 2*R*2/np.pi*np.arctan(mag)
        proj = norm*(proj/mag)
    return proj

fig = plt.figure(figsize=(10,10), facecolor='gray')
cmap_f = get_cmap('viridis')
cmap_p = get_cmap('plasma')
ax = plt.subplot(111, projection='3d')
ax.set_facecolor('w')
ax.set_xlim(-R,R)
ax.set_ylim(-R,R)
ax.set_zlim(-R,R)
plt.axis('off')
#ax.view_init(azim=90, elev=5)

# Plot sphere
phi = np.linspace(0, 2*np.pi, N)
th = np.linspace(0, np.pi, N)
x = 0.95*R*np.outer(np.cos(phi), np.sin(th))
y = 0.95*R*np.outer(np.sin(phi), np.sin(th))
z = 0.95*R*np.outer(np.ones(len(phi)), np.cos(th))
#ax.plot_surface(x,y,z, color='gray', alpha=0.9, zorder=-1)

'''
# Create trajectories
for i in range(18):
    color = TrajColorMap(10*i)
    #print(color)
    traj3d = GetLTrack(10*i+5)
    ax.plot(traj3d[0], traj3d[1], traj3d[2], '.-', c=color, alpha=1.0, zorder=1)
    traj2d = StereographicProjection(traj3d, normalized=True)
    ax.plot(traj2d[0], traj2d[1], traj2d[2], color=color, zorder=2)
plt.show()
'''


traj_animated = []
traj_color = []
# Create trajectories
for i in range(18):
    color = TrajColorMap(10*i)
    #print(color)
    traj3d = GetLTrack(10*i+5)
    traj2d = StereographicProjection(traj3d, normalized=True)
    traj_animated.append(Interpolation(traj3d, traj2d, 100))
    traj_color.append(color)

print(len(traj_animated), len(traj_animated[0]), len(traj_animated[0][0]) )
def anim(t):
    plt.cla()
    ax.set_xlim(-R,R)
    ax.set_ylim(-R,R)
    ax.set_zlim(-R,R)
    for i in range(18):
        traj_frame = traj_animated[i][t]

        #ax.plot(traj_frame[0], traj_frame[1], traj_frame[2], color = traj_color[i])
        ax.plot(traj_frame[0], traj_frame[1], traj_frame[2], color = traj_color[i])

a = FuncAnimation(fig, anim, frames=100, interval=100)
plt.show()


print("TRANSITION COLORS")
print(get_cmap('viridis')(1))
print(get_cmap('plasma')(0))
#Use line/scatter collections and colormaps???