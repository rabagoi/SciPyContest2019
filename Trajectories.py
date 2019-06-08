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

R = 4.5
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
        # Stop the integration after the vector has "looped".
        prec_angle = np.arctan2(L_Kozai[1], L_Kozai[0])
        if prec_angle < 0 and len(inc_p) > 0 and inc_p[-1] > 0:
            # Add the first data point to the list to close the "loop".
            # This does not always happen with calculated data, but it produces a better plot.
            trajectory.append(trajectory[0])
            break
        inc_p.append(prec_angle)
    return np.array(trajectory).T

# Project the 3D trajectories into 2D using a stereographic projection,
# plus a mapping to a disk of radius R to keep reasonable sizing.
def StereographicProjection(xyz, normalized=False):
    proj = [xyz[0]/(R-xyz[2]), xyz[1]/(R-xyz[2]), 0*xyz[2]]
    if normalized:
        mag = np.sqrt(proj[0]**2+proj[1]**2)
        norm = 2*R*2/np.pi*np.arctan(mag)
        proj = norm*(proj/mag)
    return proj

# Generate 2D and 3D interpolated trajectories
def CreateTrajectories(N=180):
    traj_animated = []
    traj_color = []
    # Create trajectories
    for i in range(18):
        color = TrajColorMap(10*i)
        #print(color)
        traj3d = GetLTrack(10*i+5)
        traj2d = StereographicProjection(traj3d, normalized=True)
        traj_animated.append(Interpolation(traj3d, traj2d, N))
        traj_color.append(color)

    return traj_animated, traj_color

# For standalone plots of Scenes 3 and 4
'''
fig = plt.figure(figsize=(10,10), facecolor='gray')
ax = plt.subplot(111, projection='3d')


traj_animated, traj_color = CreateTrajectories()

print(len(traj_animated), len(traj_animated[0]), len(traj_animated[0][0]) )
def anim(t):
    i_t = t
    # Transform from 3D to 2D
    plt.cla()
    #plt.axis('off')
    ax.set_xlabel("X")
    ax.set_facecolor((0.7, 0.7,0.7))
    ax.set_xlim(-R,R)
    ax.set_ylim(-R,R)
    ax.set_zlim(-R,R)
    for i in range(18):
        if t <=180:
            # Highlighting individual traces
            #if 10*i < i_t:
            traj_alpha = 0.9 if 10*i < i_t else 0.0
            color = traj_color[i] if 10*i+3 < i_t else (1,1,1)
            ax.plot(traj_animated[i][0][0], traj_animated[i][0][1], traj_animated[i][0][2], color = color, alpha=traj_alpha)

        else:
            # Transform from 3D to 2D
            traj_frame = traj_animated[i][int(t-180)]
            ax.plot(traj_frame[0], traj_frame[1], traj_frame[2], color = traj_color[i])

a = FuncAnimation(fig, anim, frames=np.append(Interpolation(0, 180, 180), Interpolation(180, 360, 180)), interval=20)
plt.show()
'''