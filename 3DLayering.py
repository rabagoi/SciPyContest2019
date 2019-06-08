import rebound
from rebound.plotting import fading_line
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import line_collection_2d_to_3d, Line3DCollection
from matplotlib.collections import LineCollection
from matplotlib.animation import FuncAnimation
from IPython.display import HTML
from os import system
from AnimationTools import Interpolation

# LinearSegmentedColormaps needed for fading line


#Matplotlib and REBOUND code to plot binary orbits.
def Plot_Binary(sim, t, is3d=False):    
    #Wrapper thing for adding a line to the line collection.
    # Changes to the 3D version for 3D plots.
    adeline = ax.add_collection3d if is3d else ax.add_collection
    
    for i in range(len(sim.particles)):
        p = sim.particles[i]
        if i < 2:
            #Plot stars
            ax.scatter3D(p.x, p.y, p.z, marker='*', s=100, color='k')
            
            #Creates orbits for binary systems
            o = np.array(p.sample_orbit(primary=sim.particles[1])) if i==0 else np.array(p.sample_orbit())   
            lc = fading_line(o[:,0], o[:,1])
            
            if sim.t < 1.1*P_b:
                binary_traj[i].append([p.x, p.y, p.z])
            bt = np.array(binary_traj[i])
            ax.plot3D(bt[:,0], bt[:,1], bt[:,2], color='gray')
            #adeline(binary_orbits[i], zs = o[:-1,2])
        else:
            # Plot circumbinary particle
            ax.scatter3D(p.x, p.y, p.z, color='k', zorder=200)
            o = np.array(p.sample_orbit())
            lc = fading_line(o[:,0], o[:,1], color='blue')
            #print("COLOR", lc.get_color())
            # Convert to Line3DCollection
            if is3d:
                pts = o.reshape(-1,1,3)
                #print(pts)
                seg = np.concatenate((pts[:-1], pts[1:]), axis=1)
                #print(seg)
                lc = Line3DCollection(seg, cmap=lc.get_cmap(), zorder=1)

            adeline(lc)
            #Calculate angular momentum vector
            #L_Kozai = np.array([p.y*p.vz - p.z*p.vy,
            #     p.z*p.vx - p.x-p.vz,
            #     p.x*p.vy - p.y*p.vx])
            
            L_Kozai = np.cross([p.x, p.y, p.z], [p.vx, p.vy, p.vz])
            #print(p.v, p.vx, p.vy, p.vz)
            
            R = 3
            #Graph a unit angular momentum vector
            L_Kozai = R*L_Kozai/np.sqrt(np.sum(L_Kozai**2))    
            #Lx = [0, L_Kozai[0]]
            #Ly = [0, L_Kozai[1]]
            #Lz = [0, L_Kozai[2]]

            #ax.plot3D(Lx, Ly, Lz, 'r-o')
            ax.quiver(0,0,0,L_Kozai[0], L_Kozai[1],L_Kozai[2], color=arrow_color[i-2])

            # plot the overall trajectory of L
            trajectory[i-2].append(L_Kozai)
            # Calculate and check precession angle
            # Can be done with angle range [-pi, pi]?
            prec_angle = np.arctan2(L_Kozai[1], L_Kozai[0])
            if prec_angle < 0 and len(inc_p) > 0 and inc_p[-1] > 0:
                print("Loop at t =", t)
                
            inc_p.append(prec_angle)
            
            traj = np.array(trajectory[i-2]).T
            ax.scatter(traj[0], traj[1], traj[2], color=traj_color[i-2], s=10, alpha=0.5, zorder=3)


    # Plot a transparent sphere
    
    N = 100
    phi = np.linspace(0, 2*np.pi, N)
    th = np.linspace(0, np.pi, N)
    spx = R*np.outer(np.cos(phi), np.sin(th))
    spy = R*np.outer(np.sin(phi), np.sin(th))
    spz = R*np.outer(np.ones(len(phi)), np.cos(th))
    
    ax.plot_surface(spx,spy,spz, color='gray', alpha=0.25, zorder=2)
    
        
    # Set plot limits
    plt.xlim(-5, 5)
    plt.ylim(-5, 5)
    plt.xlabel("X")
    plt.ylabel("Y")
    ax.set_zlim(-5, 5)
    #plt.show()
    
    #return fig
    
# Animation function
def anim(t):
    if t < 260:
        plt.cla()
        plt.axis('off')
        sim.integrate(13*t)
        Plot_Binary(sim, t, is3d=is3d)
    
    az = 280 if t<300 else t-20
    alt = 10 if t<300 else (t-300)/4 + 10
    ax.view_init(azim=az, elev=alt)
    # Add a trace for the angular momentum vector???
    
'''----------------------------------------------------------------'''

i_o = 30*np.pi/180
omega_o = np.pi/2 if i_o < np.pi/2 else -1*np.pi/2
is3d = True
subplot_type = {'projection':'3d'} if is3d else {}
fig, ax = plt.subplots(figsize=(8,8), facecolor=(0.8,0.8,0.8), subplot_kw = subplot_type)   
trajectory = [[],[]]
arrow_color = ['xkcd:dark green', 'xkcd:maroon']
traj_color = ['g', 'r']

# Outside plot test
inc_p = []

#Setup binary system
sim = rebound.Simulation()
sim.add(m=1)
sim.add(m=1, a=2, e=0.5)    #Separation = 1
#Add test particle at 60 degree inclination, Omega=90 degrees
sim.add(m=0, a=6, inc=i_o, Omega=omega_o)
sim.add(m=0, a=7, inc=2*i_o, Omega=omega_o)
sim.move_to_com()

# Record binary orbits (Incorrect orbits)
binary_orbits = []
binary_traj = [[],[]]
for i in range(2):
    p = sim.particles[i]
    #Creates orbits for binary systems
    o = np.array(p.sample_orbit(primary=sim.particles[1])) if i==0 else np.array(p.sample_orbit())   
    binary_orbits.append(fading_line(o[:,0], o[:,1]))
    binary_orbits[i].set_alpha(1)

#Plot_Binary(sim, is3d=True)
P_b = sim.particles[1].P   #Binary orbital period
print(sim.particles[-1].h)
times = np.append( Interpolation(0, 260, 1000), Interpolation(300, 390, 240) )

a = FuncAnimation(fig, anim, frames=times, interval=20)
a.save("CB_Kozai.avi", dpi=300, savefig_kwargs={'facecolor':(0.8,0.8,0.8)})

print(len(trajectory), len(inc_p))
print("Average", np.average(inc_p))
print(binary_traj)

system("vlc -I dummy ~/Music/Jingle1.wav vlc://quit")