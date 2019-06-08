import rebound
from rebound.plotting import fading_line
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Line3DCollection
from matplotlib.animation import FuncAnimation
from AnimationTools import TrajColorMap, Interpolation
from Trajectories import CreateTrajectories

matplotlib.rcParams['mathtext.fontset'] = 'cm'

'''Matplotlib and REBOUND code to plot binary orbits.'''
def Plot_Binary(sim, t, is3d=False):    
    #Wrapper function for adding a line to the line collection.
    # Changes to the 3D version for 3D plots.
    adeline = ax.add_collection3d if is3d else ax.add_collection
    
    # Alpha setup for Scenes 3 and 4
    orbit_alpha = 1.0 if t < 400 else max(0.0, 1.0-0.01*(t-400))
    for i in range(len(sim.particles)):
        p = sim.particles[i]
        if i < 2:
            #Plot stars
            ax.scatter3D(p.x, p.y, p.z, marker='o', s=100, color='w', edgecolor='gold', alpha=orbit_alpha)
            
            #Create and plot the binary orbits using the data from their first orbit
            if sim.t < 1.01*P_b:
                binary_traj[i].append([p.x, p.y, p.z])
            bt = np.array(binary_traj[i])
            ax.plot3D(bt[:,0], bt[:,1], bt[:,2], color='gray', alpha=orbit_alpha)
            #adeline(binary_orbits[i], zs = o[:-1,2])
        else:
            # Plot circumbinary particle
            p_color = TrajColorMap(p.inc*180/np.pi)
            ax.scatter3D(p.x, p.y, p.z, color=p_color, edgecolor='k', alpha=orbit_alpha, zorder=200)
            o = np.array(p.sample_orbit())
            lc = fading_line(o[:,0], o[:,1])

            # Convert to Line3DCollection
            if is3d:
                pts = o.reshape(-1,1,3)
                #print(pts)
                seg = np.concatenate((pts[:-1], pts[1:]), axis=1)
                #print(seg)
                lc = Line3DCollection(seg, color=p_color, alpha=orbit_alpha, zorder=1)

            adeline(lc)
            
            L_Kozai = np.cross([p.x, p.y, p.z], [p.vx, p.vy, p.vz])
            
            R = 4.5
            #Graph a unit angular momentum vector
            L_Kozai = R*L_Kozai/np.sqrt(np.sum(L_Kozai**2))    
            ax.quiver(0,0,0,L_Kozai[0], L_Kozai[1],L_Kozai[2], color=p_color, alpha=orbit_alpha)

            # plot the overall trajectory of L
            trajectory[i-2].append(L_Kozai)
            
            traj = np.array(trajectory[i-2]).T
            ax.scatter(traj[0], traj[1], traj[2], color=p_color, s=10, alpha=orbit_alpha, zorder=3)


    # Plot a transparent sphere
    N = 100
    phi = np.linspace(0, 2*np.pi, N)
    th = np.linspace(0, np.pi, N)
    spx = R*np.outer(np.cos(phi), np.sin(th))
    spy = R*np.outer(np.sin(phi), np.sin(th))
    spz = R*np.outer(np.ones(len(phi)), np.cos(th))
    sp_alpha =  0.1 if t < 600 else max(0.0, 0.1-0.001*(t-600))
    ax.plot_surface(spx,spy,spz, color='gray', alpha=sp_alpha, zorder=2)
    ax.plot(R*np.cos(phi), R*np.sin(phi), 'k--', alpha=sp_alpha)
    
    # Plot text
    text_alpha = 1.0 if t < 25 else max(0.0, 1.0-0.05*(t-25))
    ax.text(7, 0, 1, r"$\mathit{i = 20\degree}$", alpha=text_alpha, fontsize=14)
    ax.text(1, 0, 7, r"$\mathit{i = 70\degree}$", alpha=text_alpha, fontsize=14)

    
''' Animation function '''
def anim(t):
    plt.cla()
    Plot_Binary(sim, t, is3d=is3d)
    # Scene 1: Evolution of the Circumbinary Orbits (t=0 to t=300)
    if t < 300:
        ax.view_init(azim=280,elev=10)
        sim.integrate(13*t)
    # Scene 2: Pan Camera and fade out orbits (t=300 to t=400)
    elif t > 300 and t < 390:
        az = 280+(t-300)/3
        alt = (t-300)/5 + 10
        ax.view_init(azim=az, elev=alt)

    if t > 400:
        i_t = t-400

        # Plot text for the sweeping inclination
        inc_alpha = 1.0 if t < 580 else max(0.0, 1.0-0.01*(t-580))
        ax.text(0, 8, 8, r"$\mathit{{i = {:.0f}\degree}}$".format(min(i_t, 180)), fontsize=16, alpha=inc_alpha)
        for i in range(18):
            # Scene 3: Add traces for all inclinations
            if t < 580:
                traj_alpha = 1.0 if 10*i < i_t else 0.0                 # Hide the trajectory until the proper inclination is reached
                color = traj_color[i] if 10*i+3 < i_t else (1,1,1)      # Flash white momentarily to highlight the trajectory
                ax.plot(traj_animated[i][0][0], traj_animated[i][0][1], traj_animated[i][0][2], color = color, alpha=traj_alpha)

            # Scene 4: Project sphere onto 2D surface and pan camera
            elif t >= 600 and t < 840:
                traj_frame = traj_animated[i][int(t-600)]
                ax.plot(traj_frame[0], traj_frame[1], traj_frame[2], color = traj_color[i])
                az = 310+(t-600)*5/24
                alt = 28+(t-600)/10
                ax.view_init(azim=az, elev=alt)
            # Scene 5: Still shot of the final projection
            elif t >= 840:
                traj_frame = traj_animated[i][-1]
                ax.plot(traj_frame[0], traj_frame[1], traj_frame[2], color = traj_color[i])


    # Set plot limits and plot color
    ax.set_facecolor((0.8,0.8,0.8))
    ax.set_xlim(-5, 5)
    ax.set_ylim(-5, 5)
    ax.set_zlim(-5, 5)
    plt.axis('off')
    


'''---------------------------Main Variables and Setup-------------------------------------'''

# Variables for binary initial conditions
i_o = 10*np.pi/180
omega_o = np.pi/2 if i_o < np.pi/2 else -1*np.pi/2
is3d = True
subplot_type = {'projection':'3d'} if is3d else {}
fig, ax = plt.subplots(figsize=(8,8), facecolor=(0.8,0.8,0.8), subplot_kw = subplot_type)   
trajectory = [[],[]]


#Setup binary system
sim = rebound.Simulation()
sim.add(m=1)
sim.add(m=1, a=2, e=0.5)    #Separation = 1
#Add test particle at 60 degree inclination, Omega=90 degrees
sim.add(m=0, a=6, inc=i_o, Omega=omega_o)
sim.add(m=0, a=7, inc=np.pi/2-i_o, Omega=omega_o)
sim.move_to_com()

P_b = sim.particles[1].P   #Binary orbital period
binary_traj = [[],[]]

# Create plots for Scenes 3 and 4
traj_animated, traj_color = CreateTrajectories(240)

# Times for each scene
time_s1 = Interpolation(0,300,1000)
time_s2 = Interpolation(300,390,240)
time_s3 = Interpolation(400,580,360)
time_s4 = Interpolation(600,840,240)
time_s5 = np.linspace(900, 1000, 100)
times = np.concatenate( (time_s1,time_s2,time_s3, time_s4, time_s5) )


# Create and save animation
a = FuncAnimation(fig, anim, frames=times, interval=20)
a.save("CB_Kozai.avi", dpi=300, savefig_kwargs={'facecolor':(0.8,0.8,0.8)})