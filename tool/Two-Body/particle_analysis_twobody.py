import numpy as np
import matplotlib.pylab as plt
from matplotlib.colors import LogNorm

N = 1
center_radius = 3.0
t = []
p1_vx = []
p1_vy = []
p2_vx = []
p2_vy = []
for i in range(N):
    with open("Particle_%02d.txt"%i, "r") as f:
        lines = f.readlines()
        name = lines[0].split()[1:]
        lines = lines[1:]
        for j in range(len(lines)):
            lines[j] = lines[j].split()
            for k in range(len(lines[j])):
                lines[j][k] = float(lines[j][k])
        data = np.array(lines)
    
    v1 = data[0][4]**2 + data[0][5]**2
    v2 = data[1][4]**2 + data[1][5]**2
    print(v1, v2)
    t.append(i)
    p1_vx.append(data[0][4])
    p1_vy.append(data[0][5])
    p2_vx.append(data[1][4])
    p2_vy.append(data[1][5])

    
    fig, ax = plt.subplots(1, 1, figsize = (center_radius, center_radius))
    ax.plot(data[:, 2], data[:, 3], 'o')
    cir = plt.Circle( (center_radius, center_radius), 1, color='r', fill=False)
    ax.add_artist(cir)
    ax.set(xlim = [0, center_radius*2], ylim = [0, center_radius*2])
    plt.savefig('particle_%02d.png'%i)
    plt.close()
    

fig, ax = plt.subplots(1, 1, figsize = (center_radius, center_radius))

ax.plot(t, p1_vx, label = "P1_Vx")
ax.plot(t, p1_vy, label = "P1_Vy")
ax.plot(t, p2_vx, label = "P2_Vx")
ax.plot(t, p2_vy, label = "P2_Vy")
ax.legend(loc = 0)
plt.savefig('particle_vel.png')

