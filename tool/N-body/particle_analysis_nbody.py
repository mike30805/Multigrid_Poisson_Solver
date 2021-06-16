import numpy as np
import matplotlib.pylab as plt
from matplotlib.colors import LogNorm

N = 1
center_radius = 3.0

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

    fig, ax = plt.subplots(1, 1, figsize = (5, 5))
    ax.plot(data[:, 2], data[:, 3], 'o')

    ax.set(xlim = [0, center_radius*2], ylim = [0, center_radius*2])
    plt.savefig('particle_%02d.png'%i)
    plt.close()
    