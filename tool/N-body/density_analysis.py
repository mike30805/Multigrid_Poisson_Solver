import numpy as np
import matplotlib.pylab as plt
from matplotlib.colors import LogNorm

def plot_matrix( matrix, plot_title, savename="" ):
    fig, ax = plt.subplots(1, 1, figsize = (10, 10))
    #im0 = ax[1].imshow(data, origin = 'lower', norm=LogNorm(vmin=1.e-1, vmax=1.))
    im0 = ax.imshow(matrix, origin = 'lower')
    ax.set( title = plot_title, xlabel = 'x', ylabel = 'y' )
    fig.colorbar(im0, ax = ax)

    if savename != "":
        plt.savefig(savename)

    plt.close()

    return



def plot_data( filename, savename = "", series_3D = False ):
    data = np.loadtxt(filename)
    title = filename[:-4]

    #2-dim case
    if data.shape[0] == data.shape[1]:
        plot_matrix( data, title, savename )

    #3-dim case
    else:
        size = data.shape[1]
        data_temp = np.reshape(data, (size, size, size))
        if series_3D:
            for i in range(size):
                if savename == "":
                    save_new = ""
                else:
                    save_new = savename[:-4]
                    save_new = save_new+"_%03d.png"%i
                    plot_matrix( data_temp[i], title, save_new )
        else:
            plot_matrix( data_temp[int(size/2)], title, savename )

    return

def plot_profile( filename, savename, center=[], n_bin=100 ):
    #center is the center cell indices
    if center == []:
        print("Center need to be given!")
        return

    data = np.loadtxt(filename)
    plot_title = filename[:-4]
    
    size = data.shape[1]
    
    dens = np.zeros(n_bin)
    count = np.zeros(n_bin)

    #2-dim case
    if data.shape[0] == data.shape[1]:
        r = np.linspace(0, size*np.sqrt(2), n_bin)
        dr = r[1]-r[0]

        for i in range(size):
            for j in range(size):
                dist = np.sqrt( (center[0]-i)**2 + (center[1]-j)**2 )
                r_idx = int( dist/dr )
                dens[r_idx] += data[i, j]
                count[r_idx] += 1
    #3-dim case
    else: 
        data_temp = np.reshape(data, (size, size, size))
        r = np.linspace(0, size*np.sqrt(3), n_bin)
        dr = r[1]-r[0]

        for i in range(size):
            for j in range(size):
                for k in range(size):
                    dist = np.sqrt( (center[0]-i)**2 + (center[1]-j)**2 + (center[2]-k)**2 )
                    r_idx = int( dist/dr )
                    dens[r_idx] += data_temp[i, j, k]
                    count[r_idx] += 1

    dens = dens / count

    fig, ax = plt.subplots(1, 1, figsize = (10, 10))
    ax.plot(r, dens)
    ax.set( title = plot_title, xlabel = 'x', ylabel = 'y' )

    
    if savename != "":
        plt.savefig(savename)
    plt.close()

    return


N = 1
n_bin = 20
data_center = [10, 10, 10]
ana_obj = "Density"

for i in range(N):
    print(i)
    filename = ana_obj+"_%02d.txt"%i
    savename = filename[:-4]+".png"

    #plot_data( filename, savename, True )
    plot_data( filename, savename, False )
    savename = filename[:-4]+"_profile.png"
    plot_profile( filename, savename, data_center, n_bin )
