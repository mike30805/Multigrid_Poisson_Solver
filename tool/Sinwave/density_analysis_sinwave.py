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

def plot_profile( filename, savename):

    data = np.loadtxt(filename)
    plot_title = filename[:-4]
    size = data.shape[1]
    
    dens = np.zeros(size)

    #2-dim case
    if data.shape[0] == data.shape[1]:
        dens = data[int(size/2),:]
        r = np.linspace(0, size, size)
        
    #3-dim case
    else: 
        data_temp = np.reshape(data, (size, size, size))
        dens = data_temp[int(size/2),int(size/2),:]
        r = np.linspace(0, size, size)

    fig, ax = plt.subplots(1, 1, figsize = (10, 10))
    ax.plot(r, dens)
    ax.set( title = plot_title, xlabel = 'x', ylabel = 'y' )

    
    if savename != "":
        plt.savefig(savename)
    plt.close()

    return


N = 1
n_bin = 20
ana_obj = "Density"

for i in range(N):
    
    #Plot 2D or 3D density
    filename = ana_obj+"_%02d.txt"%i
    savename = filename[:-4]+".png"
    #plot_data( filename, savename, True )
    plot_data( filename, savename, False )
    
    #Plot density profile along the middle line
    savename = filename[:-4]+"_profile.png"
    plot_profile( filename, savename )
