from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.pyplot as plt
import numpy as np

from collections.abc import Iterable


def colorbar(mappable, pad=0.1, side="right"):
    '''
    colorbar whose height (or width) in sync with the master axe
    https://matplotlib.org/mpl_toolkits/axes_grid/users/overview.html#colorbar-whose-height-or-width-in-sync-with-the-master-axes
    '''
    ax = mappable.axes
    fig = ax.figure
    divider = make_axes_locatable(ax)
    cax = divider.append_axes(side, size="5%", pad=pad)
    return fig.colorbar(mappable, cax=cax)


def colorbar2(mappable, shift=0.05, width=0.05, ax=None, trim_left=0, trim_right=0, side="right", label=""):
    # creates a color bar that does not shrink the main plot or panel
    # only works for horizontal bars so far

    if ax is None:
        ax = mappable.axes

    # Get current figure dimensions
    try:
        fig = ax.figure
        p = np.zeros([1,4])
        p[0,:] = ax.get_position().get_points().flatten()
    except:
        fig = ax[0].figure
        p = np.zeros([ax.size,4])
        for k, a in enumerate(ax):
            p[k,:] = a.get_position().get_points().flatten()
    xmin = np.amin(p[:,0]) ; xmax = np.amax(p[:,2]) ; dx = xmax - xmin
    ymin = np.amin(p[:,1]) ; ymax = np.amax(p[:,3]) ; dy = ymax - ymin

    if side=="top":
        cax = fig.add_axes([xmin + trim_left, ymax + shift * dy, dx - trim_left - trim_right, width * dx])
        cax.xaxis.set_ticks_position('top')
        return fig.colorbar(mappable, cax=cax, orientation="horizontal", label=label)
    elif side=="right":
        cax = fig.add_axes([xmax + shift*dx, ymin, width * dx, dy])
        cax.xaxis.set_ticks_position('top')
        return fig.colorbar(mappable, cax=cax, orientation="vertical", label=label)

def shift_axes(axes,dx,dy):

    # only 1 axis, we make it iterable
    if not isinstance(axes, Iterable):
        axes = [axes]

    for ax in axes:
        pos = ax.get_position()
        pos = [pos.x0 + dx, pos.y0 + dy, pos.width, pos.height]
        ax.set_position(pos)

def L2R(L,Teff):
  '''
  L2R(Lstar/Lsun,Teff) renvoie Rstar/Rsun
   '''
  return np.sqrt(L * (5777./Teff)**4)

def R2L(R,Teff):
  '''
  R2L(Rstar/Rsun,Teff) renvoie Lstar/Lsun
  '''
  return R**2 * (Teff/5777.)**4

def pdf(filename):
    if (filename[-4:] != ".pdf"):
        filename += ".pdf"
    print(filename)
    plt.savefig(filename, bbox_inches='tight')

def get_channel_maps(vz_field, vchannels, deltav):
    vchannels.sort()

    Nvch = len(vchannels)
    channel_maps = 100 * np.ones(vz_field.shape) # Arbitrary large value
    for i in range(len(channel_maps[:,0])):
        for j in range(len(channel_maps[0,:])):
            for ch in range(Nvch):
                if (vz_field[i,j] > vchannels[ch] - deltav and vz_field[i,j] < vchannels[ch] + deltav):
                    channel_maps[i,j] = vz_field[i,j]

    return channel_maps