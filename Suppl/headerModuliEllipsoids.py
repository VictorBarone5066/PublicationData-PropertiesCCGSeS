import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

pi = 3.141592654
sin = np.sin
cos = np.cos

#eV per cubic angstrom to GPa
def EvAngstToGPa(v):
    return v*160.2176487

#Makes 3D axis all equal.  The fact that I have to do this myself instead of matplotlib doing it automatically
#is completly insane
#Input: matplotlib axis (like from pyplot.plt.gca())
def EqualAxis(ax):
    xLims = ax.get_xlim3d()
    yLims = ax.get_ylim3d()
    zLims = ax.get_zlim3d()

    xRange = abs(xLims[1] - xLims[0])
    xMid = np.mean(xLims)
    yRange = abs(yLims[1] - yLims[0])
    yMid = np.mean(yLims)
    zRange = abs(zLims[1] - zLims[0])
    zMid = np.mean(zLims)

    plRad = max([xRange, yRange, zRange])/2.
    ax.set_xlim3d([xMid - plRad, xMid + plRad])
    ax.set_ylim3d([yMid - plRad, yMid + plRad])
    ax.set_zlim3d([zMid - plRad, zMid + plRad])

    return ax

#Expects a list of a list of triples in the form pLis=[[x00, y00, z00], [x01, y01, z01], ...],
#                                                      [x10, y10, z10], [x11, y11, z11], ...]]
#Where the nth list is the nth surface to plot and the mth value of that list is an xyz pair describing a pt
#Optional: a list of plot labels as [xLabel, yLabel, zLabel] and offset in units of a, b, c
def GetOrthoVects(ax, pLis, labs=['a', 'b', 'c'], offset=0.5, fontsize_=12):
    absMaxX, absMaxY, absMaxZ = -1, -1, -1
    for all in pLis:
        for i in range(0, len(all[0])):
            if (abs(all[0][i]) > absMaxX):
                absMaxX = abs(all[0][i])
            if (abs(all[1][i]) > absMaxY):
                absMaxY = abs(all[1][i])
            if (abs(all[2][i]) > absMaxZ):
                absMaxZ = abs(all[2][i])

    ##Set axis labels up to meet bounds of vectors
    ax.text(absMaxX+offset, 0, 0, labs[0], color='k', fontsize=fontsize_) ##a direction
    ax.plot([-absMaxX, absMaxX], [0, 0], [0, 0], color='k')
    ax.text(0, absMaxY+offset, 0, labs[1], color='k', fontsize=fontsize_) ##b direction
    ax.plot([0, 0], [-absMaxY, absMaxY], [0, 0], color='k')
    ax.text(0, 0, absMaxZ+offset, labs[2], color='k', fontsize=fontsize_) ##c direction
    ax.plot([0, 0], [0, 0], [-absMaxZ, absMaxZ], color='k')

    return ax

#Gets the direction cosines for a vector defined by u*xhat + v*yhat + w*zhat
def GetDirCosines(u, v, w):
    dist = (u**(2) + v**(2) + w**(2))**(1./2.)
    #cos alpha: between vec(u, v, w) and x direction, cos beta, cos gamma
    return u/dist, v/dist, w/dist

#Make sure sij matrix is indexed such that s[1][1] = s11 instead of s[0][0] = s11.
#ref: Crystal Structure and Mechanical Properties of ThBC2 by Xinchun Zhou and Baobing Zheng
def GetDirYoungsMod(u, v, w, s):
    al, be, ga = GetDirCosines(u, v, w)
    return 1./(s[1][1]*al**(4) + s[2][2]*be**(4) + s[3][3]*ga**(4) + \
               (2*s[1][2] + s[6][6])*al**(2)*be**(2) + \
               (2*s[2][3] + s[4][4])*be**(2)*ga**(2) + \
               (2*s[1][3] + s[5][5])*al**(2)*ga**(2))

#Gives a set of x, y, z points that sample the surface of a unit sphere (i.e. r^2 = x^2 + y^2 + z^2)
def GetSphereSampling(rad=1., ptDen=50):
    x, y, z = [], [], []
    for theta in list(np.linspace(0, pi, int(np.ceil(ptDen/2.)))):
        for phi in list(np.linspace(0, 2*pi, ptDen)):
            x.append(rad*sin(theta)*cos(phi))
            y.append(rad*sin(theta)*sin(phi))
            z.append(rad*cos(theta))
    return x, y, z
