from matplotlib import rcParams
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

import headerModuliEllipsoids as hEllip
import headerElas as hElas

#Figure Options
rcParams['axes.linewidth'] = 1.5
plt.rcParams["font.weight"] = "bold"
plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams["mathtext.fontset"] = "dejavuserif"
plt.tight_layout = True
plt.rcParams["figure.autolayout"] = False

WIDTH, HEIGHT = 3.5, 3.5 ##inches.  3.5 for 1 col, 7.2 for 2 col
ZOOM_FACTOR = 65 #No good way to determine what this should be.  Just keep messing around with it
AZIMUTH, ELEVATION = -55, 10 #degrees
#AZIMUTH, ELEVATION = -55, 10 #degrees
FONTSIZE = 12
DPI = 3000

TYPE = "tetra" ##tetra or ortho


tetraFolder = "t-tElas"
orthoFolder = "o-oElas"

CIJSIJ_LOC = {"tetra": tetraFolder + "//" + "tetraCijsSijs.csv",
              "ortho": orthoFolder + "//" + "orthoCijsSijs.csv"}
SAVE_LOC = TYPE[0] + "0vs1"

#phase id, concentration, color
INCLUDE = [[TYPE, 0.0, "#ea0000"], #red
           [TYPE, 1.0, "#0080FF"]] #blue

#Get data in a good format
allData = []
xSample, ySample, zSample = hEllip.GetSphereSampling(ptDen=300) ##Unit sphere
for index, inc in enumerate(INCLUDE):
    phaseId, conc, color = inc[0], inc[1], inc[2]
    data = hElas.ParseCijSijFile(CIJSIJ_LOC[phaseId])

    for num, dat in enumerate(data):
        if(dat.x == conc):
            x_, y_, z_ = xSample[:], ySample[:], zSample[:]
            for i in range(0, len(xSample)):
                thisMod = hEllip.GetDirYoungsMod(xSample[i], ySample[i], zSample[i], data[num].sij)
                x_[i] *= thisMod
                y_[i] *= thisMod
                z_[i] *= thisMod

            allData.append([x_, y_, z_, color])

#Plot Stuff
fig = plt.figure(figsize=(WIDTH, HEIGHT))
ax = fig.add_subplot(projection="3d")
#offset 40 for tetra and 15 for ortho works well
ax = hEllip.GetOrthoVects(ax, pLis=allData, offset=10, labs=[r"$Y_\mathrm{a}$", r"$Y_\mathrm{b}$",
                          r"$Y_\mathrm{c}$"], fontsize_=FONTSIZE)
for num, dat in enumerate(allData):
    ax.scatter3D(dat[0], dat[1], dat[2], color=dat[3], s=0.01)


#ax.set_axis_off()
#ax.set_xlim(-ZOOM_FACTOR, ZOOM_FACTOR)
#ax.set_ylim(-ZOOM_FACTOR, ZOOM_FACTOR)
#ax.set_zlim(-ZOOM_FACTOR, ZOOM_FACTOR)
ax.view_init(azim=AZIMUTH, elev=ELEVATION)
ax = hEllip.EqualAxis(ax)

#fig.show()
fig.savefig(SAVE_LOC + ".png", dpi=DPI)
fig.savefig(SAVE_LOC + ".pdf")
plt.clf()
