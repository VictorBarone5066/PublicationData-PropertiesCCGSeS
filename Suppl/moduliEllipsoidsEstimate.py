import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

import headerModuliEllipsoids as hEllip
import headerElas as hElas

#Figure Options
plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams["mathtext.fontset"] = "dejavuserif"
plt.tight_layout = True
plt.rcParams["figure.autolayout"] = False

WIDTH, HEIGHT = 3.5, 3.5 ##inches.  3.5 for 1 col, 7.2 for 2 col
ZOOM_FACTOR = 0.5 #No good way to determine what this should be.  Just keep messing around with it
AZIMUTH, ELEVATION = -35, 26 #degrees
FONTSIZE = 12
DPI = 1000

tetraFolder = "t-tElas"
orthoFolder = "o-oElas"

STRESSES_LOC = {"tetra": tetraFolder + "//" + "outputstresses.csv",
                "ortho": orthoFolder + "//" + "outputstresses.csv"}
SAVE_LOC = "C://Users//baron//Desktop//allEllips"

COLORS = ["red",
          "green",
          "blue",
          "black"]

paramDict = {}
for phaseId in ["tetra", "ortho"]:
    data = hElas.ParseOutfile(STRESSES_LOC[phaseId])

    for conc, set in data.items():
        if(conc == 0.0 or conc == 1.0):
            # Data from strains in only the x, y, z directions are stored in idNums 1, 2, and 3.
            youngX = hEllip.GetDirYoungModEst(set[1])
            youngY = hEllip.GetDirYoungModEst(set[2])
            youngZ = hEllip.GetDirYoungModEst(set[3])
            print(phaseId, conc, youngX, youngY, youngZ)
            paramDict[phaseId + str(conc)] = [youngX, youngY, youngZ]

#Now, The dict is {conc: [Yx, Yy, Yz], ...}.  Get axis params
##Plot Setup
fig = plt.figure()#figsize=(WIDTH, HEIGHT))
ax = fig.add_subplot(projection='3d')
paramList = list(paramDict.values())
paramList = hEllip.SortByRadius(paramList)

ax = hEllip.GetOrthoVects(ax, paramList, labs=[r"$Y_\mathrm{a}'$", r"$Y_\mathrm{b}'$", r"$Y_\mathrm{c}'$"],
                          offset=0.15)
alphas = hEllip.GetAlphaList(paramList, alMax=0.6, alMin=0.4)
strides = hEllip.GetStrideList(paramList, stMax=9, stMin=2)

##Plot all ellipsoids...
for ind, param in enumerate(paramList):
    x, y, z = hEllip.GetEllipsoidParams(param[0], param[1], param[2])
    ax.plot_wireframe(np.array(x), np.array(y), np.array(z), rstride=strides[ind], cstride=strides[ind],
                      color=COLORS[ind], alpha=alphas[ind])

#Customization
ax.set_axis_off()
hEllip.EqualAxis(ax)

ax.set_xlim(-ZOOM_FACTOR, ZOOM_FACTOR)
ax.set_ylim(-ZOOM_FACTOR, ZOOM_FACTOR)
ax.set_zlim(-ZOOM_FACTOR, ZOOM_FACTOR)
ax.view_init(azim=AZIMUTH, elev=ELEVATION)

fig.show()
plt.savefig(SAVE_LOC + ".pdf")
#plt.clf()
