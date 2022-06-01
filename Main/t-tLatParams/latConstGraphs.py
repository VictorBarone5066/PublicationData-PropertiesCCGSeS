from matplotlib import pyplot as plt
import matplotlib.lines as lin
import numpy as np
import outputParseBase as b

plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams["mathtext.fontset"] = "dejavuserif"
plt.tight_layout = True

INFILE_LOC = "convOutput.csv"
LATTICE_PARAM_SAVE_LOC = "latParams"

#Takes the work directory name and returns the x "concentration" value
def NameToX(s):
    return float(s[1:]) / 1000

#Get all a, b, c values compared to concentration
dataLines = b.ParseOutfile(INFILE_LOC)
a, b, c, x = [], [], [], []
for d in dataLines:
    a.append(d.aV)
    b.append(d.bV)
    c.append(d.cV)
    x.append(NameToX(d.wDirectory))
#Special cases:  correct the c value and name for the 32-atom cells.  Necessary because I was lazy 
#earlier on in this project
for i in range(0, int((len(c) + len(x)) / 2)):
    if(x[i] > 1.0):
        x[i] /= 10
        c[i] /= 2

#------Plot a, b, c vs. x---------------------------------------------------------------------------
f, (axC, axAB) = plt.subplots(2, 1, sharex=True)

#x and y lables, limits
dx, dy = 0.02, 0.15 #angstroms
abMin = min([min([a, b]) for a, b in zip(a, b)])
abMax = max([max([a, b]) for a, b in zip(a, b)])

axAB.set_xticks([0.0, 0.25, 0.5, 0.75, 1.0])
axAB.set_yticks(np.arange(5.8, 6.9, 0.3))
axC.set_yticks(np.arange(8.2, 12.2, 1.0))

axAB.set_xlim(0 - dx, 1 + dx)
axAB.set_ylim(abMin - dy, abMax + dy)
axC.set_ylim(min(c) - dy, max(c) + dy)

for tick in axAB.xaxis.get_major_ticks():
    tick.label.set_fontsize(11) 
for tick in axAB.yaxis.get_major_ticks():
    tick.label.set_fontsize(11)
for tick in axC.yaxis.get_major_ticks():
    tick.label.set_fontsize(11)

#Hide borders between axes to give a broken-graph style look
axC.spines['bottom'].set_visible(False)
axAB.spines['top'].set_visible(False)
axAB.xaxis.tick_bottom()
axC.tick_params(axis='x', which="both", top=False, bottom=False)
#Include diagonal lines to represent axis break.  
#From https://matplotlib.org/3.1.0/gallery/subplots_axes_and_figures/broken_axis.html
d = .015
kwargs = dict(transform=axC.transAxes, color='k', clip_on=False, linewidth = 0.8)
axC.plot((-d, +d), (-d, +d), **kwargs)        # top-left diagonal
axC.plot((1 - d, 1 + d), (-d, +d), **kwargs)  # top-right diagonal
kwargs.update(transform=axAB.transAxes)  # switch to the bottom axes
axAB.plot((-d, +d), (1 - d, 1 + d), **kwargs)  # bottom-left diagonal
axAB.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)  # bottom-right diagonal

#Legend
mx = lin.Line2D([], [], color='m', marker='x', linewidth=0, markersize=10)
kp = lin.Line2D([], [], color='k', marker='+', linewidth=0, markersize=10)
rx = lin.Line2D([], [], color='r', marker='x', linewidth=0, markersize=10)
leg = axC.legend([kp, rx, mx], [r'$\vec{a}$', r'$\vec{b}$', r'$\vec{c}$'], loc="upper right", 
                 frameon=True, fontsize = 11)

#Lables
axAB.set_xlabel(r"$x=[\mathrm{S}]\,/\,\left([\mathrm{S}]+[\mathrm{Se}]\right)$", fontsize = 11)
f.text(0.03, 0.5, r"Lattice Parameter (Å)", rotation = 90, 
       fontsize = 12, verticalalignment='center')

#Plot the things, save
axAB.plot(x, a, 'k+', markersize = 8)
axAB.plot(x, b, 'rx')
axC.plot(x, c, 'mx')
f.savefig(LATTICE_PARAM_SAVE_LOC + ".pdf")
f.savefig(LATTICE_PARAM_SAVE_LOC + ".png", dpi=1000)
f.clf()

