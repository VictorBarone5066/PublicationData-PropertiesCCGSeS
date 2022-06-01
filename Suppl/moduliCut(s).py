from matplotlib import rcParams
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np

import headerModuliCuts as hCuts
import headerElas as hElas

#Figure Options
rcParams['axes.linewidth'] = 1.5
plt.rcParams["font.weight"] = "bold"
plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams["mathtext.fontset"] = "dejavuserif"
plt.tight_layout = True
plt.rcParams["figure.autolayout"] = False

WIDTH, HEIGHT = 7.2, 7.2 ##inches.  3.5 for 1 col, 7.2 for 2 col
FONTSIZE = 12
DPI = 1000
SAMPLING = 10000 ##total number of sampling points around the circle

tetraFolder = "t-tElas"
orthoFolder = "o-oElas"
CIJSIJ_LOC = {"tetra": tetraFolder + "//" + "tetraCijsSijs.csv",
              "ortho": orthoFolder + "//" + "orthoCijsSijs.csv"}
SAVE_LOC = "modcut"

# Linestyles (highest contrast that I could manage to make)
#color, linestyle, marker, markevery, markersize
linestyles = {0.000: ["#ff0000", "-",  None, None, 5],   #red
              0.125: ["#0058fd", "-",  's',  200,  5],   #orange
              0.250: ["#00b423", "-",  '|',  200,  10],   #lime
              0.375: ["#ff04b7", "--", 'o',  200,  5],   #d green
              0.500: ["black",   "--", '1',  200,  10],   #black
              0.625: ["#00fdea", "--", '+',  200,  10],   #cyan
              0.750: ["#93ff27", ":",  '^',  200,  5],   #blue
              0.875: ["#ffad27", ":",  'x',  200,  10],   #purple
              1.000: ["#cb34ff", ":",  'D',  200,  5]}   #pink
ALPHABET = "abcdefghijklmnopqrstuvwxyz"

def Bro(a, f, na, xl, yl, pid):
    na = na[1:]
    if(pid == "tetra"):
        na = r"$\mathrm{t}$" + na
    if(pid == "ortho"):
        na = r"$\mathrm{o}$" + na
    #a.set_xlim(-110, 110)
    #a.set_ylim(-110, 110)
    #a.set_xticks([-110, -55, 0, 55, 110])
    #a.set_yticks([-110, -55, 0, 55, 110])
    a.grid(True, alpha=0.5)
    for tick in a.xaxis.get_major_ticks():
        tick.label.set_fontsize(FONTSIZE)
    for tick in a.yaxis.get_major_ticks():
        tick.label.set_fontsize(FONTSIZE)

    a.set_title(na, fontsize=FONTSIZE)
    a.set_xlabel(xl, fontsize=FONTSIZE)
    a.set_ylabel(yl, fontsize=FONTSIZE)

# Create figure, loop through each pattern
for phaseNum, phaseId in enumerate(["tetra", "ortho"]):
    data = hElas.ParseCijSijFile(CIJSIJ_LOC[phaseId])

    """
    ***********
    NORMAL TO X
    ***********
    """
    figX, axX = plt.subplots(1, 1, sharex=False, sharey=False, figsize=(WIDTH, HEIGHT))
    Bro(axX, figX, phaseId[0] + r"$\mathrm{\!-\!Cu_2CdGe(S_xSe_{1\!-\!x})_4}: (1, 0, 0)}$", r"$Y_\mathrm{b}$ $\mathrm{(GPa)}$",
        r"$Y_\mathrm{c}$ $\mathrm{(GPa)}$", phaseId)
    for num, dat in enumerate(data):
        ##Get directional youngs mods in lists
        _, y, z = hCuts.GetCircleSampling(normal='x', ptDen=SAMPLING)
        for i in range(0, len(y)):
            thisYMod = hCuts.GetDirYoungsMod(0, y[i], z[i], dat.sij)
            y[i] *= thisYMod
            z[i] *= thisYMod
        ##Plot them
        axX.plot(y, z, linewidth=1.85, color=linestyles[dat.x][0], linestyle=linestyles[dat.x][1],
                marker=linestyles[dat.x][2], markevery=linestyles[dat.x][3],
                label=r"$x$ = " + f'{dat.x:.3f}', markersize=linestyles[dat.x][4])
    figX.legend(loc="center", bbox_to_anchor=(0.5, 0.955), fontsize=FONTSIZE, frameon=False, ncol=5)
    plt.subplots_adjust(left=0.1, right=0.96, bottom=0.093, top=0.88, hspace=0.1, wspace=0.2)
    plt.savefig(SAVE_LOC + phaseId + "100.pdf")
    plt.savefig(SAVE_LOC + phaseId + "100.png", dpi=DPI)
    figX.clf()

    """
    ***********
    NORMAL TO Y
    ***********
    """
    figY, axY = plt.subplots(1, 1, sharex=False, sharey=False, figsize=(WIDTH, HEIGHT))
    Bro(axY, figY, phaseId[0] + r"$\mathrm{\!-\!Cu_2CdGe(S_xSe_{1\!-\!x})_4}: (0, 1, 0)}$", r"$Y_\mathrm{a}$ $\mathrm{(GPa)}$",
        r"$Y_\mathrm{c}$ $\mathrm{(GPa)}$", phaseId)
    for num, dat in enumerate(data):
        ##Get directional youngs mods in lists
        x, _, z = hCuts.GetCircleSampling(normal='y', ptDen=SAMPLING)
        for i in range(0, len(x)):
            thisYMod = hCuts.GetDirYoungsMod(x[i], 0, z[i], dat.sij)
            x[i] *= thisYMod
            z[i] *= thisYMod
        ##Plot them
        axY.plot(x, z, linewidth=1.85, color=linestyles[dat.x][0], linestyle=linestyles[dat.x][1],
                marker=linestyles[dat.x][2], markevery=linestyles[dat.x][3],
                label=r"$x$ = " + f'{dat.x:.3f}', markersize=linestyles[dat.x][4])
    figY.legend(loc="center", bbox_to_anchor=(0.5, 0.955), fontsize=FONTSIZE, frameon=False, ncol=5)
    plt.subplots_adjust(left=0.1, right=0.96, bottom=0.093, top=0.88, hspace=0.1, wspace=0.2)
    plt.savefig(SAVE_LOC + phaseId + "010.pdf")
    plt.savefig(SAVE_LOC + phaseId + "010.png", dpi=DPI)
    figY.clf()

    """
    ***********
    NORMAL TO Z
    ***********
    """
    figZ, axZ = plt.subplots(1, 1, sharex=False, sharey=False, figsize=(WIDTH, HEIGHT))
    Bro(axZ, figZ, phaseId[0] + r"$\mathrm{\!-\!Cu_2CdGe(S_xSe_{1\!-\!x})_4}: (0, 0, 1)}$", r"$Y_\mathrm{a}$ $\mathrm{(GPa)}$",
        r"$Y_\mathrm{b}$ $\mathrm{(GPa)}$", phaseId)
    for num, dat in enumerate(data):
        ##Get directional youngs mods in lists
        x, y, _ = hCuts.GetCircleSampling(normal='z', ptDen=SAMPLING)
        for i in range(0, len(x)):
            thisYMod = hCuts.GetDirYoungsMod(x[i], y[i], 0, dat.sij)
            x[i] *= thisYMod
            y[i] *= thisYMod
        ##Plot them
        axZ.plot(x, y, linewidth=1.85, color=linestyles[dat.x][0], linestyle=linestyles[dat.x][1],
                 marker=linestyles[dat.x][2], markevery=linestyles[dat.x][3],
                 label=r"$x$ = " + f'{dat.x:.3f}', markersize=linestyles[dat.x][4])
    figZ.legend(loc="center", bbox_to_anchor=(0.5, 0.955), fontsize=FONTSIZE, frameon=False, ncol=5)
    plt.subplots_adjust(left=0.1, right=0.96, bottom=0.093, top=0.88, hspace=0.1, wspace=0.2)
    plt.savefig(SAVE_LOC + phaseId + "001.pdf")
    plt.savefig(SAVE_LOC + phaseId + "001.png", dpi=DPI)
    figZ.clf()
