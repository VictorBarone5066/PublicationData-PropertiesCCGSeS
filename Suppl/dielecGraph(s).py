from matplotlib import rcParams
from matplotlib import pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import headerOptics as h

#Figure Options
rcParams['axes.linewidth'] = 1.5
plt.rcParams["font.weight"] = "bold"
plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams["mathtext.fontset"] = "dejavuserif"
plt.tight_layout = True
plt.rcParams["figure.autolayout"] = False

WIDTH, HEIGHT = 7.2, 7.2  ##inches.  3.5 for 1 col, 7.2 for 2 col
FONTSIZE = 12
DPI = 1000

tetraFolder = "t-tOptics"
orthoFolder = "o-oOptics"

OUTCAR_LOC = {"tetra": tetraFolder + "//" + "OUTCARx",
              "ortho": orthoFolder + "//" + "OUTCARx"}
OPTICS_SAVE_LOC = "dielectric"

#color, linestyle, marker, markevery, markersize
linestyles = {"0000": ["#ff0000", "-",  None, None, 5],   #red
              "0125": ["#0058fd", "-",  's',  20,  5],   #orange
              "0250": ["#00b423", "-",  '|',  20,  10],   #lime
              "0375": ["#ff04b7", "--", 'o',  20,  5],   #d green
              "0500": ["black",   "--", '1',  20,  10],   #black
              "0625": ["#00fdea", "--", '+',  20,  10],   #cyan
              "0750": ["#93ff27", ":",  '^',  20,  5],   #blue
              "0875": ["#ffad27", ":",  'x',  20,  10],   #purple
              "1000": ["#cb34ff", ":",  'D',  20,  5]}   #pink
ALPHABET = "abcdefghijklmnopqrstuvwxyz"

def ThisIsDumb(ax_, fi, ph, ph_):
    ##Graph formatting
    if(ph == "real"):
        ax_.set_xlim(1, 3.5)
        ax_.set_ylim(5.5, 9.5)  ##for absorption
    if(ph == "imag"):
        ax.set_xlim(1, 3.5)
        ax.set_ylim(0.5, 6.0) ##for reflectivity
    if(ph_ == "tetra"):
        ax.set_title(r"$\mathrm{t\!-\!Cu_2CdGe(S_xSe_{1\!-\!x})_4}$", fontsize=FONTSIZE)
    if(ph_ == "ortho"):
        ax.set_title(r"$\mathrm{o\!-\!Cu_2CdGe(S_xSe_{1\!-\!x})_4}$", fontsize=FONTSIZE)

    for tick in ax_.xaxis.get_major_ticks():
        tick.label.set_fontsize(FONTSIZE)
    for tick in ax_.yaxis.get_major_ticks():
        tick.label.set_fontsize(FONTSIZE)

    # Legend Stuff
    labs, lins = [], []
    for key, val in linestyles.items():
        labs.append(r"$x$ = " + f'{float(key) / 1000.:.3f}')
        lins.append(Line2D([0], [0], color=val[0], linestyle=val[1], marker=val[2]))
    fi.legend(lins, labs, loc="center", bbox_to_anchor=(0.5, 0.965), fontsize=FONTSIZE, frameon=False,
               ncol=5)



# Create figure, loop through each pattern
for phaseNum, phaseId in enumerate(["tetra", "ortho"]):
    """
    *********
    REAL PART
    *********
    """
    fig, ax = plt.subplots(1, 1, sharex=True, figsize=(WIDTH, HEIGHT))
    include = list(np.arange(0.000, 1.125, 0.125))
    include = [(f'{i:.3f}'[0] + f'{i:.3f}'[2:]) for i in include]  ##gets rid of the decimal point
    for num, conc in enumerate(include):
        # Rename to make stuff easier to work with
        dieLis = h.GetDielecTensors(OUTCAR_LOC[phaseId] + conc)
        for i in range(0, len(dieLis)):
            if (dieLis[i][0].thry == "c-c"):
                if (dieLis[i][0].imag):
                    chgIm = dieLis[i][:]
                if (dieLis[i][0].real):
                    chgRe = dieLis[i][:]
            if (dieLis[i][0].thry == "d-d"):
                if (dieLis[i][0].imag):
                    denIm = dieLis[i][:]
                if (dieLis[i][0].real):
                    denRe = dieLis[i][:]


        ax.grid(True, alpha=0.5)
        ax.set_xlabel(r"$\mathrm{Energy\ (eV)}$", fontsize=FONTSIZE)
        ax.set_ylabel(r"$\epsilon_1$", fontsize=FONTSIZE)

        ##Graph the things
        x = [x_.energy for x_ in chgRe]
        y = [y_.x for y_ in chgRe]
        x, y = h.GetSplineInterp(x, y)
        if(linestyles[conc][-1] == None):
            ax.plot(x, y, label=r"$x$ = " + f'{float(conc)/1000:.3f}', color=linestyles[conc][0],
                          linestyle=linestyles[conc][1], linewidth=1.85)
        else:
            ax.plot(x, y, label=r"$x$ = " + f'{float(conc)/1000:.3f}', color=linestyles[conc][0],
                          linestyle=linestyles[conc][1], marker=linestyles[conc][2],
                          markevery=linestyles[conc][3], markersize=linestyles[conc][4], linewidth=1.85)

    # Full Plot Edits
    ThisIsDumb(ax, fig, "real", phaseId)
    plt.subplots_adjust(left=0.09, right=0.98, bottom=0.072, top=0.89, hspace=0.08, wspace=0.1)
    plt.savefig(OPTICS_SAVE_LOC + phaseId + "real.pdf")
    plt.savefig(OPTICS_SAVE_LOC + phaseId + "real.png", dpi=DPI)
    fig.clf()

    """
    *********
    IMAG PART
    *********
    """
    fig, ax = plt.subplots(1, 1, sharex=True, figsize=(WIDTH, HEIGHT))
    include = list(np.arange(0.000, 1.125, 0.125))
    include = [(f'{i:.3f}'[0] + f'{i:.3f}'[2:]) for i in include]  ##gets rid of the decimal point
    for num, conc in enumerate(include):
        # Rename to make stuff easier to work with
        dieLis = h.GetDielecTensors(OUTCAR_LOC[phaseId] + conc)
        for i in range(0, len(dieLis)):
            if (dieLis[i][0].thry == "c-c"):
                if (dieLis[i][0].imag):
                    chgIm = dieLis[i][:]
                if (dieLis[i][0].real):
                    chgRe = dieLis[i][:]
            if (dieLis[i][0].thry == "d-d"):
                if (dieLis[i][0].imag):
                    denIm = dieLis[i][:]
                if (dieLis[i][0].real):
                    denRe = dieLis[i][:]


        ax.grid(True, alpha=0.5)
        ax.set_xlabel(r"$\mathrm{Energy\ (eV)}$", fontsize=FONTSIZE)
        ax.set_ylabel(r"$\epsilon_2$", fontsize=FONTSIZE)

        ##Graph the things
        x = [x_.energy for x_ in chgRe]
        y = [y_.x for y_ in chgIm]
        x, y = h.GetSplineInterp(x, y)
        if(linestyles[conc][-1] == None):
            ax.plot(x, y, label=r"$x$ = " + f'{float(conc)/1000:.3f}', color=linestyles[conc][0],
                          linestyle=linestyles[conc][1], linewidth=1.85)
        else:
            ax.plot(x, y, label=r"$x$ = " + f'{float(conc)/1000:.3f}', color=linestyles[conc][0],
                          linestyle=linestyles[conc][1], marker=linestyles[conc][2],
                          markevery=linestyles[conc][3], markersize=linestyles[conc][4], linewidth=1.85)

    # Full Plot Edits
    ThisIsDumb(ax, fig, "imag", phaseId)
    plt.subplots_adjust(left=0.09, right=0.98, bottom=0.072, top=0.89, hspace=0.08, wspace=0.1)
    plt.savefig(OPTICS_SAVE_LOC + phaseId + "imag.pdf")
    plt.savefig(OPTICS_SAVE_LOC + phaseId + "imag.png", dpi=DPI)
    fig.clf()
