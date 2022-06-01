from matplotlib import rcParams
import matplotlib.pyplot as plt
import matplotlib.patches as mpat
import numpy as np
from scipy.signal import savgol_filter
import headerDos as dos

#Figure Options
rcParams['axes.linewidth'] = 1.5
plt.rcParams["font.weight"] = "bold"
plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams["mathtext.fontset"] = "dejavuserif"
plt.tight_layout = True
plt.rcParams["figure.autolayout"] = False

WIDTH, HEIGHT = 7.2, 8.0 ##inches.  3.5 for 1 col, 7.2 for 2 col
FONTSIZE = 12
DPI = 1000

tetraFolder = "t-tPdos"
orthoFolder = "o-oPdos"

DOSCAR_LOC = {"tetra": tetraFolder + "//" + "DOSCARx",
              "ortho": orthoFolder + "//" + "DOSCARx"}
POSCAR_LOC = {"tetra": tetraFolder + "//" + "POSCARx",
              "ortho": orthoFolder + "//" + "POSCARx"}
PDOS_SAVE_LOC = "pdos"


SAVGOL_WINSIZE = 11
SAVGOL_POLY = 3

#Define colors to use so that they are consistent throughout the large graph.  Ditto Hatches
ELEMS = {"Cu": ["#0e9400", "-"], #green
         "Cd": ["#ea0000", (0, (3, 5, 1, 5, 1, 5))], #red
         "Ge": ["#00008b", ":"], #blue
         'S' : ["#e127f1", "-."], #purple
         "Se": ["#ff8a00", "--"]} #orange
ALPHABET = "abcdefghijklmnopqrstuvwxyz"


#Create a figure, loop through each pattern
for phaseNum, phaseId in enumerate(["tetra", "ortho"]):
    fig, ax = plt.subplots(9, 1, sharex=False, figsize=(WIDTH, HEIGHT))
    include = list(np.arange(0.000, 1.125, 0.125))
    include = [(f'{i:.3f}'[0] + f'{i:.3f}'[2:]) for i in include] ##gets rid of the decimal point
    for num, conc in enumerate(include):

        #Setup reading from DOSCAR
        ranges = dos.ScanPoscar(POSCAR_LOC[phaseId] + conc)
        nedos = dos.GetNEDOS(DOSCAR_LOC[phaseId] + conc)

        #Get atoms in a dictionary, sum the (s, p, d) elements to make the graph managable to read
        info = {}
        for r in ranges:
            mi = dos.AtomGroup(DOSCAR_LOC[phaseId] + conc, r[1], nedos, atomType = r[0])
            info[mi.atomType] = [mi, [s + p + d for s, p, d in zip(mi.sDosSum, mi.pDosSum,
                mi.dDosSum)]]

        #Individual plot customizations
        ax[num].set_xlim(-3.0, 4.0)
        ax[num].set_ylim(0.0, 6.0)
        ax[num].set_xticks(np.arange(-3.0, 5.0, 1.0))
        ax[num].set_yticks(np.arange(0.0, 9.0, 3.0))

        #Place concentration and ref letter on appropriate plot
        ax[num].text(3.0, 4.60, r"$x=$" + f'{float(conc)/1000:.3f}', fontsize=FONTSIZE)

        #X and Y Labels
        if(num == 8):
            ax[num].set_xlabel(r"$E-E_\mathrm{F}$ $\mathrm{(eV)}$", fontsize=FONTSIZE)

        ##Tick settings
        for tick in ax[num].xaxis.get_major_ticks():
            tick.label.set_fontsize(FONTSIZE)
        for tick in ax[num].yaxis.get_major_ticks():
            tick.label.set_fontsize(FONTSIZE)
        ###Really dumb but necessary workarounds to force gridlines but not tick markers
        if(num != 8):
            ax[num].tick_params(axis='x', colors="white")
            for tick in ax[num].xaxis.get_major_ticks():
                tick.label.set_fontsize(0.01)

        #Plot DOS
        for key, val in info.items():
            scaledEnergy = [e - val[0].dosData[0]["eFermi"] for e in val[0].energy]
            thisDos = savgol_filter([v - 0.055 for v in val[1]], SAVGOL_WINSIZE, SAVGOL_POLY)
            if(num == 3):
                ax[num].plot(scaledEnergy, thisDos, label=key, color=ELEMS[key][0], linewidth=1.85,
                             linestyle=ELEMS[key][1])
            else:
                ax[num].plot(scaledEnergy, thisDos, color=ELEMS[key][0], linewidth=1.85,
                             linestyle=ELEMS[key][1])

        ax[num].axvline(x=0, color='k', linestyle="--", linewidth=1.0)
        ax[num].grid(True, alpha=0.5)

    #Main plot edits, save and finish
    plt.subplots_adjust(left=0.08, right=0.98, bottom=0.07, top=0.96, hspace=0.23, wspace=0.0)
    ##Legend stuff
    fig.legend(loc="center", ncol=len(ELEMS), bbox_to_anchor=(0.5, 0.98),
               fontsize=FONTSIZE, frameon=False)

    #Additional Labels
    ylabel=r"$\mathrm{Electronic\ LDOS\ (states/eV)}$"
    fig.text(0.025, 0.5, ylabel, rotation=90, fontsize=FONTSIZE, transform=fig.transFigure,
             verticalalignment="center", horizontalalignment="center")

    plt.savefig(PDOS_SAVE_LOC + phaseId + ".pdf")
    plt.savefig(PDOS_SAVE_LOC + phaseId + ".png", dpi=1000)
    fig.clf()
