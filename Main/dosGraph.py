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

WIDTH, HEIGHT = 7.2, 7.2 ##inches.  3.5 for 1 col, 7.2 for 2 col
FONTSIZE = 12
DPI = 1000

tetraFolder = "t-tPdos"
orthoFolder = "o-oPdos"

DOSCAR_LOC = {"tetra": tetraFolder + "//" + "DOSCARx",
              "ortho": orthoFolder + "//" + "DOSCARx"}
POSCAR_LOC = {"tetra": tetraFolder + "//" + "POSCARx",
              "ortho": orthoFolder + "//" + "POSCARx"}
PDOS_SAVE_LOC = "allDos"


SAVGOL_WINSIZE = 11
SAVGOL_POLY = 3

#Define colors to use so that they are consistent throughout the large graph.  Ditto Hatches
ELEMS = {"Cu": ["#ea0000", "-"], #red
         "Ge": ["#00008b", ":"], #blue
         'S' : ["#90ee90", "--"], #purple
         "Se": ["#013220", "--"], #d green
         "Cd": ["cyan",   "-."]} #not orange
ALPHABET = "ab.......cd..............." #lol

def IndexToIJ(conc, phasId):
    if(phasId[0] == 't'):
        if(conc[0] == '0'):
            return 0, 0
        return 1, 0
    if(conc[0] == '0'):
        return 0, 1
    return 1, 1

#Create a figure, loop through each pattern
fig, ax = plt.subplots(2, 2, sharex=False, figsize=(WIDTH, HEIGHT))
for phaseNum, phaseId in enumerate(["tetra", "ortho"]):
    include = [0.0, 1.0]
    include = [(f'{i:.3f}'[0] + f'{i:.3f}'[2:]) for i in include] ##gets rid of the decimal point
    for num, conc in enumerate(include):
        i, j = IndexToIJ(conc, phaseId)

        #Setup reading from DOSCAR
        ranges = dos.ScanPoscar(POSCAR_LOC[phaseId] + conc)
        nedos = dos.GetNEDOS(DOSCAR_LOC[phaseId] + conc)

        #Get atoms in a dictionary, sum the (s, p, d) elements to make the graph managable to read
        info = {}
        for r in ranges:
            if(r[0] != "I dont even know why this if statement exists"):
                mi = dos.AtomGroup(DOSCAR_LOC[phaseId] + conc, r[1], nedos, atomType = r[0])
                info[mi.atomType] = [mi, [s + p + d for s, p, d in zip(mi.sDosSum, mi.pDosSum,
                    mi.dDosSum)]]

        #Individual plot customizations
        ax[i][j].set_xlim(-3.0, 4.0)
        ax[i][j].set_ylim(0.0, 6.0)
        ax[i][j].set_xticks(np.arange(-3.0, 5.0, 1.0))
        ax[i][j].set_yticks(np.arange(0.0, 8.0, 2.0))

        #Place concentration and ref letter on appropriate plot
        ax[i][j].text(-2.824, 5.62, ALPHABET[num + 9*phaseNum] + ')', fontsize=FONTSIZE,
                      backgroundcolor="white")

        ##Tick settings
        for tick in ax[i][j].xaxis.get_major_ticks():
            tick.label.set_fontsize(FONTSIZE)
        for tick in ax[i][j].yaxis.get_major_ticks():
            tick.label.set_fontsize(FONTSIZE)
        ###Really dumb but necessary workarounds to force gridlines but not tick markers
        if(i != 1):
            ax[i][j].tick_params(axis='x', colors="white")
            for tick in ax[i][j].xaxis.get_major_ticks():
                tick.label.set_fontsize(0.01)
        if(j != 0):
            ax[i][j].tick_params(axis='y', colors="white", labelsize=0.1)
            for tick in ax[i][j].yaxis.get_major_ticks():
                tick.label.set_fontsize(0.01)

        #Plot DOS -0.055
        for key, val in info.items():
            scaledEnergy = [e - val[0].dosData[0]["eFermi"] for e in val[0].energy]
            thisDos = savgol_filter(val[1], SAVGOL_WINSIZE, SAVGOL_POLY)
            ax[i][j].plot(scaledEnergy, thisDos, color=ELEMS[key][0], linestyle=ELEMS[key][1], linewidth=1.85)

        ax[i][j].axvline(0, color='k', linewidth=1.5, linestyle='--')
        ax[i][j].grid(True, alpha=0.5)

#Main plot edits, save and finish
plt.subplots_adjust(left=0.08, right=0.96, bottom=0.07, top=0.985, hspace=0.1, wspace=0.1)
##Legend stuff
for key, val in ELEMS.items():
    for i in range(0, 2):
        for j in range(0, 2):
            ax[i][j].plot(0, -20, label=key, color=ELEMS[key][0], linestyle=ELEMS[key][1])

ax[0][0].legend(title=r"$\mathrm{t\!-\!Cu_2CdGeSe_4}$", title_fontsize=FONTSIZE-2, loc="upper right", fontsize=FONTSIZE-2, ncol=2)
ax[1][0].legend(title=r"$\mathrm{t\!-\!Cu_2CdGeS_4}$", title_fontsize=FONTSIZE-2, loc="upper right", fontsize=FONTSIZE-2, ncol=2)
ax[0][1].legend(title=r"$\mathrm{o\!-\!Cu_2CdGeSe_4}$", title_fontsize=FONTSIZE-2, loc="upper right", fontsize=FONTSIZE-2, ncol=2)
ax[1][1].legend(title=r"$\mathrm{o\!-\!Cu_2CdGeS_4}$", title_fontsize=FONTSIZE-2, loc="upper right", fontsize=FONTSIZE-2, ncol=2)

#fig.legend(loc="center", ncol=len(ELEMS), bbox_to_anchor=(0.5, 0.98), fontsize=FONTSIZE, frameon=False)

#Additional Labels
zeroAx, oneAx = ax[0][1].twinx(), ax[1][1].twinx()
#zeroAx.set_ylabel(r"$x=0.000$", fontsize=FONTSIZE)
#oneAx.set_ylabel(r"$x=1.000$", fontsize=FONTSIZE)
zeroAx.yaxis.set_ticks([])
oneAx.yaxis.set_ticks([])

#ax[0][0].set_title(r"$\mathrm{t\!-\!phase}$")
#ax[0][1].set_title(r"$\mathrm{o\!-\!phase}$")
xlabel=r"$\mathrm{E-E}_\mathrm{F}$ $\mathrm{(eV)}$"
fig.text(0.5, 0.025, xlabel, fontsize=FONTSIZE, transform=fig.transFigure,
         verticalalignment="center", horizontalalignment="center")
ylabel=r"$\mathrm{Electronic\ LDOS\ (states/eV)}$"
fig.text(0.025, 0.5, ylabel, rotation=90, fontsize=FONTSIZE, transform=fig.transFigure,
         verticalalignment="center", horizontalalignment="center")

plt.savefig(PDOS_SAVE_LOC + ".pdf")
plt.savefig(PDOS_SAVE_LOC + ".png", dpi=1000)
fig.clf()
