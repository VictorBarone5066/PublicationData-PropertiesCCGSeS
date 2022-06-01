from matplotlib import rcParams
from matplotlib import pyplot as plt
import matplotlib.patches as mpat
import numpy as np
from scipy.signal import savgol_filter
import headerCohp as h

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

tetraFolder = "t-tCohp"
orthoFolder = "o-oCohp"
COHPCAR_LOC = {"tetra": tetraFolder + "//" + "COHPCARx",
               "ortho": orthoFolder + "//" + "COHPCARx"}
SAVE_LOC = "allCohp"

SAVGOL_WINSIZE = 31
SAVGOL_POLY = 3
RED, BLUE, GREEN = "#ea0000", "#0080FF", "#0e9400"
#dict[[elem1, elem2]] = [colors, hatching, alpha, fillBool, linestyle]
#Amazingly, lists can't be used as hash keys but tuples can.  Why?  I don't know.
PAIRS_OF_INTEREST = {tuple(["Se", "Cu"]): ["#ea0000", "", 1., False, "-"],      #red
                     tuple(["Se", "Cd"]): ["#0e9400", "", 1., False, "--"],      #green
                     tuple(["Se", "Ge"]): ["#00008b", "", 1., False, ":"],      #blue
                     tuple(['S',  "Cu"]): ["#ea0000", "", 1., False, "-"],     #red
                     tuple(['S',  "Cd"]): ["#0e9400", "", 1., False, "--"],     #green
                     tuple(['S',  "Ge"]): ["#00008b", "", 1., False, ":"]}     #blue
ALPHABET = "ab.......cd..............." #no longer the alphabet because I needed to change stuff and I'm lazy

def IndexToIJ(conc, phasId):
    if(phasId[0] == 't'):
        if(conc[0] == '0'):
            return 0, 0
        return 1, 0
    if(conc[0] == '0'):
        return 0, 1
    return 1, 1

def bs(s):
    if (s == "Cu"):
        return r"$\mathrm{Cu}$"
    if (s == "Cd"):
        return r"$\mathrm{Cd}$"
    if (s == "Ge"):
        return r"$\mathrm{Ge}$"

usedKeys = [] ##Now, not necessairiy every potential atom pair is actually in the files
#Create a figure, loop through each pattern
fig, ax = plt.subplots(2, 2, sharex=True, figsize=(WIDTH, HEIGHT))
for phaseNum, phaseId in enumerate(["tetra", "ortho"]):
    include = [0.0, 1.0]
    include = [(f'{i:.3f}'[0] + f'{i:.3f}'[2:]) for i in include] ##gets rid of the decimal point
    for num, conc in enumerate(include):
        i, j = IndexToIJ(conc, phaseId)

        #Get the sums of the PAIRS_OF_INTEREST pairs seperatly
        data = h.ParseCohpFile(COHPCAR_LOC[phaseId] + conc)
        pairs = {} ##pairs[[elem1, elem2]] = [[energies], [pcohps], [ipcohps]]
        for d in data:
            thisKey = None
            for poi in PAIRS_OF_INTEREST.keys():
                if(set([d.elem1, d.elem2]) == set(poi)):
                    thisKey = poi[:]
            if(thisKey != None):
                if(thisKey not in pairs.keys()):
                    pairs[thisKey] = [d.energies[:], d.pCohps[:], d.ipCohps[:]]
                else:
                    for m in range(0, len(pairs[thisKey][1])):
                        pairs[thisKey][1][m] += d.pCohps[m]

        #Turn pCOHP -> -pCohp.  Ditto for ipCohps
        for vals in pairs.values():
            vals[1] = [-v for v in vals[1]]
            vals[2] = [-v for v in vals[2]]

        #Individual plot customizations
        ax[i][j].set_xlim(-6, 6)
        ax[i][j].set_ylim(-36, 36)
        ax[i][j].set_xticks(np.arange(-6, 8, 2))
        ax[i][j].set_yticks(np.arange(-36, 48, 12))

        ##Place concentration and alphabet number on appropriate plot
        ax[i][j].text(-5.65, 31.45, ALPHABET[num + 9*phaseNum] + ')', fontsize=FONTSIZE)


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


        #Plot COHPs
        ##Individual pair COHPs
        for key, vals in pairs.items():
            thisLabel = key[1]
            thisColor = PAIRS_OF_INTEREST[key][0]
            thisHatch = PAIRS_OF_INTEREST[key][1]
            thisAlpha = PAIRS_OF_INTEREST[key][2]
            thisFillB = PAIRS_OF_INTEREST[key][3]
            thisLinestyle = PAIRS_OF_INTEREST[key][4]

            if key not in usedKeys:
                usedKeys.append(key)
            ###Splits each piece into only positive and negative parts for reasons.  Yikes
            x, y = vals[0][:], list(savgol_filter(vals[1][:], SAVGOL_WINSIZE, SAVGOL_POLY))
            #y = GetRootScaling(y)
            x.append(x[-1] + 0.01)  #Adds a small endpoint to my functions
            y.append(-y[-1]/max(y)) #so that this algorithm works

            pieces, thisPiece = [], [[], []]
            for m in range(1, len(y)):
                if(y[m-1] * y[m] <= 0.0 and abs(y[m-1]) + abs(y[m]) != 0.0):
                    pieces.append(thisPiece)
                    thisPiece = [[x[m]], [0.0]]
                else:
                    thisPiece[0].append(x[m-1])
                    thisPiece[1].append(y[m-1])

            ###Plot each piece as if it were one continuous function
            for piece in pieces:
                #Implemented piece algrthm before adding filter.  Caused problems, this "fixes" them
                if(piece[0] == [] or piece[1] == []):
                    continue
                #Don't graph the unphysical tiny spikes caused by the filter's poly interp
                if(piece[0][-1] - piece[0][0] < 0.2):
                    continue
                piece[0].append(piece[0][-1] + 0.0001) #Adds an endpoint to each piece so that
                piece[1].append(0.0)                   #my pieces always end on the x axis

                #Get rid of insignificant contributions to DOS
                loCut, hiCut = -0.15, 0.15
                for ind in range(0, len(piece[1])):
                    if(piece[1][ind] < 0.0):
                        if(piece[1][ind] > loCut):
                            piece[1][ind] = 0
                        else:
                            piece[1][ind] -= loCut
                    elif(piece[1][ind] > 0.0):
                        if(piece[1][ind] < hiCut):
                            piece[1][ind] = 0
                        else:
                            piece[1][ind] -= hiCut

                ax[i][j].plot(piece[0], piece[1], color=thisColor, linestyle=thisLinestyle, linewidth=1.85)


        ax[i][j].axvline(0.0, color='k', linewidth=1.5, linestyle='--')
        ax[i][j].axhline(0.0, color='k', linewidth=1.0)
        ax[i][j].grid(True, alpha=0.5)

#Main plot edits, save and finish
plt.subplots_adjust(left=0.08, right=0.95, bottom=0.07, top=0.985, hspace=0.09, wspace=0.1)
##Legend stuff
for key, val in PAIRS_OF_INTEREST.items():
    for i in range(0, 2):
        for j in range(0, 2):
            if(key[0] == "Se"):
                continue
            if(i == 0):
                ax[i][j].plot(-200, -200, color=val[0], alpha=val[2], linestyle=val[4],
                              label=bs(key[1]) + r"$\mathrm{-Se}$")
            if(i == 1):
                ax[i][j].plot(-200, -200, color=val[0], alpha=val[2], linestyle=val[4],
                              label=bs(key[1]) + r"$\mathrm{-S}$")

ax[0][0].legend(title=r"$\mathrm{t\!-\!Cu_2CdGeSe_4}$", title_fontsize=FONTSIZE-2, loc="upper right", fontsize=FONTSIZE-2)
ax[1][0].legend(title=r"$\mathrm{t\!-\!Cu_2CdGeS_4}$", title_fontsize=FONTSIZE-2, loc="upper right", fontsize=FONTSIZE-2)
ax[0][1].legend(title=r"$\mathrm{o\!-\!Cu_2CdGeSe_4}$", title_fontsize=FONTSIZE-2, loc="upper right", fontsize=FONTSIZE-2)
ax[1][1].legend(title=r"$\mathrm{o\!-\!Cu_2CdGeS_4}$", title_fontsize=FONTSIZE-2, loc="upper right", fontsize=FONTSIZE-2)

#Additional Labels
zeroAx, oneAx = ax[0][1].twinx(), ax[1][1].twinx()
#zeroAx.set_ylabel(r"$x=0.000$", fontsize=FONTSIZE, fontweight="bold")
#oneAx.set_ylabel(r"$x=1.000$", fontsize=FONTSIZE)
zeroAx.yaxis.set_ticks([])
oneAx.yaxis.set_ticks([])

#ax[0][0].set_title("$\mathrm{t\!-\!Cu_2CdGe(S_xSe_{1\!-\!x})_4}$")
#ax[0][1].set_title("$\mathrm{o\!-\!Cu_2CdGe(S_xSe_{1\!-\!x})_4}$")
xlabel=r"$\mathrm{E-E}_\mathrm{F}$ $\mathrm{(eV)}$"
fig.text(0.5, 0.025, xlabel, fontsize=FONTSIZE, transform=fig.transFigure,
         verticalalignment="center", horizontalalignment="center")
ylabel=r"$\mathrm{-pCOHP}$"
fig.text(0.025, 0.5, ylabel, rotation=90, fontsize=FONTSIZE, transform=fig.transFigure,
         verticalalignment="center", horizontalalignment="center")

plt.savefig(SAVE_LOC + ".pdf")
plt.savefig(SAVE_LOC + ".png", dpi=DPI)
fig.clf()
