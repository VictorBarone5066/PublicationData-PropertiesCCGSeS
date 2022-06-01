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

WIDTH, HEIGHT = 7.2, 8.0 ##inches.  3.5 for 1 col, 7.2 for 2 col
FONTSIZE = 12
DPI = 1000

tetraFolder = "t-tCohp"
orthoFolder = "o-oCohp"

COHPCAR_LOC = {"tetra": tetraFolder + "//" + "COHPCARx",
               "ortho": orthoFolder + "//" + "COHPCARx"}
SAVE_LOC = "cohp"

SAVGOL_WINSIZE = 31
SAVGOL_POLY = 3

#dict[[elem1, elem2]] = [colors, hatching, alpha, fillBool]
#Amazingly, lists can't be used as hash keys but tuples can.  Why?  I don't know.
PAIRS_OF_INTEREST_SE = {tuple(["Se", "Cu"]): ["#0e9400", "",     0.4,  True, "-"], #d pink
                     tuple(["Se", "Cd"]):    ["#ea0000", "",     0.4,  True, "--"],    #d green
                     tuple(["Se", "Ge"]):    ["#00008b", "",     0.4,  True, ":"]}    #l blue
PAIRS_OF_INTEREST_S = {tuple(['S',  "Cu"]):  ["#0e9400", "xxxx", 0.65, False, "-"], #not l pink
                     tuple(['S',  "Cd"]):    ["#ea0000", "||||", 0.65, False, "--"],   #not l green
                     tuple(['S',  "Ge"]):    ["#00008b", "----", 0.65, False, ":"]}   #not l blue
ALPHABET = "abcdefghijklmnopqrstuvwxyz"

def GetRootScaling(f):
    f_ = f[:]
    for m in range(0, len(f)):
        if(f_[m] < 0.0):
            f_[m] = -(abs(f_[m]))**(1./2.)
        else:
            f_[m] = (f_[m])**(1./2.)
    return f #lol don't return root scaling.  this is a super good idea :)

##pairs[[elem1, elem2]] = [[energies], [pcohps]]
def GetPairs(dat, el):
    pai = {}
    if(el == "Se"):
        for d in dat:
            thisKey = None
            for poi in PAIRS_OF_INTEREST_SE.keys():
                if (set([d.elem1, d.elem2]) == set(poi)):
                    thisKey = poi[:]
            if (thisKey != None):
                if (thisKey not in pai.keys()):
                    pai[thisKey] = [d.energies[:], d.pCohps[:], d.ipCohps[:]]
                else:
                    for m in range(0, len(pai[thisKey][1])):
                        pai[thisKey][1][m] += d.pCohps[m]
                        pai[thisKey][2][m] += d.ipCohps[m]
    if(el == "S"):
        for d in dat:
            thisKey = None
            for poi in PAIRS_OF_INTEREST_S.keys():
                if (set([d.elem1, d.elem2]) == set(poi)):
                    thisKey = poi[:]
            if (thisKey != None):
                if (thisKey not in pai.keys()):
                    pai[thisKey] = [d.energies[:], d.pCohps[:], d.ipCohps[:]]
                else:
                    for m in range(0, len(pai[thisKey][1])):
                        pai[thisKey][1][m] += d.pCohps[m]
                        pai[thisKey][2][m] += d.ipCohps[m]
    return pai

def GetPieces(x, y):
    x.append(x[-1] + 0.01)  # Adds a small endpoint to my functions
    y.append(-y[-1] / max(y))  # so that this algorithm works

    pieces, thisPiece = [], [[], []]
    for m in range(1, len(y)):
        if (y[m - 1] * y[m] <= 0.0 and abs(y[m - 1]) + abs(y[m]) != 0.0):
            pieces.append(thisPiece)
            thisPiece = [[x[m]], [0.0]]
        else:
            thisPiece[0].append(x[m - 1])
            thisPiece[1].append(y[m - 1])

    ###Plot each piece as if it were one continuous function
    for piece in pieces:
        # Implemented piece algrthm before adding filter.  Caused problems, this "fixes" them
        if (piece[0] == [] or piece[1] == []):
            continue
        # Don't graph the unphysical tiny spikes caused by the filter's poly interp
        if (piece[0][-1] - piece[0][0] < 0.2):
            continue
        piece[0].append(piece[0][-1] + 0.0001)  # Adds an endpoint to each piece so that
        piece[1].append(0.0)  # my pieces always end on the x axis
    return pieces

def Cut(piece):
    # Get rid of insignificant contributions to DOS (+/- 0.4 seem significant, but
    # recall that we're now using root scaling.  This matches with the DOS b.g. cutoff
    loCut, hiCut = -0.16, 0.16
    for ind in range(0, len(piece[1])):
        if (piece[1][ind] < 0.0):
            if (piece[1][ind] > loCut):
                piece[1][ind] = 0
            else:
                piece[1][ind] -= loCut
        elif (piece[1][ind] > 0.0):
            if (piece[1][ind] < hiCut):
                piece[1][ind] = 0
            else:
                piece[1][ind] -= hiCut
    return

def MySuperCoolFunction(a, f, n):
    # Individual plot customizations
    a[n].set_xlim(-6, 6)
    a[n].set_ylim(-36, 36)
    a[n].set_xticks(np.arange(-6, 8, 2))
    a[n].set_yticks(np.arange(-36, 54, 18))

    ##Place concentration and alphabet number on appropriate plot
    a[n].text(4.0, 18.0, r"$x=$" + f'{float(conc) / 1000:.3f}', fontsize=FONTSIZE)

    ##Tick settings
    for tick in a[n].xaxis.get_major_ticks():
        tick.label.set_fontsize(FONTSIZE)
    for tick in a[n].yaxis.get_major_ticks():
        tick.label.set_fontsize(FONTSIZE)

    return


usedKeys = [] ##Now, not necessairiy every potential atom pair is actually in the files
#Create a figure, loop through each pattern
for phaseNum, phaseId in enumerate(["tetra", "ortho"]):
    """
    ********
    SE PAIRS
    ********
    """
    fig, ax = plt.subplots(8, 1, sharex=True, figsize=(WIDTH, HEIGHT))
    include = list(np.arange(0.000, 1.000, 0.125))
    include = [(f'{i:.3f}'[0] + f'{i:.3f}'[2:]) for i in include] ##gets rid of the decimal point
    for num, conc in enumerate(include):
        #Get the sums of the PAIRS_OF_INTEREST pairs seperatly
        data = h.ParseCohpFile(COHPCAR_LOC[phaseId] + conc)
        pairs = GetPairs(data, "Se")

        MySuperCoolFunction(ax, fig, num)

        #Turn pCOHP -> -pCohp.  Ditto for ipCohps
        for vals in pairs.values():
            vals[1] = [-v for v in vals[1]]
            vals[2] = [-v for v in vals[2]]

        #Plot COHPs
        ##Individual pair COHPs
        for key, vals in pairs.items():
            col = PAIRS_OF_INTEREST_SE[key][0]
            ls = PAIRS_OF_INTEREST_SE[key][4]

            if key not in usedKeys:
                usedKeys.append(key)
            ###Splits each piece into only positive and negative parts for reasons.  Yikes
            x, y = vals[0][:], list(savgol_filter(vals[1][:], SAVGOL_WINSIZE, SAVGOL_POLY))
            y = GetRootScaling(y)

            pieces = GetPieces(x, y)
            for num_, piece in enumerate(pieces):
                Cut(piece)
                if(num == 1 and num_ == 0):
                    ax[num].plot(piece[0], piece[1], label=key[0] + '-' + key[1], color=col, linestyle=ls,
                                 linewidth=1.85)
                else:
                    ax[num].plot(piece[0], piece[1], color=col, linestyle=ls, linewidth=1.85)


        ##The averaged out ipCOHPs - not scaled
        if(num == 1):
            ax[num].plot(data[0].energies, [-p*(len(data)-1) for p in data[0].pCohps], color='k',
                          linewidth=1.0, label=r"$\mathrm{Total\ pCOHP}$")
        else:
            ax[num].plot(data[0].energies, [-p*(len(data)-1) for p in data[0].pCohps], color='k',
                          linewidth=1.0)

        ax[num].axvline(0.0, color='k', linestyle="--", linewidth=1.0)
        ax[num].axhline(0.0, color='k', linewidth=1.0)
        ax[num].grid(True, alpha=0.5)

    #Main plot edits, save and finish
    for i in range(0, len(ax)):
        ###Really dumb but necessary workarounds to force gridlines but not tick markers
        if (num != len(ax) - 1):
            a[i].tick_params(axis='x', colors="white")
            for tick in a[n].xaxis.get_major_ticks():
                tick.label.set_fontsize(0.01)
    ax[-1].set_xlabel(r"$E-E_\mathrm{F}$ $\mathrm{(eV)}$", fontsize=FONTSIZE)
    ylabel=r"$\mathrm{-pCOHP}$"
    fig.text(0.025, 0.5, ylabel, rotation=90, fontsize=FONTSIZE, transform=fig.transFigure,
             verticalalignment="center", horizontalalignment="center")
    plt.subplots_adjust(left=0.08, right=0.98, bottom=0.07, top=0.96, hspace=0.23, wspace=0.0)
    fig.legend(loc="center", ncol=len(PAIRS_OF_INTEREST_SE) + 1, bbox_to_anchor=(0.5, 0.98),
               fontsize=FONTSIZE, frameon=False)

    plt.savefig(SAVE_LOC + phaseId + "se.pdf")
    plt.savefig(SAVE_LOC + phaseId + "se.png", dpi=DPI)
    fig.clf()

    """
    ********
    S PAIRS
    ********
    """
    fig, ax = plt.subplots(8, 1, sharex=True, figsize=(WIDTH, HEIGHT))
    include = list(np.arange(0.125, 1.125, 0.125))
    include = [(f'{i:.3f}'[0] + f'{i:.3f}'[2:]) for i in include] ##gets rid of the decimal point
    for num, conc in enumerate(include):
        #Get the sums of the PAIRS_OF_INTEREST pairs seperatly
        data = h.ParseCohpFile(COHPCAR_LOC[phaseId] + conc)
        pairs = GetPairs(data, "S")

        MySuperCoolFunction(ax, fig, num)

        #Turn pCOHP -> -pCohp.  Ditto for ipCohps
        for vals in pairs.values():
            vals[1] = [-v for v in vals[1]]
            vals[2] = [-v for v in vals[2]]

        #Plot COHPs
        ##Individual pair COHPs
        for key, vals in pairs.items():
            col = PAIRS_OF_INTEREST_S[key][0]
            ls = PAIRS_OF_INTEREST_S[key][4]

            if key not in usedKeys:
                usedKeys.append(key)
            ###Splits each piece into only positive and negative parts for reasons.  Yikes
            x, y = vals[0][:], list(savgol_filter(vals[1][:], SAVGOL_WINSIZE, SAVGOL_POLY))
            y = GetRootScaling(y)

            pieces = GetPieces(x, y)
            for num_, piece in enumerate(pieces):
                Cut(piece)
                if(num == 1 and num_ == 0):
                    ax[num].plot(piece[0], piece[1], label=key[0] + '-' + key[1], color=col, linestyle=ls,
                                 linewidth=1.85)
                else:
                    ax[num].plot(piece[0], piece[1], color=col, linestyle=ls, linewidth=1.85)


        ##The averaged out ipCOHPs - not scaled
        if(num == 1):
            ax[num].plot(data[0].energies, [-p*(len(data)-1) for p in data[0].pCohps], color='k',
                          linewidth=1.0, label=r"$\mathrm{Total\ pCOHP}$")
        else:
            ax[num].plot(data[0].energies, [-p*(len(data)-1) for p in data[0].pCohps], color='k',
                          linewidth=1.0)

        ax[num].axvline(0.0, color='k', linestyle="--", linewidth=1.0)
        ax[num].axhline(0.0, color='k', linewidth=1.0)
        ax[num].grid(True, alpha=0.5)

    #Main plot edits, save and finish
    for i in range(0, len(ax)):
        ###Really dumb but necessary workarounds to force gridlines but not tick markers
        if (num != len(ax) - 1):
            a[i].tick_params(axis='x', colors="white")
            for tick in a[n].xaxis.get_major_ticks():
                tick.label.set_fontsize(0.01)
    ax[-1].set_xlabel(r"$E-E_\mathrm{F}$ $\mathrm{(eV)}$", fontsize=FONTSIZE)
    ylabel=r"$\mathrm{-pCOHP}$"
    fig.text(0.025, 0.5, ylabel, rotation=90, fontsize=FONTSIZE, transform=fig.transFigure,
             verticalalignment="center", horizontalalignment="center")
    plt.subplots_adjust(left=0.08, right=0.98, bottom=0.07, top=0.96, hspace=0.23, wspace=0.0)
    fig.legend(loc="center", ncol=len(PAIRS_OF_INTEREST_S) + 1, bbox_to_anchor=(0.5, 0.98),
               fontsize=FONTSIZE, frameon=False)

    plt.savefig(SAVE_LOC + phaseId + "s.pdf")
    plt.savefig(SAVE_LOC + phaseId + "s.png", dpi=DPI)
    fig.clf()
