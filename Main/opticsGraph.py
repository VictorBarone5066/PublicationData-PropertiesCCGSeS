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

#Get solar spectrum
solX, solY = [], []
with open("am1.5spectrum.csv", 'r') as infile:
    for n, line in enumerate(infile):
        if(n < 2):
            continue

        solX.append(1239.841875/float(line.split()[0])) ##magic number = inv nm to eV
        solY.append(float(line.split()[3]))
    infile.close()



WIDTH, HEIGHT = 7.2, 7.2 ##inches.  3.5 for 1 col, 7.2 for 2 col
FONTSIZE = 12
DPI = 1000

tetraFolder = "t-tOptics"
orthoFolder = "o-oOptics"

OUTCAR_LOC = {"tetra": tetraFolder + "//" + "OUTCARx",
              "ortho": orthoFolder + "//" + "OUTCARx"}
OPTICS_SAVE_LOC = "allOptics"

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
linestyles = {"0000": ["#ff0000", "-",  None, None, 5],   #red
              "0250": ["#00b423", "-",  '|',  20,  10],   #lime
              "0500": ["black",   "--", '1',  20,  10],   #black
              "0750": ["#93ff27", ":",  '^',  20,  5],   #blue
              "1000": ["#cb34ff", ":",  'D',  20,  5]}   #pink
ALPHABET = "abcdefghijklmnopqrstuvwxyz"

#actually returns column, then row bc im dumb
def GetRowAndCol(typeStr, phas):
    if(typeStr[0] == 'a'): ##absorbance
        if(phas[0] == 't'): ###tetra
            return 0, 0
        return 1, 0 ###ortho
    if(typeStr[0] == 'r'): ##reflectivity
        if(phas[0] == 't'): ###tetra
            return 0, 1
        return 1, 1 ###ortho
    return 0.5, 0.5 ##will raise index error if this ever happens

#Create figure, loop through each pattern
fig, ax = plt.subplots(2, 2, sharex=True, figsize=(WIDTH, HEIGHT))
for phaseNum, phaseId in enumerate(["tetra", "ortho"]):
    include = [0.0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.0]
    include = [0.0, 0.25, 0.5, 0.75, 1.0]
    include = [(f'{i:.3f}'[0] + f'{i:.3f}'[2:]) for i in include] ##gets rid of the decimal point
    for num, conc in enumerate(include):
        #Rename to make stuff easier to work with
        dieLis = h.GetDielecTensors(OUTCAR_LOC[phaseId] + conc)
        for i in range(0, len(dieLis)):
            if(dieLis[i][0].thry == "c-c"):
                if(dieLis[i][0].imag):
                    chgIm = dieLis[i][:]
                if(dieLis[i][0].real):
                    chgRe = dieLis[i][:]
            if(dieLis[i][0].thry == "d-d"):
                if(dieLis[i][0].imag):
                    denIm = dieLis[i][:]
                if(dieLis[i][0].real):
                    denRe = dieLis[i][:]

        #Absorbance Graphs
        i, j = GetRowAndCol("abs", phaseId)

        ##Graph formatting
        ax[j][i].set_xlim(0.5, 5)
        ax[j][i].set_ylim(0, 80) ##for absorption
        axT = ax[j][i].twinx() ##for solar spectrum
        axT.set_ylim(0, 1.5)

        for tick in ax[j][i].xaxis.get_major_ticks():
            tick.label.set_fontsize(FONTSIZE)
        for tick in ax[j][i].yaxis.get_major_ticks():
            tick.label.set_fontsize(FONTSIZE)
        for tick in axT.yaxis.get_major_ticks():
            tick.label.set_fontsize(FONTSIZE+8)

        if(j != 1):
            ax[j][i].tick_params(axis='x', colors="white")
            for tick in ax[j][i].xaxis.get_major_ticks():
                tick.label.set_fontsize(0.01)
        if(i != 0):
            ax[j][i].tick_params(axis='y', colors="white")
            for tick in ax[j][i].yaxis.get_major_ticks():
                tick.label.set_fontsize(0.01)
        if(i != 1):
            axT.tick_params(axis='y', colors="white")
            for tick in axT.yaxis.get_major_ticks():
                tick.label.set_fontsize(0.01)
        if(j == 0 and i == 0):
            ax[j][i].set_title(r"$\mathrm{t\!-\!Cu_2CdGe(S_xSe_{1\!-\!x})_4}$")
        if(j == 0 and i == 1):
            ax[j][i].set_title(r"$\mathrm{o\!-\!Cu_2CdGe(S_xSe_{1\!-\!x})_4}$")

        if(conc == "1000"):
            ax[j][i].text(0.5 + 0.05, 75.5, ALPHABET[2*phaseNum] + ')', fontsize=FONTSIZE)
        if(j == 0 and i == 0):
            ax[j][i].set_ylabel(r"$\mathrm{Absorption\ Coefficient\ (μm}^{-1})$", fontsize=FONTSIZE)
        if(j == 0 and i == 1):
            axT.set_ylabel(r"Spectral Irradiance (W$\,$m$^{-2}$$\,$nm$^{-1}$)", fontsize=FONTSIZE+2,
                           fontweight='bold', fontfamily="Times New Roman")

        ##Plot absorption
        x = [x_.energy for x_ in chgRe]
        y = []
        for rea, img in zip(chgRe, chgIm):
            y.append(h.AbsCoeff(rea, img)[0]*10**(-6)) #now in microm^-1
        x, y = h.GetSplineInterp(x, y)
        if(conc == '0000'):
            axT.plot(solX, solY, color='k', alpha=0.5, linewidth=0.5)
            axT.fill(solX, solY, color='grey', alpha=0.40)
            if(linestyles[conc][-1] == None):
                ax[j][i].plot(x, y, label=r"$x$ = " + f'{float(conc)/1000:.3f}', color=linestyles[conc][0],
                              linestyle=linestyles[conc][1], linewidth=1.85)
            else:
                ax[j][i].plot(x, y, label=r"$x$ = " + f'{float(conc)/1000:.3f}', color=linestyles[conc][0],
                              linestyle=linestyles[conc][1], marker=linestyles[conc][2],
                              markevery=linestyles[conc][3], markersize=linestyles[conc][4], linewidth=1.85)

        else:
            if(linestyles[conc][-1] == None):
                ax[j][i].plot(x, y, color=linestyles[conc][0],
                              linestyle=linestyles[conc][1], linewidth=1.85)
            else:
                ax[j][i].plot(x, y, color=linestyles[conc][0],
                              linestyle=linestyles[conc][1], marker=linestyles[conc][2],
                              markevery=linestyles[conc][3], markersize=linestyles[conc][4], linewidth=1.85)


        ax[j][i].grid(True, alpha=0.5)

        #Reflect Graphs
        i, j = GetRowAndCol("ref", phaseId)

        ##Graph formatting
        ax[j][i].set_xlim(0.5, 5)
        ax[j][i].set_ylim(0.19, 0.33) ##for reflectivity
        axT = ax[j][i].twinx() ##for solar spectrum
        axT.set_ylim(0, 1.5)

        for tick in ax[j][i].xaxis.get_major_ticks():
            tick.label.set_fontsize(FONTSIZE)
        for tick in ax[j][i].yaxis.get_major_ticks():
            tick.label.set_fontsize(FONTSIZE)
        for tick in axT.yaxis.get_major_ticks():
            tick.label.set_fontsize(FONTSIZE)

        if(j != 1):
            ax[j][i].tick_params(axis='x', colors="white")
            for tick in ax[j][i].xaxis.get_major_ticks():
                tick.label.set_fontsize(0.01)
        if(i != 0):
            ax[j][i].tick_params(axis='y', colors="white")
            for tick in ax[j][i].yaxis.get_major_ticks():
                tick.label.set_fontsize(0.01)
        if(i != 1):
            axT.tick_params(axis='y', colors="white")
            for tick in axT.yaxis.get_major_ticks():
                tick.label.set_fontsize(0.01)

        if(conc == "1000"):
            ax[j][i].text(0.5 + 0.05, 0.3217, ALPHABET[1 + 2*phaseNum] + ')', fontsize=FONTSIZE)
        if(j == 1 and i == 0):
            ax[j][i].set_ylabel(r"$\mathrm{Reflectivity}$", fontsize=FONTSIZE)
        if(j == 1 and i == 1):
            axT.set_ylabel(r"Spectral Irradiance (W$\,$m$^{-2}$$\,$nm$^{-1}$)", fontsize=FONTSIZE+2,
                           fontweight='bold', fontfamily="Times New Roman")


        ##Plot reflectivity
        if(conc == "0000"):
            axT.plot(solX, solY, color='k', alpha=0.5, linewidth=0.5)
            axT.fill(solX, solY, color='grey', alpha=0.40)

        x = [x_.energy for x_ in chgRe]
        y = []
        for rea, img in zip(chgRe, chgIm):
            y.append(h.Reflec(rea, img)[0])
        x, y = h.GetSplineInterp(x, y)
        if(linestyles[conc][-1] == None):
            ax[j][i].plot(x, y, label=r"$x$ = " + f'{float(conc)/1000:.3f}', color=linestyles[conc][0],
                          linestyle=linestyles[conc][1], linewidth=1.85)
        else:
            ax[j][i].plot(x, y, label=r"$x$ = " + f'{float(conc)/1000:.3f}', color=linestyles[conc][0],
                          linestyle=linestyles[conc][1], marker=linestyles[conc][2],
                          markevery=linestyles[conc][3], markersize=linestyles[conc][4], linewidth=1.85)


        ax[j][i].grid(True, alpha=0.5)

#Full Plot Edits
plt.subplots_adjust(left=0.09, right=0.92, bottom=0.072, top=0.935, hspace=0.08, wspace=0.1)

xlabel=r"$\mathrm{Energy\ (eV)}$"
fig.text(0.5, 0.025, xlabel, fontsize=FONTSIZE, transform=fig.transFigure,
         verticalalignment="center", horizontalalignment="center")

#Legend Stuff
labs, lins = [], []
for key, val in linestyles.items():
    labs.append(r"$x$ = " + f'{float(key) / 1000.:.3f}')
    lins.append(Line2D([0], [0], color=val[0], linestyle=val[1]))
fig.legend(lins, labs, loc="center", bbox_to_anchor=(0.5, 0.985), fontsize=FONTSIZE, frameon=False,
           ncol=5)

##Stupid way to do legends:
#for i in range(0, 2):
#    for j in range(0, 2):
#        labs, lins = [], []
#        for key, val in linestyles.items():
#            labs.append(r"$x$ = " + f'{float(key)/1000.:.3f}')
#            lins.append(Line2D([0], [0], color=val[0], linestyle=val[1], marker=val[2]))
#        if(j == 0):
#            ax[i][j].legend(lins, labs, labelspacing=0.4,borderpad=0.4,fontsize=FONTSIZE-3, ncol=2,
#                            title=r"$\mathrm{t\!-\!Cu_2CdGe(S_xSe_{1\!-\!x})_4}$", title_fontsize=FONTSIZE-3)
#        if(j == 1):
#            ax[i][j].legend(lins, labs, labelspacing=0.4,borderpad=0.4,fontsize=FONTSIZE-3, ncol=2,
#                            title=r"$\mathrm{o\!-\!Cu_2CdGe(S_xSe_{1\!-\!x})_4}$", title_fontsize=FONTSIZE-3)


plt.savefig(OPTICS_SAVE_LOC + ".pdf")
plt.savefig(OPTICS_SAVE_LOC + ".png", dpi=DPI)
fig.clf()
