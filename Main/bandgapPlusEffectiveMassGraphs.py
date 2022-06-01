from matplotlib import rcParams
from matplotlib import pyplot as plt
import numpy as np
import headerDensityOfStates as dos

import headerPoscar as POH
import headerBandPlot as PBH
import headerAemt as EMH
from scipy.interpolate import CubicSpline

#Figure Options
rcParams['axes.linewidth'] = 1.5
plt.rcParams["font.weight"] = "bold"
plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams["mathtext.fontset"] = "dejavuserif"
plt.tight_layout = True
plt.rcParams["figure.autolayout"] = False

WIDTH, HEIGHT = 4.0, 5.2 ##inches.  3.5 for 1 col, 7.2 for 2 col
FONTSIZE = 12
DPI = 1000


tetraFolder = "t-tBandgap"
orthoFolder = "o-oBandgap"
DOSCAR_LOC =  {"tetra": "t-tBandgap" + "//" + "DOSCARx",
               "ortho": "o-oBandgap" + "//" + "DOSCARx"}
EIGENVAL_LOC = {"tetra": "t-tEffMass" + "//" + "EIGENVALx",
              "ortho": "o-oEffMass" + "//" + "EIGENVALx"}
POSCAR_LOC = {"tetra": "t-tEffMass" + "//" + "POSCARx",
               "ortho": "o-oEffMass" + "//" + "POSCARx"}
NEW_MASS_LOC = {"tetra": "t-tEffMassNew" + "//" + "outTetra-x",
                "ortho": "o-oEffMassNew" + "//" + "outOrtho-x"}
OUTFILE_LOC = "bgAndEffMassGraphs"

TOL = 0.100 #in whatever units the DOS is written as

MASS_TYPE = "cond" ##cond or doss

#color, linewidth, linestyle, label, fitOrder, markerstyle
BG_LINE_DETAILS = {"tetra":   ["#ea0000", 1.85, '-',  r"$\mathrm{t\!-\!Cu_2CdGe(S_xSe_{1\!-\!x})_4}$", 1, 'o'], #red
                "ortho":   ["#0080FF", 1.85, '--', r"$\mathrm{o\!-\!Cu_2CdGe(S_xSe_{1\!-\!x})_4}$", 1, '^']} #blue
#-------------------------------------color, linewidth, linestyle, label, interp poly order, marker
ME_LINE_DETAILS = {"tetra":   {"hole": ["#ea0000", 1.85, '-',  r"$\mathrm{t\!-\!Cu_2CdGe(S_xSe_{1\!-\!x})_4}$", 1, 'o'], #red
                            "elec": ["#ea0000", 1.85, '-',  r"$\mathrm{t\!-\!Cu_2CdGe(S_xSe_{1\!-\!x})_4}$", 1, 'o']  #red
                            },
                "ortho":   {"hole": ["#0080FF", 1.85, '--',  r"$\mathrm{o\!-\!Cu_2CdGe(S_xSe_{1\!-\!x})_4}$", 1, '^'], #blue
                            "elec": ["#0080FF", 1.85, '--',  r"$\mathrm{o\!-\!Cu_2CdGe(S_xSe_{1\!-\!x})_4}$", 1, '^']  #blue
                            }
                }


def GetLinestyle(phase, con):
    return BG_LINE_DETAILS[phase]

def FermiLoc(nrgs, ef):
    bst, loc = abs(nrgs[0] - ef), 0
    for i in range(0, len(nrgs)):
        if(abs(nrgs[i] - ef) < bst):
            bst = abs(nrgs[i] - ef)
            loc = i
    return loc

def PolyFit(x, y, deg = None, den = 0.01):
    if(deg == None):
        deg = len(x) - 1
    r = np.polyfit(x, y, deg)
    x_ = np.arange(x[0], x[-1], den)
    y_ = []
    for x__ in x_:
        y__ = 0
        for i in range(0, deg + 1):
            y__ += r[i]*x__**(deg - i)
        y_.append(y__)
    return list(x_),  y_

def GetAvgEffMass(lis):
    mEffAvg = 0.
    for i in range(2, len(lis)):
        mEffAvg += lis[i][0]

    return mEffAvg / (float(len(lis)) - 2.)

def GetSplineInterp(x, y, ptsForSpline=1000):
    spline = CubicSpline(x, y)
    xG = np.linspace(x[0], x[-1], ptsForSpline)
    yG = np.array(spline(xG))

    return list(xG), list(yG)

data = {} #data[phase] = [[concs], [bandgaps], [line details]]
#Create figure, loop through each pattern.  Initialize the data dict
fig, ax = plt.subplots(3, 1, sharex=True, figsize=(WIDTH, HEIGHT))
for phaseId in ["tetra", "ortho"]:
    include = list(np.arange(0.000, 1.125, 0.125))
    include = [(f'{i:.3f}'[0] + f'{i:.3f}'[2:]) for i in include] ##gets rid of the decimal point
    for num, conc in enumerate(include):

        #Find fermi level / its location, get energies and dos values in a more workable format
        eFermi = dos.GetEnergyInfo(DOSCAR_LOC[phaseId] + conc, 0)["eFermi"]
        vals = dos.GetAtomDosInfo(DOSCAR_LOC[phaseId] + conc, 0, spin = False)
        energies = vals["energy"]
        eFermiLoc = FermiLoc(energies, eFermi)
        doss = vals["dos"]

        #Verify that there are no states at eFermi
        if(doss[eFermiLoc] >= TOL):
            print("State(s) at the fermi level (" + str(eFermi) + ") eV!")
            exit

        #Get bandgap value
        ##Search backwards from eFermi
        i = eFermiLoc
        while(i > 0):
            if(doss[i] >= TOL):
                vbm = 1./2. * (energies[i + 1] + energies[i])
                break
            i = i - 1
        ##Serch forwards from eFermi
        i = eFermiLoc
        while(i < len(energies)):
            if(doss[i] >= TOL):
                cbm = 1./2. * (energies[i] + energies[i - 1])
                break
            i = i + 1

        #Add conc, bandgap to relevant dict
        style = GetLinestyle(phaseId, float(conc)/1000.)
        if(style[3] not in data.keys()):
            data[style[3]] = [[float(conc)/1000.], [cbm - vbm], style]
        else:
            data[style[3]][0].append(float(conc)/1000.) ##append x
            data[style[3]][1].append(cbm - vbm) ##append y
        print(phaseId, conc, np.round(cbm-vbm, decimals=3))

#For each seperate bit of data, plot data points and poly fits + add legend
for key, val in data.items():
    if(key != "else"):
        xP, yP = PolyFit(val[0], val[1], deg=val[2][4])
        ax[0].plot(xP, yP, color=val[2][0], linewidth=val[2][1], linestyle=val[2][2])
    ax[0].plot(val[0], val[1], color=val[2][0], marker=val[2][5], linewidth=0, markersize=6.5)

#Plot Configuration
ax[0].set_ylim(1.1, 1.9)
ax[0].grid(True, alpha=0.5)
for tick in ax[0].yaxis.get_major_ticks():
    tick.label.set_fontsize(FONTSIZE)
ax[0].set_ylabel(r"$E_\mathrm{g}$ $\mathrm{(eV)}$", fontsize=FONTSIZE)
ax[0].set_yticks([1.2, 1.5, 1.8])






massData = {}
for phaseId in ["tetra", "ortho"]:
    massData[phaseId] = {}
    massData[phaseId]["hole"] = {}
    massData[phaseId]["elec"] = {}

    include = list(np.arange(0.000, 1.125, 0.125))
    include = [(f'{i:.3f}'[0] + f'{i:.3f}'[2:]) for i in include] ##remove decim pt
    for conc in include:
        thisN, thisP = EMH.InfileToEigDat(infileLoc=NEW_MASS_LOC[phaseId] + conc + ".aem", n=True, p=True)
        nMu = list(thisN.keys())[0]
        pMu = list(thisP.keys())[0]
        ##hole
        massData[phaseId]["hole"][conc] = thisP[pMu][MASS_TYPE]["mass"][0]/1.1
        ##electron
        massData[phaseId]["elec"][conc] = thisN[nMu][MASS_TYPE]["mass"][0]/1.1

#Plot x and y
laTet, laOrt = False, False #To make labels not act obnoxious
for phaseKey, phaseDict in massData.items():
    for massKey, massDict in phaseDict.items():

        x, y = [], []
        for concKey, concList in massDict.items():
            x.append(float(concKey)/1000)
            y.append(concList)
        print(phaseKey, massKey, x, y)

        color = ME_LINE_DETAILS[phaseKey][massKey][0]
        width = ME_LINE_DETAILS[phaseKey][massKey][1]
        style = ME_LINE_DETAILS[phaseKey][massKey][2]
        label = ME_LINE_DETAILS[phaseKey][massKey][3]
        marke = ME_LINE_DETAILS[phaseKey][massKey][5]

        if(massKey == "hole"):
            xP, yP = GetSplineInterp(x, y)
            ax[1].plot(x, [y_ for y_ in y], color=color, marker=marke, linewidth=0., markersize=6.5)
            ax[1].plot(xP, [yp_ for yp_ in yP], color=color, linewidth=width, linestyle=style)
            if (laTet == False and phaseKey == "tetra"):
                laTet = True
                ax[1].plot(-200, -200, color=color, marker=marke, linewidth=width, markersize=6.5,
                           linestyle=style,
                           label=label)
        if(massKey == "elec"):
            xP, yP = GetSplineInterp(x, y)
            ax[2].plot(x, y, color=color, marker=marke, linewidth=0., markersize=6.5)
            ax[2].plot(xP, yP, color=color, linewidth=width, linestyle=style)
            if (laOrt == False and phaseKey == "ortho"):
                laOrt = True
                ax[2].plot(-200, -200, color=color, marker=marke, linewidth=width, markersize=6.5,
                           linestyle=style,
                           label=label)

#Final plot customizations
for i in range(0, 3):
    ax[i].set_xlim(-0.020, 1.020)
    ax[i].grid(True, alpha=0.5)

    ###Really dumb but necessary workarounds to force gridlines but not tick markers
    ax[0].tick_params(axis='x', colors="white")
    ax[1].tick_params(axis='x', colors="white")
    for tick in ax[i].xaxis.get_major_ticks():
        tick.label.set_fontsize(FONTSIZE)
    for tick in ax[i].yaxis.get_major_ticks():
        tick.label.set_fontsize(FONTSIZE)
    for tick in ax[0].xaxis.get_major_ticks():
        tick.label.set_fontsize(0.01)
    for tick in ax[1].xaxis.get_major_ticks():
        tick.label.set_fontsize(0.01)
    ax[i].set_xticks([0.0, 0.250, 0.50, 0.750, 1.0])

ax[1].set_yticks([1.9, 2.9, 3.9, 4.9, 5.9])
ax[2].set_yticks([1.0, 1.5, 2.0, 2.5, 3.0])
ax[1].set_ylim(1.0, 6.1)
ax[2].set_ylim(0.9, 3.1)

ax[2].set_xlabel(r"$x=[\mathrm{S}]\,/\,\left([\mathrm{S}]+[\mathrm{Se}]\right)$", fontsize=FONTSIZE)
ax[1].set_ylabel(r"$\langle m^*_\mathrm{h}$$\mathrm{/}$$m_\mathrm{0} \rangle$", fontsize=FONTSIZE)
ax[2].set_ylabel(r"$\langle m^*_\mathrm{e}$$\mathrm{/}$$m_\mathrm{0} \rangle$", fontsize=FONTSIZE)

fig.legend(loc="center", bbox_to_anchor=(0.5, 0.95), fontsize=FONTSIZE, frameon=False)

plt.subplots_adjust(left=0.17, right=0.965, bottom=0.1, top=0.910, hspace=0.05)
fig.savefig(OUTFILE_LOC + ".pdf")
fig.savefig(OUTFILE_LOC + ".png", dpi=DPI)
