from matplotlib import rcParams
from matplotlib import pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import headerXrd as b

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

tetraFolder = "t-tXrd"
orthoFolder = "o-oXrd"

XRD_DATA_LOC = {"tetra": tetraFolder + "//" + "xrdx",
               "ortho": orthoFolder + "//" + "xrdx"}
XRD_GRAPH_SAVE_LOC = "allXrd"

P1_P2_CUTOFF = 0.5
#color, linewidth, linestyle, label
LINE_DETAILS = {"tetra":   ["#ea0000", 0.8, '-', r"$\mathrm{t\!-\!Cu_2CdGe(S_xSe_{1\!-\!x})_4}$"], #red
                "ortho":   ["#0080FF", 0.8, '--', r"$\mathrm{o\!-\!Cu_2CdGe(S_xSe_{1\!-\!x})_4}$"]} #blue

def GetLinestyle(phase, con):
    return LINE_DETAILS[phase]


#Create a figure, loop through each pattern
include = list(np.arange(0.000, 1.125, 0.125))
fig, ax = plt.subplots(len(include), 1, sharex=True, figsize=(WIDTH, HEIGHT))
legendHandles = {}
for phaseId in ["tetra", "ortho"]:
    include = list(np.arange(0.000, 1.125, 0.125))
    include = [(f'{i:.3f}'[0] + f'{i:.3f}'[2:]) for i in include] ##gets rid of the decimal point
    for num, conc in enumerate(include):
        #Parse the correct output file, get the lists of angles and intensities
        #Also make a list of significant peaks to label the miller indices
        xrdData = b.ParseOutfile(XRD_DATA_LOC[phaseId] + conc + ".txt")
        x, y, mill = [], [], []
        for d in xrdData:
            x.append(d.ang)
            y.append(d.intsty)
            if(d.intsty > 10):
                mill.append([d.ang, d.intsty, d.h, d.k, d.l])
        #Get rid of copied indices.  Dont include negative indices.
        millTmp = []
        for m in mill:
            if (m not in millTmp and m[2:] not in [m_[2:] for m_ in millTmp] and
                (m[2]>=0 and m[3]>=0 and m[4]>=0)):
                millTmp.append(m)
        mill = millTmp[:]

        #Fill in zero-intensities to give graph a sharper look
        dx = 0.2
        extras = []
        for i in range(0, len(x) - 1):
            if((x[i+1] - x[i]) > dx):
                tmp = np.arange(x[i], x[i+1], dx)
                for t in tmp:
                    extras.append(t)
        for e in extras:
            x.append(e)
            y.append(0.0)
        x, y = zip(*sorted(zip(x, y))) ##literally magic.  Sorting parallel arrays at the same time

        #Axes labels.  Why I cant use this method for the y axis is a mystery to me
        if(conc == "1000"):
            ax[num].set_xlabel(r"$2\theta$ $\mathrm{(degrees)}$", fontsize=FONTSIZE)

        #Plotting
        ##Place concentration on appropriate plot
        if(phaseId == "tetra"):
            ax[num].text(22, 86, r"$x=$" + f'{float(conc)/1000.:.3f}', fontsize=FONTSIZE)

        ##Remove ticks unless they're necessary.  Set font size
        if(num != len(include) - 1):
            ax[num].axes.xaxis.set_visible(False)
        ax[num].axes.yaxis.set_visible(False)
        for tick in ax[num].xaxis.get_major_ticks():
            tick.label.set_fontsize(FONTSIZE)
        ##Limits, plot points
        ax[num].set_xlim(20, 60)#(26.5, 28.5)
        ax[num].set_ylim(0, 110)
        style = GetLinestyle(phaseId, float(conc)/1000.)
        ax[num].plot(x, y, color=style[0], linewidth=style[1], linestyle=style[2])
        #add to legend dict
        if(style[3] not in legendHandles.keys() and style[3] != ''):
            legendHandles[style[3]] = style

#Main plot edits, save and finish
labs, lins = [], []
for key, val in legendHandles.items():
    labs.append(key)
    lins.append(Line2D([0], [0], color=val[0], linewidth=val[1], linestyle=val[2]))
fig.legend(lins, labs, loc="center", ncol=len(labs), bbox_to_anchor=(0.5, 0.98),
           fontsize=FONTSIZE, frameon=False)

plt.subplots_adjust(left=0.05, right=0.98, bottom=0.07, top=0.96, hspace=0.001)
ylabel="Intensity"
plt.text(18.5, 110*(len(include) - 1)/2. - len(ylabel), ylabel, rotation=90, fontsize=FONTSIZE)
plt.savefig(XRD_GRAPH_SAVE_LOC + ".pdf")
plt.savefig(XRD_GRAPH_SAVE_LOC + ".png", dpi=DPI)
fig.clf()
