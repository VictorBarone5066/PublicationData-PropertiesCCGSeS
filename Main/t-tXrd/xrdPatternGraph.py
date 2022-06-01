from matplotlib import pyplot as plt
import numpy as np
import xrdPatternParseBase as b

#Figure Options
plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams["mathtext.fontset"] = "dejavuserif"
plt.tight_layout = True
plt.rcParams["figure.autolayout"] = False

WIDTH, HEIGHT = 7.2, 7.2 ##inches.  3.5 for 1 col, 7.2 for 2 col
FONTSIZE = 12
DPI = 1000

GENERAL_LOC = "C:\\Users\\baron\\Desktop\\Cu2CdGe(SxSe1-x)4\\ortho-ortho\\data\\xrd\\"
CONTCAR_LOC = GENERAL_LOC + "..\\..\\convergedModels\\"
XRD_GRAPH_SAVE_LOC = GENERAL_LOC + "xrd"

#Create a figure, loop through each pattern
include = list(np.arange(0.000, 1.125, 0.125))
fig, ax = plt.subplots(len(include), 1, sharex=True, figsize=(WIDTH, HEIGHT))
include = [(f'{i:.3f}'[0] + f'{i:.3f}'[2:]) for i in include] ##gets rid of the decimal point
for num, conc in enumerate(include):
    #Parse the correct output file, get the lists of angles and intensities
    #Also make a list of significant peaks to label the miller indices
    xrdData = b.ParseOutfile(GENERAL_LOC + "xrdx" + conc + ".txt")
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
        ax[num].set_xlabel(r"$2\theta$ (degrees)", fontsize=FONTSIZE)

    #Plotting
    ##Place concentration on appropriate plot
    ax[num].text(22, 88, r"$x=$" + f'{float(conc)/1000:.3f}', fontsize=FONTSIZE)
    ##Remove ticks unless they're necessary.  Set font size
    if(num != len(include) - 1):
        ax[num].axes.xaxis.set_visible(False)
    ax[num].axes.yaxis.set_visible(False)
    for tick in ax[num].xaxis.get_major_ticks():
        tick.label.set_fontsize(FONTSIZE)
    ##Limits, plot points
    ax[num].set_xlim(20, 60)#(26.5, 28.5)
    ax[num].set_ylim(0, 110)
    ax[num].plot(x, y, 'k')

#Main plot edits, save and finish
plt.subplots_adjust(left=0.05, right=0.98, bottom=0.07, top=0.98, hspace=0.001)
ylabel="Intensity (arb.)"
plt.text(18.5, 110*(len(include) - 1)/2. - len(ylabel), ylabel, rotation=90, fontsize=FONTSIZE)
plt.savefig(XRD_GRAPH_SAVE_LOC + ".pdf")
plt.savefig(XRD_GRAPH_SAVE_LOC + ".png", dpi=1000)
fig.clf()







