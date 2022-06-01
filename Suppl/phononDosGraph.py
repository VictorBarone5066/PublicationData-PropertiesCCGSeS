from matplotlib import pyplot as plt
from matplotlib import rcParams
import numpy as np
from matplotlib.lines import Line2D

#Figure Options
rcParams['axes.linewidth'] = 1.5
plt.rcParams["font.weight"] = "bold"
plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams["mathtext.fontset"] = "dejavuserif"
plt.tight_layout = True
plt.rcParams["figure.autolayout"] = False

WIDTH, HEIGHT = 3.5, 3.5 ##inches.  3.5 for 1 col, 7.2 for 2 col
FONTSIZE = 9
DPI = 1000

tetraFolder = "t-tPhononDos"
orthoFolder = "o-oPhononDos"

OUTCAR_LOC = {"tetra": tetraFolder + "//" + "dosx",
              "ortho": orthoFolder + "//" + "dosx"}
SAVE_LOC = "phDos"

WHICH = "ortho"

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

fig, ax = plt.subplots(1, 1, sharex=True, figsize=(WIDTH, HEIGHT))
include = list(np.arange(0.000, 1.125, 0.125))
include = [(f'{i:.3f}'[0] + f'{i:.3f}'[2:]) for i in include]  ##gets rid of the decimal point

xMin = 987654321
for conc in include:
    x, y = [], []  ##frequency (THz), DOS

    ##Read infile
    with open(OUTCAR_LOC[WHICH] + conc, 'r') as infile:
        for n, lin in enumerate(infile):
            if(n == 0):
                continue
            line = lin.split()
            x.append(float(line[0]) - 0.25)
            y.append(float(line[1]) - -0.5*0.125)

            if(x[-1] < xMin):
                xMin = x[-1]

        infile.close()

    ax.plot(x, y, color=linestyles[conc][0], linestyle=linestyles[conc][1],
                  marker=linestyles[conc][2], markevery=linestyles[conc][3],
                  markersize=linestyles[conc][4])
    ax.grid(True, alpha=0.5)

ax.set_xlabel("Frequency (THz)", fontsize=FONTSIZE, weight="bold")
ax.set_ylabel("DOS (States/THz)", fontsize=FONTSIZE, weight="bold")

labs, lins = [], []
for key, val in linestyles.items():
    labs.append(r"$x$ = " + f'{float(key) / 1000.:.3f}')
    lins.append(Line2D([0], [0], color=val[0], linestyle=val[1], marker=val[2]))
fig.legend(lins, labs, loc="center", bbox_to_anchor=(0.5, 0.905), fontsize=FONTSIZE, frameon=False,
          ncol=3)

plt.xlim(min(-2.2, xMin), 4)
plt.ylim(-0.01, 25)

plt.subplots_adjust(left=0.12, right=0.98, bottom=0.11, top=0.81, hspace=0.08, wspace=0.1)
plt.savefig(WHICH + SAVE_LOC + ".pdf")
plt.savefig(WHICH + SAVE_LOC + ".png", dpi=DPI)
