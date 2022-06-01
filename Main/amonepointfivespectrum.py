from matplotlib import pyplot as plt
import numpy as np

#Figure Options
plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams["mathtext.fontset"] = "dejavuserif"
plt.tight_layout = True
plt.rcParams["figure.autolayout"] = False

WIDTH, HEIGHT = 10.5, 4.5 ##inches.  3.5 for 1 col, 7.2 for 2 col
FONTSIZE = 12
DPI = 1500
fig, ax = plt.subplots(1, 1, figsize=(WIDTH, HEIGHT))

x, y = [], []
with open("am1.5spectrum.csv", 'r') as infile:
    for n, line in enumerate(infile):
        if(n < 2):
            continue

        x.append(1239.841875/float(line.split()[0])) ##magic number = inv nm to eV
        y.append(float(line.split()[3]))
    infile.close()


#Plot Configuration
ax.plot(x, y, 'k')
ax.fill(x, y, color='gray', alpha=0.5)

ax.axvline(x=1.8, color='r')
ax.axvline(x=2.8, color='r')

ax.set_xlim(1., 3.5)
ax.set_ylim(0.0, 1.5)

ax.grid(True, alpha=0.5)

for tick in ax.xaxis.get_major_ticks():
    tick.label.set_fontsize(FONTSIZE)
for tick in ax.yaxis.get_major_ticks():
    tick.label.set_fontsize(FONTSIZE)
#ax.set_xticks([0.0, 0.250, 0.50, 0.750, 1.0])

plt.subplots_adjust(left=0.07, right=0.98, bottom=0.14, top=0.97)

ax.set_xlabel(r"Energy (eV)", fontsize=FONTSIZE)
ax.set_ylabel(r"Spectral Irradiance (W$\,$m$^{-2}$$\,$nm$^{-1}$)", fontsize=FONTSIZE)

#plt.show()
plt.savefig("sir" + ".pdf")
plt.savefig("sir" + ".png", dpi=DPI)
