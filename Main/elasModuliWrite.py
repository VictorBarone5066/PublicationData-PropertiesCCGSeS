import headerElas as h
from headerPoscar import Poscar

tetraFolder = "t-tElas"
orthoFolder = "o-oElas"

CIJSIJ_LOC = {"tetra": tetraFolder + "//" + "tetraCijsSijs.csv",
              "ortho": orthoFolder + "//" + "orthoCijsSijs.csv"}
CONTCAR_LOC = {"tetra": tetraFolder + "//" + "CONTCARx",
               "ortho": orthoFolder + "//" + "CONTCARx"}
MODULI_LOC = "elasModuli.csv"

OUTFILE_LOC = "allElasModuli.csv"

#Atomic masses in kg/mol
MASSES = {"Cu": 063.546E-3,
          "Cd": 112.410E-3,
          "Ge": 072.630E-3,
          "Se": 087.971E-3,
          "S":  032.060E-3}

#Gives total mass of a cell in kg
def GetCellMass(conc):
    return (2.*MASSES["Cd"] + 4.*MASSES["Cu"] + 2.*MASSES["Ge"] + \
           8.*(conc*MASSES["S"] + (1.-conc)*MASSES["Se"]))/6.02214179E23

#Gives total volume in m^3
def GetCellVolume(phaseId, conc):
    p = Poscar(CONTCAR_LOC[phaseId] + f'{conc:.3f}'[0] + f'{conc:.3f}'[2:]) ##turns decimal conc into string
    return p.volume*(1E-10)**(3.)

with open(OUTFILE_LOC, 'w') as outfile:
    for phaseId in ["tetra", "ortho"]:
        moduliList = [] #[x, bulk, shear, young, poisson, c11-c22, c13-c23, c44-c55]
        data = h.ParseCijSijFile(CIJSIJ_LOC[phaseId])
        for dat in data:
            c, s = dat.cij, dat.sij

            #shear moduli
            reussShear = 15.0/(4.0*(s[1][1]+s[2][2]+s[3][3]) + \
                               3.0*(s[4][4]+s[5][5]+s[6][6]) - \
                               4.0*(s[1][2]+s[1][3]+s[2][3]))
            voigtShear = 1.0/15.0*(c[1][1]+c[2][2]+c[3][3]-c[1][2]-c[1][3]-c[2][3]) + \
                         1.0/05.0*(c[4][4]+c[5][5]+c[6][6])

            #bulk moduli
            reussBulk = 1.0/((s[1][1]+s[2][2]+s[3][3]) + 2*(s[1][2]+s[1][3]+s[2][3]))
            voigtBulk = 1.0/9.0*(c[1][1]+c[2][2]+c[3][3]) + 2.0/9.0*(c[1][2]+c[1][3]+c[2][3])

            #Good estimates are averages
            shear = 1.0/2.0*(reussShear + voigtShear) #GPa
            bulk = 1.0/2.0*(reussBulk + voigtBulk) #GPa

            #Youngs mod, Poissons ratio
            young = (9.0*bulk*shear)/(3.0*bulk + shear) #GPa
            poisson = (3.0*bulk - 2.0*shear)/(2.0*(3.0*bulk + shear)) #unitless

            #Vickers Hardness.
            #Emperical formula from "Microscopic theory of hardness and design of novel superhard crystals"
            #(Tian et al) 2012.  If Hardness is over 5 GPa, this estimate is good.
            hard = 0.92*(shear/bulk)**(1.137) * shear**(0.708)

            #Debye Temp.  ref: Z T Y Liu et al 2014 J. Phys.: Condens. Matter 26 025404
            h_ = 6.62606896E-34 #plank const in Js
            kB = 1.3806504E-23 #boltzman const in J/K
            Na = 6.02214179E23 #avagadros number in 1/mol
            nAtoms = 16. #number of atoms per cell
            massDen = GetCellMass(dat.x)/GetCellVolume(phaseId, dat.x) #cell mass density in kg/m^3
            vTra = (shear*1E9/massDen)**(1./2.) #transverse speed in m/s
            vLon = (((3*bulk + 4*shear)*1E9)/(3*massDen))**(1./2.) #longitudional speed in m/s
            vMean = ((1./3.)*(2./vTra**(3) + 1/vLon**(3)))**(-1./3.) #mean sound speed in m/s
            tDeb = h_/kB * ((3*nAtoms)/(4*3.14159)*(massDen/GetCellMass(dat.x)))**(1./3.) * vMean

            #Elastic constant differences
            d1122 = c[1][1]-c[2][2]
            d1323 = c[1][3]-c[2][3]
            d4455 = c[4][4]-c[5][5]

            moduliList.append([dat.x, bulk, shear, young, poisson, d1122, d1323, d4455, hard, tDeb])

        #Write output to csv
        outfile.write(phaseId + "\n")
        outfile.write("x,B,G,Y,\\nu,C_11-C_22,C_13-C_23,C_44-C_55,H,tDeb\n")
        for mod in moduliList:
            line = ""
            for elem in mod:
                line += str(elem) + ','
            outfile.write(line + "\n")
    outfile.close()
